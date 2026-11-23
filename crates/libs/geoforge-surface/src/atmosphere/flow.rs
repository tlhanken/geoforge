//! Flows that conserve mass on the atmosphere grid: a wind or current field adjusted, as
//! little as it can be, so that what flows into each cell flows out again and nothing crosses
//! a closed edge (a coast for the sea).
//!
//! The field lives on a staggered (Arakawa C) grid: an eastward speed on each cell's east
//! face and a northward speed on its north face, with the cells' depth on the faces. The
//! correction is the gradient of a potential found by conjugate gradients: the
//! least-squares nearest field with the divergence asked for (Helmholtz projection).

use super::AtmoGrid;

/// The divergence a flow should keep.
#[derive(Debug, Clone, Copy)]
pub enum Keep {
    /// None: whatever flows in flows out (a gyre).
    Nothing,
    /// Its own, but none against a closed edge: what ran into the coast is turned along it.
    /// The rest of a surface current's divergence is real, water sinking or welling up
    /// (Ekman pumping, coastal and equatorial upwelling, the overturning's sinking).
    Interior,
    /// The divergence the same flow would have over a flat layer this deep, m: rising and
    /// sinking set by the large-scale pattern stay, but the layer thinning over high ground
    /// must be flowed around or squeezed through (air near the ground).
    FlatLayer(f64),
}

/// Geometry of the grid's faces, m.
struct Faces {
    /// Length of an east face (the same in every row) and of row `y`'s north face.
    east_len: f64,
    north_len: Vec<f64>,
    /// Distance between neighbours across an east face (by row) and a north face.
    east_gap: Vec<f64>,
    north_gap: f64,
}

impl Faces {
    fn of(g: &AtmoGrid, radius_km: f64) -> Self {
        let r = radius_km * 1000.0;
        let d = g.res.to_radians();
        let north_edge = |y: usize| (90.0 - y as f64 * g.res).to_radians();
        Self {
            east_len: r * d,
            north_len: (0..g.h).map(|y| r * d * north_edge(y).cos().max(0.0)).collect(),
            east_gap: (0..g.h).map(|y| r * d * g.lat(y).to_radians().cos().max(1e-3)).collect(),
            north_gap: r * d,
        }
    }
}

/// Neighbour across the east face of cell `i`, and across its north face (none in the top
/// row).
fn east_of(g: &AtmoGrid, i: usize) -> usize {
    (i / g.w) * g.w + (i % g.w + 1) % g.w
}
fn north_of(g: &AtmoGrid, i: usize) -> Option<usize> {
    (i >= g.w).then(|| i - g.w)
}

/// `flow` (m/s, east and north, per cell) made to conserve mass over the `open` cells of a
/// layer `depth` m deep (per cell): through every face between open cells the transport is
/// adjusted, through every face touching a closed cell it is zero.
#[must_use]
pub fn conserve(g: &AtmoGrid, radius_km: f64, open: &[bool], depth: &[f64], flow: &[(f64, f64)], keep: Keep) -> Vec<(f64, f64)> {
    let f = Faces::of(g, radius_km);
    let n = g.len();
    let face_depth = |a: usize, b: usize| if open[a] && open[b] { depth[a].min(depth[b]).max(0.0) } else { 0.0 };
    // Speeds on the faces, and their conductance (depth × length / gap). Stable air near the
    // ground does not climb over a range: through a face thinned by high ground the speed
    // falls with the layer's depth and resists change as much (it is blocked), so the solve
    // carries the rest round.
    let blocked = |a: usize, b: usize| match keep {
        Keep::Nothing | Keep::Interior => 1.0,
        Keep::FlatLayer(h) => (face_depth(a, b) / h).clamp(0.0, 1.0),
    };
    let mut ue = vec![0.0; n];
    let mut vn = vec![0.0; n];
    let mut ce = vec![0.0; n];
    let mut cn = vec![0.0; n];
    for i in 0..n {
        let y = i / g.w;
        let e = east_of(g, i);
        if open[i] && open[e] {
            ue[i] = 0.5 * (flow[i].0 + flow[e].0) * blocked(i, e);
            ce[i] = face_depth(i, e) * blocked(i, e) * f.east_len / f.east_gap[y];
        }
        if let Some(nb) = north_of(g, i) {
            if open[i] && open[nb] {
                vn[i] = 0.5 * (flow[i].1 + flow[nb].1) * blocked(i, nb);
                cn[i] = face_depth(i, nb) * blocked(i, nb) * f.north_len[y] / f.north_gap;
            }
        }
    }
    // Net outflow of each cell, by transport (depth × speed × face length).
    let outflow = |ue: &[f64], vn: &[f64], depth_of: &dyn Fn(usize, usize) -> f64| -> Vec<f64> {
        let mut out = vec![0.0; n];
        for i in 0..n {
            let y = i / g.w;
            let e = east_of(g, i);
            let q = depth_of(i, e) * ue[i] * f.east_len;
            out[i] += q;
            out[e] -= q;
            if let Some(nb) = north_of(g, i) {
                let q = depth_of(i, nb) * vn[i] * f.north_len[y];
                out[i] += q;
                out[nb] -= q;
            }
        }
        out
    };
    let have = outflow(&ue, &vn, &face_depth);
    let want = match keep {
        Keep::Nothing => vec![0.0; n],
        Keep::Interior => {
            // As if the coast let the flow through: each open cell's own speed across the face.
            let mut out = have.clone();
            for i in 0..n {
                let y = i / g.w;
                let e = east_of(g, i);
                if open[i] != open[e] {
                    let (o, q) = if open[i] { (i, depth[i] * flow[i].0) } else { (e, -depth[e] * flow[e].0) };
                    out[o] += q * f.east_len;
                }
                if let Some(nb) = north_of(g, i) {
                    if open[i] != open[nb] {
                        let (o, q) = if open[i] { (i, depth[i] * flow[i].1) } else { (nb, -depth[nb] * flow[nb].1) };
                        out[o] += q * f.north_len[y];
                    }
                }
            }
            out
        }
        Keep::FlatLayer(h) => {
            // The unblocked flow over flat ground.
            let flat = |a: usize, b: usize| if open[a] && open[b] { h } else { 0.0 };
            let ue0: Vec<f64> = (0..n).map(|i| 0.5 * (flow[i].0 + flow[east_of(g, i)].0)).collect();
            let vn0: Vec<f64> = (0..n).map(|i| north_of(g, i).map_or(0.0, |nb| 0.5 * (flow[i].1 + flow[nb].1))).collect();
            // Over high ground only the layer that is left rises or sinks with the pattern:
            // a plateau's thin layer cannot draw in or put out what a full column over the
            // plains does, which had blasted air off the plateaus.
            let mut out = outflow(&ue0, &vn0, &flat);
            for (o, d) in out.iter_mut().zip(depth) {
                *o *= (d / h).clamp(0.0, 1.0);
            }
            out
        }
    };
    // Solve L φ = have − want, (L φ)_i = Σ c (φ_i − φ_j) over open faces; moving c (φ_j − φ_i)
    // across each face from i to j then leaves every cell the outflow wanted.
    let mut rhs: Vec<f64> = (0..n).map(|i| have[i] - want[i]).collect();
    if matches!(keep, Keep::Interior | Keep::FlatLayer(_)) {
        // What ran into the coasts of a closed sea all told sinks evenly through it.
        balance(g, radius_km, &ce, &cn, open, &mut rhs);
    }
    let phi = solve(g, &ce, &cn, &rhs);
    for i in 0..n {
        let y = i / g.w;
        let e = east_of(g, i);
        if ce[i] > 0.0 {
            let dpth = face_depth(i, e);
            ue[i] += (phi[e] - phi[i]) * ce[i] / (dpth * f.east_len).max(1e-12);
        }
        if let Some(nb) = north_of(g, i) {
            if cn[i] > 0.0 {
                let dpth = face_depth(i, nb);
                vn[i] += (phi[nb] - phi[i]) * cn[i] / (dpth * f.north_len[y]).max(1e-12);
            }
        }
    }
    // Back to the cells: the mean of each cell's two faces in each direction.
    (0..n)
        .map(|i| {
            if !open[i] {
                return (0.0, 0.0);
            }
            let (x, y) = (i % g.w, i / g.w);
            let west = y * g.w + (x + g.w - 1) % g.w;
            let south = if y + 1 < g.h { vn[i + g.w] } else { 0.0 };
            (0.5 * (ue[i] + ue[west]), 0.5 * (vn[i] + south))
        })
        .collect()
}

/// Make `b` sum to zero over each region joined by open faces, spreading the excess by area.
fn balance(g: &AtmoGrid, radius_km: f64, ce: &[f64], cn: &[f64], open: &[bool], b: &mut [f64]) {
    let n = g.len();
    let r = radius_km * 1000.0;
    let d = g.res.to_radians();
    let area = |y: usize| r * r * d * ((90.0 - y as f64 * g.res).to_radians().sin() - (90.0 - (y + 1) as f64 * g.res).to_radians().sin()).abs();
    let mut region = vec![usize::MAX; n];
    for s in (0..n).filter(|&i| open[i]) {
        if region[s] != usize::MAX {
            continue;
        }
        let mut stack = vec![s];
        let mut cells = Vec::new();
        region[s] = s;
        while let Some(i) = stack.pop() {
            cells.push(i);
            let (x, y) = (i % g.w, i / g.w);
            let west = y * g.w + (x + g.w - 1) % g.w;
            let links = [(east_of(g, i), ce[i]), (west, ce[west]), (north_of(g, i).unwrap_or(i), cn[i]), (if y + 1 < g.h { i + g.w } else { i }, if y + 1 < g.h { cn[i + g.w] } else { 0.0 })];
            for (j, c) in links {
                if c > 0.0 && region[j] == usize::MAX {
                    region[j] = s;
                    stack.push(j);
                }
            }
        }
        let total: f64 = cells.iter().map(|&i| b[i]).sum();
        let a: f64 = cells.iter().map(|&i| area(i / g.w)).sum();
        for &i in &cells {
            b[i] -= total * area(i / g.w) / a;
        }
    }
}

/// Conjugate gradients on L φ = b, L the weighted graph Laplacian of the open faces
/// (conductances `ce`, `cn`), Jacobi-preconditioned. L is singular on each connected
/// region; `b` sums to zero on each, so the iteration converges within it.
fn solve(g: &AtmoGrid, ce: &[f64], cn: &[f64], b: &[f64]) -> Vec<f64> {
    let n = g.len();
    let apply = |x: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0; n];
        for i in 0..n {
            let e = east_of(g, i);
            if ce[i] > 0.0 {
                let d = ce[i] * (x[i] - x[e]);
                out[i] += d;
                out[e] -= d;
            }
            if let Some(nb) = north_of(g, i) {
                if cn[i] > 0.0 {
                    let d = cn[i] * (x[i] - x[nb]);
                    out[i] += d;
                    out[nb] -= d;
                }
            }
        }
        out
    };
    let mut diag = vec![0.0; n];
    for i in 0..n {
        let e = east_of(g, i);
        diag[i] += ce[i];
        diag[e] += ce[i];
        if let Some(nb) = north_of(g, i) {
            diag[i] += cn[i];
            diag[nb] += cn[i];
        }
    }
    let precond = |r: &[f64]| -> Vec<f64> { r.iter().zip(&diag).map(|(r, d)| if *d > 0.0 { r / d } else { 0.0 }).collect() };
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let mut x = vec![0.0; n];
    // Cells with no open face take no part.
    let mut r: Vec<f64> = b.iter().zip(&diag).map(|(b, d)| if *d > 0.0 { *b } else { 0.0 }).collect();
    let mut z = precond(&r);
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let tol = 1e-10 * dot(&r, &r).sqrt().max(1e-30);
    for _ in 0..4 * n.min(5000) {
        let ap = apply(&p);
        let pap = dot(&p, &ap);
        if pap <= 0.0 {
            break;
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        if dot(&r, &r).sqrt() <= tol {
            break;
        }
        z = precond(&r);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    x
}

/// Net outflow of each open cell per unit area, s⁻¹ × depth (m/s for a 1 m layer): zero
/// everywhere for a flow that conserves mass. Speeds on the faces are the mean of the two
/// cells', as [`conserve`] takes them.
#[must_use]
pub fn divergence(g: &AtmoGrid, radius_km: f64, open: &[bool], depth: &[f64], flow: &[(f64, f64)]) -> Vec<f64> {
    let f = Faces::of(g, radius_km);
    let n = g.len();
    let mut out = vec![0.0; n];
    let face_depth = |a: usize, b: usize| if open[a] && open[b] { depth[a].min(depth[b]).max(0.0) } else { 0.0 };
    for i in 0..n {
        let y = i / g.w;
        let e = east_of(g, i);
        let q = face_depth(i, e) * 0.5 * (flow[i].0 + flow[e].0) * f.east_len;
        out[i] += q;
        out[e] -= q;
        if let Some(nb) = north_of(g, i) {
            let q = face_depth(i, nb) * 0.5 * (flow[i].1 + flow[nb].1) * f.north_len[y];
            out[i] += q;
            out[nb] -= q;
        }
    }
    let r = radius_km * 1000.0;
    let d = g.res.to_radians();
    (0..n)
        .map(|i| {
            let y = i / g.w;
            let (top, bottom) = ((90.0 - y as f64 * g.res).to_radians(), (90.0 - (y + 1) as f64 * g.res).to_radians());
            let area = r * r * d * (top.sin() - bottom.sin()).abs();
            if open[i] {
                out[i] / area
            } else {
                0.0
            }
        })
        .collect()
}

/// The flow of a streamfunction `psi` (per cell, speed × km): u = −∂ψ/∂y, v = ∂ψ/∂x. Read
/// at the cells' corners, with ψ = 0 at every corner touching a closed cell, it is exactly
/// what flows into each cell flowing out, and no flow crosses a closed edge: where a gyre's
/// streamlines meet a shore they run along it.
#[must_use]
pub fn from_streamfunction(g: &AtmoGrid, radius_km: f64, open: &[bool], psi: &[f64]) -> Vec<(f64, f64)> {
    let (ue, vn) = streamfunction_faces(g, radius_km, open, psi);
    (0..g.len())
        .map(|i| {
            if !open[i] {
                return (0.0, 0.0);
            }
            let (x, y) = (i % g.w, i / g.w);
            let west = y * g.w + (x + g.w - 1) % g.w;
            let south = if y + 1 < g.h { vn[i + g.w] } else { 0.0 };
            (0.5 * (ue[i] + ue[west]), 0.5 * (vn[i] + south))
        })
        .collect()
}

/// Speeds through each cell's east face and north face of the flow of `psi`.
fn streamfunction_faces(g: &AtmoGrid, radius_km: f64, open: &[bool], psi: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let (w, h) = (g.w, g.h);
    let d = g.res.to_radians() * radius_km;
    // ψ at each cell's north-east corner; the poles' edges are closed.
    let ne: Vec<f64> = (0..g.len())
        .map(|i| {
            let (x, y) = (i % w, i / w);
            if y == 0 {
                return 0.0;
            }
            let xe = (x + 1) % w;
            let cells = [i, y * w + xe, i - w, (y - 1) * w + xe];
            if cells.iter().all(|&c| open[c]) {
                cells.iter().map(|&c| psi[c]).sum::<f64>() / 4.0
            } else {
                0.0
            }
        })
        .collect();
    let corner = |x: usize, y: usize| if y < h { ne[y * w + x % w] } else { 0.0 };
    let north_len = |y: usize| d * (90.0 - y as f64 * g.res).to_radians().cos().max(1e-3);
    let ue = (0..g.len()).map(|i| (corner(i % w, i / w + 1) - corner(i % w, i / w)) / d).collect();
    let vn = (0..g.len()).map(|i| (corner(i % w, i / w) - corner(i % w + w - 1, i / w)) / north_len(i / w)).collect();
    (ue, vn)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One closed basin (x 40–120 open) and a flow that runs into its walls and piles up.
    #[test]
    fn a_current_ends_with_nothing_piling_up_and_nothing_through_the_coast() {
        let g = AtmoGrid::new(2.0);
        let open: Vec<bool> = (0..g.len()).map(|i| (20..60).contains(&(i % g.w)) && (10..80).contains(&(i / g.w))).collect();
        let flow: Vec<(f64, f64)> = (0..g.len()).map(|i| (0.3 + 0.1 * ((i % 7) as f64), 0.2 * (g.lat(i / g.w).to_radians().sin()))).collect();
        let depth = vec![1.0; g.len()];
        let before = divergence(&g, 6371.0, &open, &depth, &flow);
        let fixed = conserve(&g, 6371.0, &open, &depth, &flow, Keep::Nothing);
        let after = divergence(&g, 6371.0, &open, &depth, &fixed);
        let worst = |d: &[f64]| d.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(worst(&before) > 1e-7, "the flow piles up before: {}", worst(&before));
        // Cell speeds are face means, so the check on them is loose; the faces are exact.
        assert!(worst(&after) < 0.05 * worst(&before), "and much less after: {} vs {}", worst(&after), worst(&before));
        // No flow into the walls: at the basin's east edge the eastward speed is half or less
        // of the interior's (its east face is shut).
        let y = 45;
        let (edge, inner) = (fixed[y * g.w + 59].0, fixed[y * g.w + 40].0);
        assert!(edge.abs() <= 0.6 * inner.abs().max(1e-6) || edge.abs() < 0.02, "edge {edge} interior {inner}");
        assert!(fixed.iter().zip(&open).all(|(f, o)| *o || *f == (0.0, 0.0)));
    }

    /// Air flowing east over flat ground keeps going; a high block in its way is flowed round
    /// (north and south of it the wind turns), while the large-scale divergence stays.
    #[test]
    fn wind_goes_round_high_ground() {
        let g = AtmoGrid::new(2.0);
        let open = vec![true; g.len()];
        let flow: Vec<(f64, f64)> = (0..g.len()).map(|_| (10.0, 0.0)).collect();
        let flat = conserve(&g, 6371.0, &open, &vec![3000.0; g.len()], &flow, Keep::FlatLayer(3000.0));
        assert!(flat.iter().all(|f| (f.0 - 10.0).abs() < 1e-6 && f.1.abs() < 1e-6), "flat ground changes nothing");
        // A 2.6 km block, 20° by 20°, at 30–50°N 80–100°E; the layer above it is 400 m deep.
        let block = |i: usize| (20..30).contains(&(i / g.w)) && (130..140).contains(&(i % g.w));
        let depth: Vec<f64> = (0..g.len()).map(|i| if block(i) { 400.0 } else { 3000.0 }).collect();
        let round = conserve(&g, 6371.0, &open, &depth, &flow, Keep::FlatLayer(3000.0));
        let at = |lat: f64, lon: f64| round[g.index(lat, lon)];
        // Upstream of the block the air splits: northward on its north half, southward on its south.
        assert!(at(48.0, 78.0).1 > 0.5, "turns north {:?}", at(48.0, 78.0));
        assert!(at(32.0, 78.0).1 < -0.5, "turns south {:?}", at(32.0, 78.0));
        // Over the block the wind is slow; round its ends it is fast.
        assert!(at(40.0, 90.0).0 < 5.0, "slow over the block {:?}", at(40.0, 90.0));
        assert!(at(54.0, 90.0).0 > 10.5, "faster past its end {:?}", at(54.0, 90.0));
    }

    /// A gyre drawn from a streamfunction: what flows into every cell flows out, and nothing
    /// flows through a coast, however the coast runs.
    #[test]
    fn a_streamfunction_gyre_closes_and_keeps_off_the_coasts() {
        let g = AtmoGrid { w: 90, h: 45, res: 4.0 };
        let r = 6371.0;
        // A basin from 20°W to 60°E with a cape and an island, between 50°S and 60°N.
        let open: Vec<bool> = (0..g.len())
            .map(|i| {
                let (lon, lat) = (-180.0 + (i % g.w) as f64 * 4.0 + 2.0, g.lat(i / g.w));
                let cape = (lon - 10.0).abs() < 6.0 && lat > 20.0 && lat < 34.0;
                let island = (lon - 40.0).abs() < 5.0 && (lat + 10.0).abs() < 5.0;
                lon > -20.0 && lon < 60.0 && lat > -50.0 && lat < 60.0 && !cape && !island
            })
            .collect();
        let psi: Vec<f64> = (0..g.len()).map(|i| (g.lat(i / g.w).to_radians() * 3.0).sin() * 1000.0 * (1.0 + (i % 7) as f64 / 7.0)).collect();
        let (ue, vn) = streamfunction_faces(&g, r, &open, &psi);
        let f = Faces::of(&g, r);
        let mut worst: f64 = 0.0;
        let mut scale: f64 = 0.0;
        for i in 0..g.len() {
            let (x, y) = (i % g.w, i / g.w);
            let e = east_of(&g, i);
            let west = y * g.w + (x + g.w - 1) % g.w;
            let south = if y + 1 < g.h { vn[i + g.w] * f.north_len[y + 1] } else { 0.0 };
            if open[i] != open[e] {
                assert!(ue[i].abs() < 1e-9, "flow through the coast east of {x},{y}: {}", ue[i]);
            }
            if let Some(nb) = north_of(&g, i) {
                if open[i] != open[nb] {
                    assert!(vn[i].abs() < 1e-9, "flow through the coast north of {x},{y}: {}", vn[i]);
                }
            }
            if open[i] {
                let out = (ue[i] - ue[west]) * f.east_len + vn[i] * f.north_len[y] - south;
                worst = worst.max(out.abs());
                scale = scale.max(ue[i].abs() * f.east_len);
            }
        }
        assert!(scale > 0.0 && worst < 1e-9 * scale, "net outflow {worst} against transports of {scale}");
    }
}
