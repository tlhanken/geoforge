//! Constellations drawn from a sky's brightest stars, for a world that has drawn none.
//!
//! People join bright stars that lie close together into figures. Here the stars brighter
//! than 4.0 (some 500 in the Sun's sky) are joined by the shortest lines that connect them
//! all (a minimum spanning tree over each star's eight nearest), the tree is cut wherever a
//! line would be longer than 15°, and a figure of more than ten stars is split at its longest
//! lines until none is. Figures of three stars or more are kept, and pairs that hold a star
//! brighter than 1.5; a first-magnitude star left alone is joined to the two nearest stars to
//! the fifth magnitude, as Fomalhaut is to the Southern Fish. A figure of five stars or more
//! closes one loop, with its shortest line
//! not yet drawn, if that is no longer than its longest (Orion's belt and shoulders, the
//! Great Square). Each is named from its brightest star. The IAU's 88 figures, for
//! comparison, use stars to about the fifth magnitude and cover the whole sky.

use crate::sky::Neighbour;
use geoforge_core::galactic::{Constellation, StarId};
use geoforge_core::seed::Seed;

/// Faintest star used.
const FAINTEST: f64 = 4.0;
/// Nearest stars considered for each line.
const NEIGHBOURS: usize = 8;
/// Longest line, degrees.
const LONGEST_DEG: f64 = 15.0;
/// Most stars in a figure.
const MOST: usize = 10;
/// A pair is kept if it holds a star brighter than this.
const BRIGHT_PAIR: f64 = 1.5;
/// A lone bright star is joined to stars this bright.
const LONE_FAINTEST: f64 = 5.0;

/// The figures of `stars` (a sky's, any order), named from `seed`.
#[must_use]
pub fn automatic(stars: &[Neighbour], seed: Seed) -> Vec<Constellation> {
    let mut bright: Vec<&Neighbour> = stars.iter().filter(|s| s.v <= FAINTEST).collect();
    bright.sort_by(|a, b| a.v.total_cmp(&b.v));
    let n = bright.len();
    let angle = |i: usize, j: usize| dot(bright[i].direction, bright[j].direction).clamp(-1.0, 1.0).acos().to_degrees();

    // Candidate lines: each star to its nearest.
    let mut edges: Vec<(f64, usize, usize)> = Vec::new();
    for i in 0..n {
        let mut near: Vec<(f64, usize)> = (0..n).filter(|&j| j != i).map(|j| (angle(i, j), j)).collect();
        near.sort_by(|a, b| a.0.total_cmp(&b.0));
        edges.extend(near.into_iter().take(NEIGHBOURS).filter(|(a, _)| *a <= LONGEST_DEG).map(|(a, j)| (a, i.min(j), i.max(j))));
    }
    edges.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    edges.dedup_by(|a, b| a.1 == b.1 && a.2 == b.2);

    // The spanning forest (Kruskal).
    let mut parent: Vec<usize> = (0..n).collect();
    let mut tree = Vec::new();
    for &(a, i, j) in &edges {
        let (ri, rj) = (root(&mut parent, i), root(&mut parent, j));
        if ri != rj {
            parent[ri] = rj;
            tree.push((a, i, j));
        }
    }

    // Figures: components, split at their longest lines while too big.
    let mut work = vec![tree];
    let mut figures = Vec::new();
    while let Some(lines) = work.pop() {
        for part in components(&lines) {
            let stars = members(&part);
            if stars.len() > MOST {
                let mut part = part;
                part.sort_by(|a, b| a.0.total_cmp(&b.0));
                part.pop();
                work.push(part);
            } else {
                figures.push(part);
            }
        }
    }

    let mut out: Vec<(f64, Constellation)> = figures
        .into_iter()
        .filter_map(|mut lines| {
            let stars = members(&lines);
            let brightest = *stars.iter().min()?;
            if stars.len() < 3 && bright[brightest].v > BRIGHT_PAIR {
                return None;
            }
            if stars.len() >= 5 {
                let longest = lines.iter().map(|l| l.0).fold(0.0, f64::max);
                let closing = edges.iter().find(|(a, i, j)| *a <= longest && stars.contains(i) && stars.contains(j) && !lines.iter().any(|l| (l.1, l.2) == (*i, *j)));
                lines.extend(closing.copied());
            }
            let name = geoforge_names::name(bright[brightest].id.seed(seed).child(b"constellation").0);
            let lines: Vec<[StarId; 2]> = lines.iter().map(|&(_, i, j)| [bright[i].id, bright[j].id]).collect();
            Some((bright[brightest].v, Constellation { name, lines }))
        })
        .collect();
    // First-magnitude stars still alone: joined to their two nearest stars to the fifth.
    let used: Vec<StarId> = out.iter().flat_map(|(_, c)| c.lines.iter().flatten().copied()).collect();
    for lone in bright.iter().filter(|s| s.v < BRIGHT_PAIR && !used.contains(&s.id)) {
        let mut near: Vec<(f64, &Neighbour)> = stars
            .iter()
            .filter(|s| s.v <= LONE_FAINTEST && s.id != lone.id && !used.contains(&s.id))
            .map(|s| (dot(s.direction, lone.direction).clamp(-1.0, 1.0).acos().to_degrees(), s))
            .filter(|(a, _)| *a <= LONGEST_DEG)
            .collect();
        near.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut lines = Vec::new();
        let mut last = lone.id;
        for (_, s) in near.into_iter().take(2) {
            lines.push([last, s.id]);
            last = s.id;
        }
        if !lines.is_empty() {
            out.push((lone.v, Constellation { name: geoforge_names::name(lone.id.seed(seed).child(b"constellation").0), lines }));
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out.into_iter().map(|(_, c)| c).collect()
}

fn root(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

/// The stars of `lines`, sorted (by brightness, as their indices are).
fn members(lines: &[(f64, usize, usize)]) -> Vec<usize> {
    let mut s: Vec<usize> = lines.iter().flat_map(|l| [l.1, l.2]).collect();
    s.sort_unstable();
    s.dedup();
    s
}

/// `lines` split into connected parts.
fn components(lines: &[(f64, usize, usize)]) -> Vec<Vec<(f64, usize, usize)>> {
    let stars = members(lines);
    let index = |s: usize| stars.binary_search(&s).unwrap_or(0);
    let mut parent: Vec<usize> = (0..stars.len()).collect();
    for l in lines {
        let (a, b) = (root(&mut parent, index(l.1)), root(&mut parent, index(l.2)));
        parent[a] = b;
    }
    let mut groups: std::collections::BTreeMap<usize, Vec<(f64, usize, usize)>> = std::collections::BTreeMap::new();
    for &l in lines {
        let r = root(&mut parent, index(l.1));
        groups.entry(r).or_default().push(l);
    }
    groups.into_values().collect()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::galaxy::Galaxy;
    use crate::sky::{neighbourhood, Observer, NAKED_EYE};
    use crate::starfield::Starfield;
    use geoforge_core::galactic::StarEdits;

    #[test]
    fn a_sun_like_sky_gets_figures_like_ours() {
        let field = Starfield::new(&Galaxy::milky_way(), Seed(7), StarEdits::default());
        let o = Observer { position_pc: [8200.0, 0.0, 25.0], to_system: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], exclude: None, clear_pc: 0.5 };
        let sky = neighbourhood(&field, &o, NAKED_EYE);
        let figures = automatic(&sky.stars, Seed(3));
        assert!((40..120).contains(&figures.len()), "{} figures", figures.len());
        let dir = |id: &StarId| sky.stars.iter().find(|s| s.id == *id).unwrap().direction;
        for f in &figures {
            let mut stars: Vec<StarId> = f.lines.iter().flatten().copied().collect();
            stars.sort();
            stars.dedup();
            assert!((2..=MOST).contains(&stars.len()), "{}: {} stars", f.name, stars.len());
            assert!(f.lines.iter().all(|[a, b]| dot(dir(a), dir(b)).acos().to_degrees() <= LONGEST_DEG + 1e-9));
        }
        // Every first-magnitude star anchors a figure.
        let used: Vec<StarId> = figures.iter().flat_map(|f| f.lines.iter().flatten().copied()).collect();
        assert!(sky.stars.iter().filter(|s| s.v < BRIGHT_PAIR).all(|s| used.contains(&s.id)), "first-magnitude stars left out");
        // The same sky, the same figures.
        assert_eq!(automatic(&sky.stars, Seed(3)), figures);
    }
}
