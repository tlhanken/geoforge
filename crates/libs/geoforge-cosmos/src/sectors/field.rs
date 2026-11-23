//! Field stars of one cube: a Poisson count of each bin's stars on the populations' mass in
//! the cube, each placed by the density inside it and drawn from its population's kinds.

use super::bins;
use super::{GalacticStar, Sectors};
use crate::draw::Draw;
use crate::galaxy::shape::draw_height;
use crate::galaxy::Cell;
use geoforge_core::galactic::{StarId, StarKind};
use geoforge_core::seed::Seed;

/// The draws of cube `cube` of bin `bin`.
pub(crate) fn cube_seed(base: Seed, label: &[u8], bin: usize, cube: [i32; 3]) -> Seed {
    cube.iter().fold(base.child(label).child_u64(b"bin", bin as u64), |s, &v| s.child_u64(b"axis", u64::from(v as u32)))
}

/// A point drawn from `cells` by mass: evenly across the cell, and for a disk's cells by
/// the exponential in height.
pub(crate) fn place(d: &mut Draw, cells: &[Cell], total: f64) -> [f64; 3] {
    let mut u = d.u() * total;
    let cell = cells
        .iter()
        .find(|c| {
            u -= c.mass;
            u <= 0.0
        })
        .unwrap_or(&cells[cells.len() - 1]);
    let x = cell.lo[0] + d.u() * (cell.hi[0] - cell.lo[0]);
    let y = cell.lo[1] + d.u() * (cell.hi[1] - cell.lo[1]);
    let z = match cell.height_pc {
        Some(h) => draw_height(d, cell.lo[2], cell.hi[2], h),
        None => cell.lo[2] + d.u() * (cell.hi[2] - cell.lo[2]),
    };
    [x, y, z]
}

impl Sectors {
    /// The field stars of `bins` (all at `level`) in cube `cube`.
    pub(crate) fn field_cube(&self, level: u32, cube: [i32; 3], bins: &[usize]) -> Vec<GalacticStar> {
        let (lo, hi) = bins::corners(cube, level);
        let pops: Vec<(usize, f64, &crate::galaxy::Prepared)> = self.field.iter().map(|p| (p.component, p.age, &p.shape)).collect();
        let cells: Vec<(Vec<Cell>, f64)> = self
            .cells(&pops, lo, hi)
            .into_iter()
            .zip(&self.field)
            .map(|(cells, pop)| {
                let mass = cells.iter().map(|c| c.mass).sum::<f64>();
                (cells, mass * pop.mass_msun)
            })
            .collect();
        let mut out = Vec::new();
        for &b in bins {
            let rates: Vec<f64> = self.field.iter().zip(&cells).map(|(pop, (_, m))| m * pop.table.per_msun[b]).collect();
            let total: f64 = rates.iter().sum();
            if total <= 0.0 {
                continue;
            }
            let mut d = Draw::new(cube_seed(self.seed(), b"field", b, cube));
            for n in 0..d.poisson(total) {
                let mut u = d.u() * total;
                let p = rates
                    .iter()
                    .position(|&r| {
                        u -= r;
                        u <= 0.0
                    })
                    .unwrap_or(rates.len() - 1);
                let pop = &self.field[p];
                let (pcells, pmass) = (&cells[p].0, cells[p].1 / pop.mass_msun);
                let position = place(&mut d, pcells, pmass);
                let Some(kind) = pop.table.bins[b].pick(d.u()) else { continue };
                let metallicity_dex = self.galaxy().metallicity(pop.component, position, kind.age_gyr, kind.scatter);
                let kind = StarKind { initial_mass: kind.initial_mass, age_gyr: kind.age_gyr, metallicity_dex, companion_mass: kind.companion };
                let id = StarId::Field { bin: b as u8, cube, n: n as u32 };
                out.push(GalacticStar::new(id, position, kind, Some(self.galaxy().components[pop.component].kind)));
            }
        }
        out
    }
}
