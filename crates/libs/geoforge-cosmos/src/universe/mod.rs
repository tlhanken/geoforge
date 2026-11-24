//! A universe of galaxies around a home galaxy, made as it is looked at.
//!
//! The levels above one galaxy: halos of dark matter by mass, placed cube by cube and
//! clustered ([`halos`]), and the galaxies in each, a central and its satellites
//! ([`members`]). Every galaxy has a stable identity ([`GalaxyId`]) that seeds it, so any
//! observer sees the same galaxy and a query generates only the cubes it reaches.
//!
//! Each galaxy has a model to draw its stars from ([`UniverseGalaxy::model`]): its type's
//! preset scaled to its mass and size. A system's galaxy is the one its address names
//! ([`crate::resolve`]); a system made before addresses is placed in the galaxy like its
//! model nearest its universe's origin ([`Universe::home`]). The Solar System's sky has the
//! real Local Group and nearby bright galaxies instead ([`sky`]).
//!
//! Not yet: superclusters' filaments beyond the lognormal field.

pub mod halos;
pub mod members;
pub mod sky;

use crate::galaxy::presets::Morphology;
use crate::galaxy::Galaxy;
use geoforge_core::address::{Address, GalaxyRef};
use geoforge_core::body::GalacticSite;
use geoforge_core::galactic::{GalaxyId, UniverseSite};
use geoforge_core::seed::Seed;
use halos::{bin_range, DensityField, Halo, BINS};
pub use members::UniverseGalaxy;
use members::{central, of_halo, stellar_mass, virial_radius_mpc};

/// A universe from a seed.
#[derive(Debug, Clone)]
pub struct Universe {
    seed: Seed,
    field: DensityField,
}

/// Broad classes a home galaxy is matched by.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Class {
    Spiral,
    Lenticular,
    Elliptical,
    Dwarf,
}

impl Class {
    fn of(m: Morphology) -> Self {
        match m {
            Morphology::Elliptical => Self::Elliptical,
            Morphology::Lenticular => Self::Lenticular,
            Morphology::DwarfIrregular | Morphology::DwarfSpheroidal => Self::Dwarf,
            Morphology::Spiral { .. } | Morphology::Starburst => Self::Spiral,
        }
    }
}

impl Universe {
    #[must_use]
    pub fn new(seed: Seed) -> Self {
        Self { seed, field: DensityField::new(seed) }
    }

    /// The generated universe around a site's galaxy, or `None` for the real one.
    #[must_use]
    pub fn of_site(site: &GalacticSite, system_seed: Seed) -> Option<Self> {
        match site.universe {
            Some(UniverseSite::Real) => None,
            Some(UniverseSite::Seeded(s)) => Some(Self::new(s)),
            None => Some(Self::new(site.galaxy_seed.unwrap_or_else(|| system_seed.child(b"galaxy")).child(b"universe"))),
        }
    }

    #[must_use]
    pub fn seed(&self) -> Seed {
        self.seed
    }

    /// The galaxy `id`, if it exists.
    #[must_use]
    pub fn galaxy(&self, id: &GalaxyId) -> Option<UniverseGalaxy> {
        let halo = halos::in_cube(self.seed, &self.field, usize::from(id.bin), id.cube).into_iter().find(|h| h.n == id.n)?;
        of_halo(self.seed, &halo).into_iter().find(|g| g.id == *id)
    }

    /// Every galaxy within `radius` Mpc of `centre` bright enough to reach magnitude
    /// `v_limit` from there (dwarfs are only looked for nearby).
    #[must_use]
    pub fn galaxies_seen(&self, centre: [f64; 3], radius: f64, v_limit: f64) -> Vec<UniverseGalaxy> {
        let mut out = Vec::new();
        for b in 0..BINS {
            let (_, hi) = bin_range(b);
            // The brightest galaxy such a halo could hold: 3σ heavy and bursting.
            let l = stellar_mass(hi) * 10f64.powf(0.45) / 0.5;
            let reach_mpc = 10f64.powf((v_limit - (4.83 - 2.5 * l.log10()) + 5.0) / 5.0) / 1e6;
            let r = radius.min(reach_mpc) + virial_radius_mpc(hi);
            for h in halos::within(self.seed, &self.field, centre, r, std::iter::once(b)) {
                out.extend(of_halo(self.seed, &h).into_iter().filter(|g| halos::dist(g.position_mpc, centre) <= radius));
            }
        }
        out
    }

    /// The galaxy a system at `placement` is in: the galaxy its address names, or (for a
    /// system made before addresses) the one like its galaxy `like` ([`Self::home`]).
    #[must_use]
    pub fn home_at(&self, placement: Option<&Address>, like: &Galaxy) -> UniverseGalaxy {
        match placement.and_then(|a| a.galaxy) {
            Some(GalaxyRef::Generated(id)) => self.galaxy(&id).unwrap_or_else(|| self.home(like)),
            _ => self.home(like),
        }
    }

    /// The galaxy like `like` nearest the universe's origin: of its type within 0.25 dex of
    /// its stellar mass (failing any within 128 Mpc, of its class, then of its mass alone).
    /// A universe's default galaxy is the one like the Milky Way.
    #[must_use]
    pub fn home(&self, like: &Galaxy) -> UniverseGalaxy {
        let target = like.mass_msun().max(1e5).log10();
        let class = Class::of(like.morphology);
        let bins: Vec<usize> = (0..BINS)
            .filter(|&b| {
                let (lo, hi) = bin_range(b);
                (stellar_mass((lo * hi).sqrt()).log10() - target).abs() < 0.8
            })
            .collect();
        let near = |r: f64, alike: &dyn Fn(&UniverseGalaxy) -> bool| -> Option<UniverseGalaxy> {
            halos::within(self.seed, &self.field, [0.0; 3], r, bins.iter().copied())
                .iter()
                .map(|h: &Halo| central(self.seed, h))
                .filter(|g| (g.stellar_msun.log10() - target).abs() < 0.25 && alike(g))
                .min_by(|a, b| halos::dist(a.position_mpc, [0.0; 3]).total_cmp(&halos::dist(b.position_mpc, [0.0; 3])))
        };
        let same = |g: &UniverseGalaxy| g.morphology == like.morphology;
        let kin = |g: &UniverseGalaxy| Class::of(g.morphology) == class;
        [8.0, 16.0, 32.0, 64.0, 128.0].into_iter().find_map(|r| near(r, &same)).or_else(|| near(128.0, &kin)).or_else(|| near(128.0, &|_| true)).unwrap_or_else(|| {
            let halo = Halo { bin: 0, cube: [0, 0, 0], n: u32::MAX, mass_msun: 1e12, position_mpc: [0.0; 3] };
            central(self.seed, &halo)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_milky_way_finds_its_place() {
        let u = Universe::new(Seed(21));
        let home = u.home(&Galaxy::milky_way());
        assert_eq!(home.morphology, Galaxy::milky_way().morphology, "the same type");
        assert!((home.stellar_msun.log10() - Galaxy::milky_way().mass_msun().log10()).abs() < 0.25);
        assert_eq!(u.galaxy(&home.id).as_ref(), Some(&home), "found again by its identity");
        // Its satellites are there when looked for.
        let around = u.galaxies_seen(home.position_mpc, 0.3, 30.0);
        assert!(around.iter().any(|g| g.id.cube == home.id.cube && g.id.n == home.id.n && g.id.member > 0), "satellites");
    }
}
