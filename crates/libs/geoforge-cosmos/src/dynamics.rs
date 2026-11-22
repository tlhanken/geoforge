//! What gravity does to newly formed planets: spacing by mutual Hill radii, scattering of
//! pairs too close to survive, and belts bounded by a giant's resonances.
//!
//! A pair's mutual Hill radius is ((m₁+m₂)/3M)^(1/3)·(a₁+a₂)/2. Two planets are safe beyond
//! 2√3 ≈ 3.5 of them (Gladman 1993), but systems of several need ~8–10 to last billions of
//! years (Chambers et al. 1996; Jupiter and Saturn sit at ~8). Neighbours closer than 8 are
//! moved apart to 8; closer than 5, they scatter, as unstable systems do within millions of
//! years (Chatterjee et al. 2008): the lighter planet is thrown out and the survivor is left
//! on an eccentric orbit (e 0.1–0.5, as the eccentric giants Kepler and radial velocities
//! find). Eccentricities are not counted in the spacing.

use crate::draw::Draw;
use crate::planets::Formed;
use geoforge_core::body::EARTH_MASS_KG;
use geoforge_core::body::SOLAR_MASS_KG;

/// Neighbours are kept at least this many mutual Hill radii apart.
pub(crate) const HILL_SPACING: f64 = 8.0;
/// Closer than this, a pair scatters.
pub(crate) const HILL_SCATTER: f64 = 5.0;

/// Separation of two orbits in mutual Hill radii, around a star of `star_mass` M☉.
pub(crate) fn hill_separation(a: &Formed, b: &Formed, star_mass: f64) -> f64 {
    (b.a_au - a.a_au).abs() / (hill_factor(a, b, star_mass) * 0.5 * (a.a_au + b.a_au))
}

fn hill_factor(a: &Formed, b: &Formed, star_mass: f64) -> f64 {
    ((a.mass_earth + b.mass_earth) * EARTH_MASS_KG / (3.0 * star_mass * SOLAR_MASS_KG)).cbrt()
}

/// Something with an orbit to settle: a formed planet, or a draft that carries one.
pub(crate) trait HasOrbit {
    fn orbit(&self) -> &Formed;
    fn orbit_mut(&mut self) -> &mut Formed;
}

impl HasOrbit for Formed {
    fn orbit(&self) -> &Formed {
        self
    }
    fn orbit_mut(&mut self) -> &mut Formed {
        self
    }
}

/// Space `planets` (any order; returned by distance) around a star of `star_mass` M☉. The
/// planet at distance `fixed` never moves or scatters: its neighbours give way. Returns how
/// many planets were thrown out.
pub(crate) fn settle<T: HasOrbit>(d: &mut Draw, planets: &mut Vec<T>, star_mass: f64, fixed: Option<f64>) -> usize {
    let by_distance = |p: &mut Vec<T>| p.sort_by(|a, b| a.orbit().a_au.total_cmp(&b.orbit().a_au));
    let is_fixed = |p: &Formed| fixed.is_some_and(|f| (p.a_au - f).abs() < 1e-12);
    by_distance(planets);
    let mut ejected = 0;
    let mut i = 0;
    // Planets inside the fixed one give way inward, the rest outward, so each only ever moves
    // one way and the spacing settles.
    while i + 1 < planets.len() {
        let (inner, outer) = (*planets[i].orbit(), *planets[i + 1].orbit());
        let gap = hill_separation(&inner, &outer, star_mass);
        if gap >= HILL_SPACING * (1.0 - 1e-9) {
            i += 1;
            continue;
        }
        // Each side moves by this fraction of the pair's mean distance to reach the spacing.
        let spread = 0.5 * HILL_SPACING * hill_factor(&inner, &outer, star_mass);
        let pinned = is_fixed(&inner) || is_fixed(&outer);
        let inward = fixed.is_some_and(|f| outer.a_au <= f + 1e-12);
        if !pinned && (gap < HILL_SCATTER || spread >= 1.0) {
            // The lighter one is thrown out; the survivor's orbit turns eccentric.
            let (gone, kept) = if inner.mass_earth < outer.mass_earth { (i, i + 1) } else { (i + 1, i) };
            let e = &mut planets[kept].orbit_mut().e;
            *e = e.max(d.range(0.1, 0.5));
            planets.remove(gone);
            ejected += 1;
        } else if spread >= 1.0 {
            // No room at any distance: the one giving way is lost.
            planets.remove(if inward { i } else { i + 1 });
            ejected += 1;
        } else if inward {
            planets[i].orbit_mut().a_au = outer.a_au * (1.0 - spread) / (1.0 + spread);
        } else {
            planets[i + 1].orbit_mut().a_au = inner.a_au * (1.0 + spread) / (1.0 - spread);
        }
        // A moved planet may have passed another: re-sort and recheck from the pair before.
        by_distance(planets);
        i = i.saturating_sub(1);
    }
    ejected
}

/// Distance at which an orbit's period is `ratio` times a giant's at `a_giant`.
pub(crate) fn resonance_au(a_giant: f64, ratio: f64) -> f64 {
    a_giant * ratio.powf(2.0 / 3.0)
}

/// An asteroid belt inside a giant at `a_giant`: from its 4:1 resonance (or clear of the last
/// rocky planet at `last_rock`) to its 2:1, where the giant's kicks clear the orbits
/// (Jupiter's 2.06–3.28 AU). `None` if the rocky planets fill the gap.
pub(crate) fn asteroid_belt(a_giant: f64, last_rock: f64) -> Option<(f64, f64)> {
    let inner = resonance_au(a_giant, 0.25).max(last_rock * 1.3);
    let outer = resonance_au(a_giant, 0.5);
    (outer > inner * 1.2).then_some((inner, outer))
}

/// A Kuiper belt outside the outermost giant at `a_giant`: from its 3:2 resonance to its 2:1
/// (Neptune's 39.4–47.8 AU), where the classical belt ends.
pub(crate) fn kuiper_belt(a_giant: f64) -> (f64, f64) {
    (resonance_au(a_giant, 1.5), resonance_au(a_giant, 2.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::body::BodyClass;

    fn planet(a_au: f64, mass_earth: f64) -> Formed {
        Formed { a_au, e: 0.02, class: if mass_earth > 50.0 { BodyClass::GasGiant } else { BodyClass::Terrestrial }, mass_earth }
    }

    #[test]
    fn solar_system_numbers() {
        let (jupiter, saturn) = (planet(5.20, 317.8), planet(9.58, 95.2));
        let gap = hill_separation(&jupiter, &saturn, 1.0);
        assert!((7.0..9.5).contains(&gap), "Jupiter–Saturn {gap}");
        let (a, b) = asteroid_belt(5.2, 1.52).unwrap();
        assert!((a - 2.06).abs() < 0.02 && (b - 3.28).abs() < 0.02, "{a}–{b}");
        let (k1, k2) = kuiper_belt(30.1);
        assert!((k1 - 39.4).abs() < 0.2 && (k2 - 47.8).abs() < 0.2, "{k1}–{k2}");
    }

    #[test]
    fn crowded_giants_space_out_or_scatter() {
        let mut d = Draw::new(geoforge_core::Seed(1));
        // Two Jupiters at a distance ratio of 1.7: pushed apart. At 1.3: one is thrown out.
        let mut spaced = vec![planet(1.0, 318.0), planet(1.7, 318.0), planet(0.4, 1.0)];
        assert_eq!(settle(&mut d, &mut spaced, 1.0, None), 0);
        for w in spaced.windows(2) {
            assert!(hill_separation(&w[0], &w[1], 1.0) >= HILL_SPACING - 1e-9, "{spaced:?}");
        }
        let mut crowded = vec![planet(5.0, 318.0), planet(6.5, 100.0)];
        assert_eq!(settle(&mut d, &mut crowded, 1.0, None), 1);
        assert_eq!(crowded.len(), 1);
        assert!(crowded[0].mass_earth > 300.0 && crowded[0].e >= 0.1, "the heavier survives, eccentric: {crowded:?}");
    }

    #[test]
    fn the_fixed_planet_stays() {
        let mut d = Draw::new(geoforge_core::Seed(2));
        // An Earth with a Jupiter just outside and another just inside: they give way.
        let mut p = vec![planet(0.8, 318.0), planet(1.0, 1.0), planet(1.25, 318.0)];
        settle(&mut d, &mut p, 1.0, Some(1.0));
        let earth = p.iter().find(|x| x.mass_earth == 1.0).unwrap();
        assert!((earth.a_au - 1.0).abs() < 1e-12 && earth.e == 0.02);
        for w in p.windows(2) {
            assert!(hill_separation(&w[0], &w[1], 1.0) >= HILL_SPACING - 1e-9, "{p:?}");
        }
    }
}
