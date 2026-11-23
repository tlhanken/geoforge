//! What gravity does to newly formed planets: spacing by mutual Hill radii, scattering of
//! pairs too close to survive, giants thrown so close to their star that tides make hot
//! Jupiters of them, belts bounded by a giant's resonances, and Trojans at its L4 and L5.
//!
//! A pair's mutual Hill radius is ((m₁+m₂)/3M)^(1/3)·(a₁+a₂)/2. Two planets are safe beyond
//! 2√3 ≈ 3.5 of them (Gladman 1993), but systems of several need ~8–10 to last billions of
//! years (Chambers et al. 1996; Jupiter and Saturn sit at ~8). Neighbours closer than 8 are
//! moved apart to 8; closer than 5, they scatter, as unstable systems do within millions of
//! years (Chatterjee et al. 2008). Eccentric neighbours also need room between the inner one's
//! farthest point and the outer one's nearest: a_out(1−e_out)/a_in(1+e_in) > 2.4·μ^(1/3)·
//! (a_out/a_in)^(1/2) + 1.15, μ the heavier planet's share of the star's mass (Petrovich 2015);
//! a pair that fails it scatters too. What a scattering does turns on the Safronov number
//! Θ = (m/M)·(a/R), how hard the heavier planet can kick against how hard the star holds: over
//! 1 (giants far out; Jupiter's is ~10) the lighter planet is thrown out and the survivor left
//! eccentric, by a Rayleigh law of scale 0.3 (Jurić & Tremaine 2008), as the eccentric giants
//! radial velocities find; under 1 (rocky planets close in; Earth's is 0.07) the two collide
//! and merge on a near-circular orbit, as compact systems do (Goldreich, Lithwick & Sari 2004;
//! Pu & Wu 2015).
//!
//! When two giants scatter, the survivor is sometimes driven (by the secular or Kozai–Lidov
//! pull of what is left) so eccentric that it passes within a few Roche radii of its star.
//! Tides there take its orbital energy but not its angular momentum, so it settles on a circle
//! at twice its nearest distance (Rasio & Ford 1996; Ford & Rasio 2006): a hot Jupiter of one to
//! five days. Simulations disagree on how often (Nagasawa, Ida & Bessho 2008 find ~30% with
//! three giants, Beaugé & Nesvorný 2012 a few percent); 7% of giant scatterings here, which
//! with disk migration makes ~1% of Sun-like stars hosts, as observed. The giant wrecks the
//! inner system on its way in: hot Jupiters are lonely (Steffen et al. 2012).

use crate::draw::Draw;
use crate::planets::{Belt, Formed};
use crate::system::planet_radius_earth;
use geoforge_core::body::{BodyClass, AU_KM, EARTH_MASS_KG, EARTH_RADIUS_KM, SOLAR_MASS_KG};

/// Neighbours are kept at least this many mutual Hill radii apart.
pub(crate) const HILL_SPACING: f64 = 8.0;
/// Closer than this, a pair scatters.
pub(crate) const HILL_SCATTER: f64 = 5.0;
/// Scale of the Rayleigh law of a scattering survivor's eccentricity.
const SCATTERED_E: f64 = 0.3;
/// Share of giant–giant scatterings that end in a hot Jupiter.
const HOT_JUPITER_SHARE: f64 = 0.07;
/// A giant driven inward reaches between one and this many Roche radii of its star.
const PERIASTRON_ROCHE: f64 = 3.0;

/// Petrovich's (2015) test for a pair of orbits, `inner` inside `outer`.
pub(crate) fn petrovich_stable(inner: &Formed, outer: &Formed, star_mass: f64) -> bool {
    let mu = inner.mass_earth.max(outer.mass_earth) * EARTH_MASS_KG / (star_mass * SOLAR_MASS_KG);
    outer.a_au * (1.0 - outer.e) / (inner.a_au * (1.0 + inner.e)) > 2.4 * mu.cbrt() * (outer.a_au / inner.a_au).sqrt() + 1.15
}

/// Safronov number of a planet around a star of `star_mass` M☉.
pub(crate) fn safronov(p: &Formed, star_mass: f64) -> f64 {
    let r_au = planet_radius_earth(p.class, p.mass_earth) * EARTH_RADIUS_KM / AU_KM;
    p.mass_earth * EARTH_MASS_KG / (star_mass * SOLAR_MASS_KG) * p.a_au / r_au
}

/// Roche limit of a planet against its star of `star_mass` M☉, AU (Jupiter's: 0.010 AU):
/// 2.16·R·(M/m)^(1/3) (Faber, Rasio & Willems 2005).
pub(crate) fn roche_au(p: &Formed, star_mass: f64) -> f64 {
    let r_au = planet_radius_earth(p.class, p.mass_earth) * EARTH_RADIUS_KM / AU_KM;
    2.16 * r_au * (star_mass * SOLAR_MASS_KG / (p.mass_earth * EARTH_MASS_KG)).cbrt()
}

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
/// planet at distance `fixed` never moves or scatters: its neighbours give way, and no giant
/// is sent through it toward the star. Returns how many planets were lost (thrown out or
/// merged).
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
        let crossing = !petrovich_stable(&inner, &outer, star_mass);
        if gap >= HILL_SPACING * (1.0 - 1e-9) && !crossing {
            i += 1;
            continue;
        }
        // Each side moves by this fraction of the pair's mean distance to reach the spacing.
        let spread = 0.5 * HILL_SPACING * hill_factor(&inner, &outer, star_mass);
        let pinned = is_fixed(&inner) || is_fixed(&outer);
        let inward = fixed.is_some_and(|f| outer.a_au <= f + 1e-12);
        if !pinned && (crossing || gap < HILL_SCATTER || spread >= 1.0) && safronov(if inner.mass_earth < outer.mass_earth { &outer } else { &inner }, star_mass) < 1.0 {
            // Too weak to throw each other out: they collide and merge.
            let (gone, kept) = if inner.mass_earth < outer.mass_earth { (i, i + 1) } else { (i + 1, i) };
            let lost = planets[gone].orbit().mass_earth;
            let merged = planets[kept].orbit_mut();
            merged.mass_earth += lost;
            merged.e = d.range(0.0, 0.1);
            planets.remove(gone);
            ejected += 1;
        } else if !pinned && (crossing || gap < HILL_SCATTER || spread >= 1.0) {
            // The lighter one is thrown out; the survivor's orbit turns eccentric.
            let (gone, kept) = if inner.mass_earth < outer.mass_earth { (i, i + 1) } else { (i + 1, i) };
            let giants = inner.class == BodyClass::GasGiant && outer.class == BodyClass::GasGiant;
            let e = d.rayleigh(SCATTERED_E).min(0.95);
            let hot = giants && fixed.is_none() && d.u() < HOT_JUPITER_SHARE;
            let survivor = planets[kept].orbit_mut();
            survivor.e = survivor.e.max(e);
            let from = survivor.a_au;
            if hot {
                // Driven to within a few Roche radii, then circularised at twice that.
                let q = roche_au(survivor, star_mass) * PERIASTRON_ROCHE.powf(d.u());
                (survivor.a_au, survivor.e) = (2.0 * q, d.range(0.0, 0.05));
            }
            planets.remove(gone);
            ejected += 1;
            if hot {
                // Everything it crossed on the way in is gone.
                let (before, kept) = (planets.len(), if gone < kept { kept - 1 } else { kept });
                let mut k = 0;
                planets.retain(|p| {
                    k += 1;
                    k - 1 == kept || p.orbit().a_au >= from
                });
                ejected += before - planets.len();
                by_distance(planets);
                i = 0;
                continue;
            }
        } else if crossing && pinned {
            // Too eccentric to share the system with the fixed planet: it is the one lost.
            planets.remove(if is_fixed(&inner) { i + 1 } else { i });
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

/// Least mass of a planet with Trojans, Earth masses.
const TROJAN_HOST_EARTH: f64 = 10.0;
/// A heavier giant this close (mutual Hill radii) leaves no stable Trojans.
const TROJAN_CLEARING_HILL: f64 = 15.0;

/// Trojan swarms, by index into `planets` (sorted by distance): at each giant of at least ten
/// Earth masses on a near-circular orbit (e < 0.1) with no heavier giant within 15 mutual Hill
/// radii, whose pull makes the Trojan region chaotic (Saturn's and Uranus's are nearly empty,
/// Nesvorný & Dones 2002; Jupiter's and Neptune's are full). A swarm spans the giant's orbit
/// ± a·√(8μ/3), the widest tadpole orbits (Jupiter's 4.94–5.46 AU).
pub(crate) fn trojans(planets: &[Formed], star_mass: f64) -> Vec<(usize, Belt)> {
    let giant = |p: &Formed| matches!(p.class, BodyClass::GasGiant | BodyClass::IceGiant) && p.mass_earth >= TROJAN_HOST_EARTH;
    planets
        .iter()
        .enumerate()
        .filter(|(_, p)| giant(p) && p.e < 0.1)
        .filter(|(i, p)| !planets.iter().enumerate().any(|(j, q)| j != *i && giant(q) && q.mass_earth > p.mass_earth && hill_separation(p, q, star_mass) < TROJAN_CLEARING_HILL))
        .map(|(i, p)| {
            let half = p.a_au * (8.0 * p.mass_earth * EARTH_MASS_KG / (3.0 * star_mass * SOLAR_MASS_KG)).sqrt();
            (i, (p.a_au - half, p.a_au + half))
        })
        .collect()
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
    fn eccentric_neighbours_need_room() {
        // Jupiter and Saturn as they are pass Petrovich's test; Saturn at e 0.5 does not.
        let (jupiter, saturn) = (Formed { e: 0.048, ..planet(5.20, 317.8) }, Formed { e: 0.054, ..planet(9.58, 95.2) });
        assert!(petrovich_stable(&jupiter, &saturn, 1.0));
        assert!(!petrovich_stable(&jupiter, &Formed { e: 0.5, ..saturn }, 1.0));
        let roche = roche_au(&jupiter, 1.0);
        assert!((0.009..0.012).contains(&roche), "Jupiter's Roche limit {roche} AU");
    }

    /// One giant scattering in fourteen makes a hot Jupiter of 1–5 days, which clears what it
    /// crossed; the rest leave an eccentric giant.
    #[test]
    fn scattered_giants_go_hot_or_eccentric() {
        let (mut hot, mut e, n) = (0, Vec::new(), 4000);
        for s in 0..n {
            let mut d = Draw::new(geoforge_core::Seed(s));
            let mut p = vec![planet(0.3, 3.0), planet(3.0, 318.0), planet(3.6, 150.0)];
            settle(&mut d, &mut p, 1.0, None);
            let giant = p.iter().find(|x| x.mass_earth > 300.0).unwrap();
            if giant.a_au < 0.1 {
                hot += 1;
                let days = giant.a_au.powf(1.5) * 365.25;
                assert!((0.7..6.0).contains(&days) && giant.e < 0.06, "{giant:?}: {days} d");
                assert_eq!(p.len(), 1, "the inner planet is gone: {p:?}");
            } else {
                e.push(giant.e);
            }
        }
        let share = f64::from(hot) / n as f64;
        assert!((0.04..0.1).contains(&share), "hot share {share}");
        let mean = e.iter().sum::<f64>() / e.len() as f64;
        assert!((0.3..0.45).contains(&mean), "survivors' mean e {mean} (Rayleigh 0.3: 0.38)");
    }

    #[test]
    fn close_rocky_planets_merge() {
        let (earth, jupiter) = (planet(1.0, 1.0), planet(5.2, 317.8));
        assert!(safronov(&earth, 1.0) < 0.1 && (5.0..15.0).contains(&safronov(&jupiter, 1.0)));
        let mut d = Draw::new(geoforge_core::Seed(3));
        let mut p = vec![planet(0.10, 2.0), planet(0.105, 1.0)];
        assert_eq!(settle(&mut d, &mut p, 1.0, None), 1);
        assert!((p[0].mass_earth - 3.0).abs() < 1e-12 && p[0].e <= 0.1, "{p:?}");
    }

    #[test]
    fn trojans_where_the_solar_system_has_them() {
        let ice = |a, m| Formed { class: BodyClass::IceGiant, ..planet(a, m) };
        let giants = [planet(5.20, 317.8), planet(9.58, 95.2), ice(19.2, 14.5), ice(30.1, 17.1)];
        let swarms = trojans(&giants, 1.0);
        let at: Vec<usize> = swarms.iter().map(|s| s.0).collect();
        assert_eq!(at, vec![0, 3], "Jupiter's and Neptune's, not Saturn's or Uranus's");
        let (lo, hi) = swarms[0].1;
        assert!((4.9..5.0).contains(&lo) && (5.4..5.5).contains(&hi), "{lo}–{hi}");
        assert!(trojans(&[Formed { e: 0.3, ..giants[0] }], 1.0).is_empty(), "none on an eccentric orbit");
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
