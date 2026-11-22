//! What a body is made of, from where it formed and what it receives: its shells, air,
//! water and tectonic regime.
//!
//! Rocky worlds (formed inside the frost line), by mass `m` (Earth = 1) and starlight `S`
//! (Earth = 1):
//!
//! | Setting | World | Like |
//! |---|---|---|
//! | `m < 0.05·√S` | airless rock, stagnant lid | Mercury, the Moon |
//! | `S > 1.1`, `m ≥ 0.3` | runaway greenhouse: ~90 bar CO₂, dry, episodic lid | Venus |
//! | `0.36 ≤ S ≤ 1.1`, `m ≥ 0.3` | ocean, ~1 bar N₂–O₂, plate tectonics | Earth |
//! | otherwise | thin CO₂, dry, stagnant lid | Mars |
//!
//! Icy worlds (formed beyond it) of planetary mass keep an ocean under an ice shell when
//! tides or their own heat keep it liquid (Europa), and the largest cold ones can hold a
//! nitrogen–methane atmosphere (Titan); the rest are frozen, cratered ice (Callisto).
//! Small bodies are inactive rubble of rock or ice.
//!
//! On Earth-like worlds the carbonate–silicate cycle sets CO₂: [`regulate_co2`].

use crate::draw::Draw;
use geoforge_core::body::{Archetype, BodyClass, BodyContext, BodyExterior, BodyInterior, Composition, Phase, Shell, ShellKind, Species, TectonicRegime, EARTH_MASS_KG};
use geoforge_core::derived::{MAXIMUM_GREENHOUSE_FLUX, RUNAWAY_GREENHOUSE_FLUX};
use geoforge_core::energy_balance::{Radiation, EARTH_HEAT_TRANSPORT, WEATHERING_E_FOLD_K};
use geoforge_core::seed::Seed;
use std::sync::OnceLock;

/// Where a body formed and what it receives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Setting {
    /// Mean starlight, W/m² (for a moon, at its planet's distance).
    pub flux_w_m2: f64,
    /// Formed beyond the frost line, so rich in ice.
    pub icy: bool,
    /// A moon close enough to a giant to be tidally heated.
    pub tidally_heated: bool,
}

/// Worlds placed exactly on a habitable-zone edge count as inside it.
const EDGE: f64 = 1e-9;

/// The interior of a body (its own layer).
#[must_use]
pub fn interior_for(seed: Seed, ext: &BodyExterior, setting: &Setting) -> BodyInterior {
    let mut d = Draw::new(seed.child(b"interior"));
    match ext.class {
        BodyClass::GasGiant => Archetype::GasGiant.interior(),
        BodyClass::IceGiant => Archetype::IceGiant.interior(),
        BodyClass::SmallBody => small_body(ext, setting),
        BodyClass::Terrestrial if setting.icy => icy_world(&mut d, ext, setting),
        BodyClass::Terrestrial => rocky_world(&mut d, ext, setting),
    }
}

fn rocky_world(d: &mut Draw, ext: &BodyExterior, setting: &Setting) -> BodyInterior {
    use Phase::{Liquid, Solid};
    use ShellKind::{Core, Crust, Mantle, Ocean};
    let m = ext.mass_kg / EARTH_MASS_KG;
    let s = setting.flux_w_m2 / 1361.0;
    let r = ext.radius_km;
    let core = r * d.range(0.5, 0.58);
    let crust = d.range(25.0, 45.0);
    let mut shells = vec![
        Shell::new(Core, Solid, core, Composition::from_weights(&[(Species::Iron, 0.85), (Species::Nickel, 0.05), (Species::Sulfur, 0.10)])),
        Shell::new(Mantle, Solid, r - core - crust, Composition::single(Species::Silicate)),
        Shell::new(Crust, Solid, crust, Composition::from_weights(&[(Species::Basalt, 0.6), (Species::Granite, 0.4)])),
    ];
    let holds_air = m >= 0.05 * s.sqrt();
    let big = m >= 0.3;
    let regime = if !holds_air {
        TectonicRegime::StagnantLid
    } else if big && s > RUNAWAY_GREENHOUSE_FLUX + EDGE {
        // Its oceans boiled away; the CO₂ they would have locked in rock stays in the air.
        let p = 90.0 * (m / 0.815) * d.range(0.5, 1.5);
        shells.push(Shell::atmosphere(p, 250.0, Composition::from_moles(&[(Species::CarbonDioxide, 0.965), (Species::Nitrogen, 0.035)])));
        TectonicRegime::Episodic
    } else if big && s >= MAXIMUM_GREENHOUSE_FLUX - EDGE {
        shells.push(Shell::new(Ocean, Liquid, d.range(2.35, 2.9), Composition::single(Species::Water)));
        let air = Composition::from_moles(&[(Species::Nitrogen, 0.78), (Species::Oxygen, 0.21), (Species::Argon, 0.0093), (Species::CarbonDioxide, 0.0004)]);
        shells.push(Shell::atmosphere(d.range(0.7, 1.6), 100.0, air));
        TectonicRegime::Plates
    } else {
        // Too small to keep its water, or too cold for it to stay liquid: a thin CO₂ remnant.
        let p = d.range(0.003, 0.02) * if big { 20.0 } else { 1.0 };
        shells.push(Shell::atmosphere(p, 11.0, Composition::from_moles(&[(Species::CarbonDioxide, 0.95), (Species::Nitrogen, 0.03), (Species::Argon, 0.02)])));
        TectonicRegime::StagnantLid
    };
    BodyInterior { shells, regime }
}

fn icy_world(d: &mut Draw, ext: &BodyExterior, setting: &Setting) -> BodyInterior {
    use Phase::{Liquid, Solid};
    use ShellKind::{Core, Crust, IceShell, Mantle, Ocean};
    let m = ext.mass_kg / EARTH_MASS_KG;
    let r = ext.radius_km;
    let rock = r * d.range(0.55, 0.75);
    let ice_total = r - rock;
    let ocean = m >= 0.005 && (setting.tidally_heated || (m >= 0.015 && d.u() < 0.5));
    let mut shells = vec![Shell::new(Core, Solid, rock, Composition::from_weights(&[(Species::Silicate, 0.8), (Species::Iron, 0.2)]))];
    let regime = if ocean {
        let lid = d.range(10.0, 150.0).min(0.4 * ice_total);
        let water = d.range(60.0, 150.0).min(0.5 * ice_total);
        shells.push(Shell::new(Mantle, Solid, ice_total - lid - water, Composition::single(Species::Ice)));
        shells.push(Shell::new(Ocean, Liquid, water, Composition::single(Species::Water)));
        shells.push(Shell::new(IceShell, Solid, lid, Composition::single(Species::Ice)));
        TectonicRegime::IceShell
    } else {
        // Frozen through: a dead, cratered crust of ice and rock.
        shells.push(Shell::new(Mantle, Solid, ice_total - 50.0_f64.min(0.3 * ice_total), Composition::single(Species::Ice)));
        shells.push(Shell::new(Crust, Solid, 50.0_f64.min(0.3 * ice_total), Composition::from_weights(&[(Species::Ice, 0.7), (Species::Silicate, 0.3)])));
        TectonicRegime::Inactive
    };
    // Titan: big and cold enough to keep nitrogen, with methane from its interior.
    if m >= 0.02 && setting.flux_w_m2 < 30.0 && d.u() < 0.4 {
        shells.push(Shell::atmosphere(d.range(0.8, 2.0), 600.0, Composition::from_moles(&[(Species::Nitrogen, 0.95), (Species::Methane, 0.05)])));
    }
    BodyInterior { shells, regime }
}

fn small_body(ext: &BodyExterior, setting: &Setting) -> BodyInterior {
    let crust = if setting.icy { Composition::from_weights(&[(Species::Ice, 0.6), (Species::Silicate, 0.4)]) } else { Composition::from_weights(&[(Species::Silicate, 0.85), (Species::Iron, 0.15)]) };
    BodyInterior {
        shells: vec![Shell::new(ShellKind::Mantle, Phase::Solid, ext.radius_km * 0.8, Composition::single(Species::Silicate)), Shell::new(ShellKind::Crust, Phase::Solid, ext.radius_km * 0.2, crust)],
        regime: TectonicRegime::Inactive,
    }
}

/// Weathering also speeds up with CO₂ itself, as its partial pressure to this power (Walker,
/// Hays & Kasting 1981).
pub const WEATHERING_CO2_EXPONENT: f64 = 0.3;

/// Earth's CO₂, bar: the reference the thermostat balances against.
const EARTH_CO2_BAR: f64 = 4.2e-4;

/// The carbonate–silicate thermostat: on a world with oceans and plate tectonics, weathering
/// draws CO₂ down and volcanoes build it up until the two balance. With Earth's outgassing,
/// CO₂ settles where the world weathers as fast as Earth (Walker, Hays & Kasting 1981):
/// `(pCO₂/pCO₂⊕)^0.3 · e^((Tw − Tw⊕)/13.7 K) = 1`, with `Tw` the weathering temperature of
/// the latitude-band energy balance ([`Radiation::weathering_temperature_c`]). Weathering
/// runs on warm, unfrozen ground, so a dim world builds up CO₂ until it thaws, and settles
/// somewhat colder than Earth (the extra CO₂ weathers faster too); a bright one draws CO₂
/// down and runs warmer. At the ice-albedo tipping point weathering jumps; the thermostat
/// settles on the warm side, as Earth's climbed out of its snowballs.
///
/// So the way a world turns matters. A tidally locked world weathers only on its day side
/// (Edson et al. 2012): with its night half frozen it builds up CO₂ before the day side
/// weathers as much as Earth does. A fast spinner holds its heat near the equator: about the
/// same CO₂, but colder poles and a colder mean; a slow one evens out. Assumes land
/// everywhere alike: a continent under the star would draw CO₂ down further.
/// Sets the atmosphere's CO₂; other worlds are left alone.
pub fn regulate_co2(body: &mut BodyContext) {
    let i = &body.interior;
    if i.regime != TectonicRegime::Plates || !i.has_surface_liquid_water() || !i.has_atmosphere() {
        return;
    }
    let locked = body.is_tidally_locked();
    static EARTH: OnceLock<f64> = OnceLock::new();
    let earth = *EARTH.get_or_init(|| {
        let e = Archetype::Earth.context(body.key, "Earth");
        Radiation::for_body(&e, EARTH_HEAT_TRANSPORT, 0.0).weathering_temperature_c(false)
    });
    // Weathering relative to Earth's, in e-folds; it rises with CO₂, so bisect on CO₂'s
    // logarithm between 1 ppm and 1 bar.
    let excess = |b: &BodyContext, ln_co2: f64| {
        let tw = Radiation::for_body(b, EARTH_HEAT_TRANSPORT, 0.0).weathering_temperature_c(locked);
        (tw - earth) / WEATHERING_E_FOLD_K + WEATHERING_CO2_EXPONENT * (ln_co2 - EARTH_CO2_BAR.ln())
    };
    let (mut lo, mut hi) = ((1e-6_f64).ln(), 0.0_f64);
    while hi - lo > 2e-3 {
        let mid = 0.5 * (lo + hi);
        set_co2(body, mid.exp());
        if excess(body, mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    // The upper end weathers at least as fast as Earth: across a jump, the warm side.
    set_co2(body, hi.exp());
}

/// Put `p_bar` of CO₂ in the air, keeping the other gases' proportions and total pressure
/// of the rest.
fn set_co2(body: &mut BodyContext, p_bar: f64) {
    let Some(air) = body.interior.shells.iter_mut().rev().find(|s| s.kind == ShellKind::Atmosphere) else { return };
    let total = air.pressure_bar.unwrap_or(0.0);
    let others: Vec<(Species, f64)> = air.composition.parts().iter().filter(|(s, _)| *s != Species::CarbonDioxide).map(|&(s, _)| (s, air.composition.mole_fraction(s))).collect();
    let rest: f64 = others.iter().map(|o| o.1).sum();
    let p_rest = total * rest;
    let new_total = p_rest + p_bar;
    let mut moles: Vec<(Species, f64)> = others.iter().map(|&(s, x)| (s, x / rest * p_rest / new_total)).collect();
    moles.push((Species::CarbonDioxide, p_bar / new_total));
    air.composition = Composition::from_moles(&moles);
    air.pressure_bar = Some(new_total);
}

/// One line on what a body is: regime, air and water.
#[must_use]
pub fn describe(body: &BodyContext) -> String {
    let i = &body.interior;
    if !i.has_solid_surface() {
        return "no solid surface".into();
    }
    let regime = match i.regime {
        TectonicRegime::Plates => "plate tectonics",
        TectonicRegime::StagnantLid => "stagnant lid",
        TectonicRegime::Episodic => "episodic lid",
        TectonicRegime::IceShell => "ice shell over an ocean",
        TectonicRegime::Inactive => "inactive",
    };
    let air = if i.has_atmosphere() {
        let main = i.shell(ShellKind::Atmosphere).and_then(|a| a.composition.parts().iter().map(|&(s, _)| (s, a.composition.mole_fraction(s))).max_by(|x, y| x.1.total_cmp(&y.1)));
        format!("{:.3} bar {}", i.surface_pressure_bar(), main.map_or("?", |(s, _)| s.formula()))
    } else {
        "airless".into()
    };
    let water = if i.has_surface_liquid_water() {
        format!(", ocean {:.1} km, CO₂ {:.0} ppm", i.surface_ocean_gel_km(), 1e6 * i.partial_pressure_bar(Species::CarbonDioxide) / i.surface_pressure_bar())
    } else {
        String::new()
    };
    format!("{regime}, {air}{water}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::body::EARTH_RADIUS_KM;
    use geoforge_core::Key;

    fn ext(m: f64) -> BodyExterior {
        BodyExterior { class: BodyClass::Terrestrial, mass_kg: m * EARTH_MASS_KG, radius_km: m.powf(0.27) * EARTH_RADIUS_KM, rotation_h: 24.0, obliquity_deg: 20.0 }
    }

    fn rocky(m: f64, s: f64) -> BodyInterior {
        interior_for(Seed(1), &ext(m), &Setting { flux_w_m2: s * 1361.0, icy: false, tidally_heated: false })
    }

    #[test]
    fn rocky_worlds_by_mass_and_light() {
        let earth = rocky(1.0, 1.0);
        assert_eq!(earth.regime, TectonicRegime::Plates);
        assert!(earth.has_surface_liquid_water());
        let venus = rocky(0.8, 1.9);
        assert_eq!(venus.regime, TectonicRegime::Episodic);
        assert!(venus.surface_pressure_bar() > 30.0 && !venus.has_surface_liquid_water());
        let mars = rocky(0.1, 0.43);
        assert_eq!(mars.regime, TectonicRegime::StagnantLid);
        assert!(mars.has_atmosphere() && mars.surface_pressure_bar() < 0.05);
        let mercury = rocky(0.055, 6.7);
        assert!(!mercury.has_atmosphere());
        let frozen = rocky(1.0, 0.25);
        assert!(!frozen.has_surface_liquid_water() && frozen.regime == TectonicRegime::StagnantLid);
    }

    #[test]
    fn icy_worlds() {
        let icy = |m: f64, heated: bool, seed: u64| interior_for(Seed(seed), &ext(m), &Setting { flux_w_m2: 15.0, icy: true, tidally_heated: heated });
        let europa = icy(0.008, true, 1);
        assert_eq!(europa.regime, TectonicRegime::IceShell);
        assert!(europa.has_solid_surface() && !europa.has_lithosphere());
        let callisto = icy(0.003, false, 1);
        assert_eq!(callisto.regime, TectonicRegime::Inactive);
        assert!(callisto.has_lithosphere(), "frozen crust: cratered by the lid generators");
        // Some big cold moons keep a Titan-like atmosphere.
        assert!((0..40).any(|s| icy(0.0225, false, s).has_atmosphere()));
    }

    #[test]
    fn thermostat_by_spin_and_lock() {
        let world = |rot: Option<f64>, flux: f64| {
            let mut b = Archetype::Earth.context(Key(1), "x");
            b.irradiance.flux_w_m2 = flux * 1361.0;
            b.exterior.rotation_h = rot.unwrap_or(b.irradiance.year_days * 24.0);
            regulate_co2(&mut b);
            let lats: Vec<f64> = (0..36).map(|i| -87.5 + 5.0 * f64::from(i)).collect();
            let r = Radiation::for_body(&b, EARTH_HEAT_TRANSPORT, 0.0);
            let (t, mean) = if b.is_tidally_locked() { r.substellar_temperature(&lats) } else { r.zonal_temperature(&lats) };
            let ppm = 1e6 * b.interior.partial_pressure_bar(Species::CarbonDioxide) / b.interior.surface_pressure_bar();
            (ppm, mean, t.iter().copied().fold(f64::MIN, f64::max), t.iter().copied().fold(f64::MAX, f64::min))
        };
        let earth = world(Some(23.934), 1.0);
        let fast = world(Some(6.0), 1.0);
        let slow = world(Some(48.0), 1.0);
        assert!(fast.1 < earth.1 - 4.0 && slow.1 > earth.1 && fast.3 < earth.3 - 15.0, "spin: {fast:?} {earth:?} {slow:?}");
        assert!((fast.0 / earth.0 - 1.0).abs() < 0.4, "about the same CO₂: {} vs {}", fast.0, earth.0);
        let locked = world(None, 1.0);
        assert!(locked.0 > 5.0 * earth.0 && locked.2 > 30.0 && locked.3 < -20.0 && locked.1 < earth.1, "locked: {locked:?}");
    }

    /// Dim worlds thaw: on the fine bands the climate layer uses, the regulated world is not
    /// left frozen over, even near the ice-albedo tipping point.
    #[test]
    fn thermostat_thaws_dim_worlds() {
        let lats: Vec<f64> = (0..180).map(|i| -89.5 + f64::from(i)).collect();
        for flux in [0.8, 0.85, 0.88, 0.9, 0.95, 1.05, 1.1] {
            let mut b = Archetype::Earth.context(Key(1), "x");
            b.irradiance.flux_w_m2 = flux * 1361.0;
            regulate_co2(&mut b);
            let (t, mean) = Radiation::for_body(&b, EARTH_HEAT_TRANSPORT, 0.0).zonal_temperature(&lats);
            let ppm = 1e6 * b.interior.partial_pressure_bar(Species::CarbonDioxide) / b.interior.surface_pressure_bar();
            assert!(mean > 0.0 && t[90] > 15.0, "{flux}× Earth's light: mean {mean:.1} °C, equator {:.1}, {ppm:.0} ppm", t[90]);
        }
    }

    #[test]
    fn thermostat_compensates_starlight() {
        let mut dim = Archetype::Earth.context(Key(1), "Dim");
        dim.irradiance.flux_w_m2 = 0.9 * 1361.0;
        let w = |b: &BodyContext| Radiation::for_body(b, EARTH_HEAT_TRANSPORT, 0.0).weathering_temperature_c(b.is_tidally_locked());
        let cold = w(&dim);
        regulate_co2(&mut dim);
        let warmed = w(&dim);
        let earth = w(&Archetype::Earth.context(Key(1), "Earth"));
        assert!(warmed > cold + 5.0, "{cold} → {warmed}");
        let co2 = dim.interior.partial_pressure_bar(Species::CarbonDioxide);
        assert!(co2 > 4.2e-4 * 1.5 && warmed < earth, "{co2} bar, {warmed} vs Earth's {earth}");
        let balance = (warmed - earth) / WEATHERING_E_FOLD_K + WEATHERING_CO2_EXPONENT * (co2 / 4.2e-4).ln();
        assert!(balance.abs() < 0.02, "weathers as fast as Earth: {balance}");
        // Earth itself needs no change.
        let mut earth_body = Archetype::Earth.context(Key(1), "Earth");
        regulate_co2(&mut earth_body);
        let ppm = 1e6 * earth_body.interior.partial_pressure_bar(Species::CarbonDioxide) / earth_body.interior.surface_pressure_bar();
        assert!((ppm / 420.0 - 1.0).abs() < 0.05, "{ppm} ppm");
        // Not on a dry world.
        let mut venus = Archetype::Venus.context(Key(1), "Venus");
        let before = venus.clone();
        regulate_co2(&mut venus);
        assert_eq!(venus, before);
    }
}
