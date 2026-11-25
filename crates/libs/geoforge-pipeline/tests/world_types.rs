//! Every kind of body runs end to end: each layer is filled by the implementation its regime
//! calls for, or skipped with a reason, and the same seed gives the same world.

use geoforge_core::grid::{Grid, GridSpec};
use geoforge_core::product::{LayerId, ProductId};
use geoforge_core::{Archetype, Key, Seed};
use geoforge_pipeline::{self as pipe, GenerateOptions, Implementation as I, NewWorld};

fn world(name: &str, body: Option<geoforge_core::BodyContext>, focus: Option<&str>) -> geoforge_store::World {
    let dir = std::env::temp_dir().join(format!("geoforge-types-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    pipe::create_world(
        &dir,
        &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(11), grid: GridSpec::Equirect { width: 128, height: 64 }, body, system: None, focus: focus.map(str::to_string) },
    )
    .unwrap()
}

/// The implementation filling each layer, `None` where the layer does not apply.
fn fillers(w: &geoforge_store::World) -> Vec<Option<I>> {
    pipe::plan(w).unwrap().layers.iter().map(|lp| lp.not_applicable.is_none().then_some(lp.implementation).flatten()).collect()
}

#[test]
fn regimes_pick_their_generators() {
    use I::*;
    let archetype = |a: Archetype| world(a.slug(), Some(a.context(Key(1), a.name())), None);
    let cases = [
        (Archetype::Earth, vec![Some(Tectonics), Some(Geology), Some(Elevation), Some(Circulation), Some(Climate), Some(Erosion), Some(Hydrology), Some(Biomes), Some(GiantDeposits), None]),
        (Archetype::Mars, vec![Some(LidTectonics), Some(LidGeology), Some(LidElevation), Some(Circulation), Some(Climate), Some(Erosion), None, None, None, None]),
        (Archetype::Venus, vec![Some(LidTectonics), Some(LidGeology), Some(LidElevation), Some(Circulation), Some(Climate), Some(Erosion), None, None, None, None]),
        (Archetype::Moon, vec![Some(LidTectonics), Some(LidGeology), Some(LidElevation), None, None, Some(Erosion), None, None, None, None]),
        (Archetype::Europa, vec![Some(IceShellTectonics), Some(IceShellGeology), Some(IceShellElevation), None, None, Some(Erosion), None, None, None, None]),
        (Archetype::Titan, vec![Some(IceShellTectonics), Some(IceShellGeology), Some(IceShellElevation), Some(Circulation), Some(Climate), Some(Erosion), Some(MethaneCycle), None, None, None]),
        (Archetype::GasGiant, [vec![None; 9], vec![Some(GiantBands)]].concat()),
        (Archetype::IceGiant, [vec![None; 9], vec![Some(GiantBands)]].concat()),
    ];
    assert_eq!(LayerId::ALL.len(), 10);
    for (a, want) in cases {
        assert_eq!(fillers(&archetype(a)), want, "{a}");
    }
}

#[test]
fn lid_worlds_generate_deterministically() {
    for a in [Archetype::Mars, Archetype::Venus] {
        let mut w1 = world(&format!("{}-1", a.slug()), Some(a.context(Key(1), a.name())), None);
        let mut w2 = world(&format!("{}-2", a.slug()), Some(a.context(Key(1), a.name())), None);
        let (p1, _) = pipe::generate(&mut w1, &GenerateOptions::default(), |_| {}).unwrap();
        let (p2, _) = pipe::generate(&mut w2, &GenerateOptions::default(), |_| {}).unwrap();
        for p in [ProductId::Provinces, ProductId::Lithology, ProductId::Elevation, ProductId::Temperature] {
            assert!(p1.get(p).is_some(), "{a}: {} missing", p.name());
            assert_eq!(p1.get(p), p2.get(p), "{a}: {} differs", p.name());
        }
        // A dry world gets no rain.
        assert!(p1.precipitation().unwrap().values.iter().all(|&v| v == 0.0), "{a}");
    }
}

#[test]
fn ice_shells_generate_deterministically() {
    let a = Archetype::Europa;
    let mut w1 = world("europa-1", Some(a.context(Key(1), a.name())), None);
    let mut w2 = world("europa-2", Some(a.context(Key(1), a.name())), None);
    let (p1, _) = pipe::generate(&mut w1, &GenerateOptions::default(), |_| {}).unwrap();
    let (p2, _) = pipe::generate(&mut w2, &GenerateOptions::default(), |_| {}).unwrap();
    for p in [ProductId::Provinces, ProductId::Lithology, ProductId::Elevation] {
        assert!(p1.get(p).is_some(), "{}", p.name());
        assert_eq!(p1.get(p), p2.get(p), "{} differs", p.name());
    }
}

/// Titan's methane gathers in lakes near the poles, about 1% of the ground (Titan: 1.1%,
/// all poleward of 55°), and the same seed gives the same lakes.
#[test]
fn titans_lakes_gather_at_the_poles() {
    let a = Archetype::Titan;
    let mut w1 = world("titan-1", Some(a.context(Key(1), a.name())), None);
    let mut w2 = world("titan-2", Some(a.context(Key(1), a.name())), None);
    let (p1, _) = pipe::generate(&mut w1, &GenerateOptions::default(), |_| {}).unwrap();
    let (p2, _) = pipe::generate(&mut w2, &GenerateOptions::default(), |_| {}).unwrap();
    assert_eq!(p1.lakes(), p2.lakes());
    let grid = GridSpec::Equirect { width: 128, height: 64 }.build().unwrap();
    let lakes = &p1.lakes().expect("lakes").owner.values;
    let area = |c: usize| grid.area_sr(c as geoforge_core::grid::Cell);
    let wet: Vec<usize> = (0..grid.len()).filter(|&c| lakes[c] != 0).collect();
    let share = wet.iter().map(|&c| area(c)).sum::<f64>() / (4.0 * std::f64::consts::PI);
    assert!((0.002..0.03).contains(&share), "lakes cover {share}");
    assert!(wet.iter().all(|&c| grid.center(c as geoforge_core::grid::Cell).to_lat_lon().lat_deg.abs() > 45.0), "lakes only near the poles");
    assert!(p1.rivers().is_some_and(|r| !r.edges.is_empty()), "rivers run into them");
}

#[test]
fn any_body_of_a_system_can_be_the_focus() {
    let g = geoforge_cosmos::generate_system(Seed(11));
    let moon = g.system.nodes.iter().find(|n| n.parent.is_some_and(|p| p != g.system.nodes[0].key)).map(|n| n.name.clone());
    let belt = g.system.nodes.iter().find(|n| n.name.ends_with("Outer Belt")).map(|n| format!("{} #1", n.name)).unwrap();
    for name in moon.iter().chain([&belt]) {
        let mut w = world(&name.replace(' ', "_").replace('#', ""), None, Some(name));
        assert_eq!(&w.load_body().unwrap().name, name);
        let (_, reports) = pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
        assert_eq!(reports.len(), LayerId::ALL.len());
    }
    let dir = std::env::temp_dir().join(format!("geoforge-types-nobody-{}", std::process::id()));
    let err = pipe::create_world(
        &dir,
        &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(11), grid: GridSpec::Equirect { width: 64, height: 32 }, body: None, system: None, focus: Some("Nowhere".into()) },
    );
    assert!(err.is_err());
}

/// A giant gets bands and nothing else, the same each time.
#[test]
fn giants_get_bands() {
    let a = Archetype::GasGiant;
    let mut w1 = world("giant-1", Some(a.context(Key(1), a.name())), None);
    let mut w2 = world("giant-2", Some(a.context(Key(1), a.name())), None);
    let (p1, _) = pipe::generate(&mut w1, &GenerateOptions::default(), |_| {}).unwrap();
    let (p2, _) = pipe::generate(&mut w2, &GenerateOptions::default(), |_| {}).unwrap();
    let bands = p1.cloud_bands().unwrap();
    assert_eq!(Some(bands), p2.cloud_bands());
    assert!(p1.zonal_wind().unwrap().values.iter().any(|&u| u > 50.0) && p1.elevation().is_none());
    let (lo, hi) = bands.values.iter().fold((1.0f32, 0.0f32), |(l, h), &v| (l.min(v), h.max(v)));
    assert!(lo < 0.3 && hi > 0.8, "dark belts and bright zones: {lo}..{hi}");
}

/// Giant deposits: the same each time, on land, each where its kind's prospectivity is.
#[test]
fn giant_deposits_follow_the_ground() {
    let a = Archetype::Earth;
    let mut w1 = world("deposits-1", Some(a.context(Key(1), a.name())), None);
    let mut w2 = world("deposits-2", Some(a.context(Key(1), a.name())), None);
    let (p1, _) = pipe::generate(&mut w1, &GenerateOptions::default(), |_| {}).unwrap();
    let (p2, _) = pipe::generate(&mut w2, &GenerateOptions::default(), |_| {}).unwrap();
    let deposits = p1.deposits().unwrap();
    assert_eq!(Some(deposits), p2.deposits());
    assert!(!deposits.rows.is_empty());
    let body = w1.load_body().unwrap();
    let grid = w1.manifest.grid.build().unwrap();
    let g = geoforge_analysis::ground::Ground::read(&p1, &grid, body.exterior.radius_km, true).unwrap();
    for d in deposits.rows.values() {
        let c = geoforge_core::Grid::locate(&grid, geoforge_core::UnitVec3::from_lat_lon(d.lat_deg, d.lon_deg)) as usize;
        let kind = geoforge_analysis::resources::Deposit::from_slug(&d.kind).unwrap();
        assert!(!g.sea[c], "{} at sea", d.kind);
        assert!(geoforge_analysis::resources::prospectivity(&g, kind)[c] > 0.0, "{} ({}) where none forms", d.kind, d.setting);
        assert!(d.ore_mt > 0.0 && d.depth_m >= 0.0);
    }
}

/// Flat slabs lift block ranges of continental basement far inland of a convergent margin,
/// never at the coast or on ocean floor, and some Earth-like seeds have them.
#[test]
fn flat_slabs_lift_block_ranges_inland() {
    use geoforge_core::surface::{BoundaryKind, CrustType, GeologicProvince};
    let mut found = 0;
    for seed in 1..=6 {
        let dir = std::env::temp_dir().join(format!("geoforge-flat-slab-{seed}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let earth = Archetype::Earth.context(Key(seed), "Earth");
        let mut w = pipe::create_world(
            &dir,
            &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(seed), grid: GridSpec::Equirect { width: 256, height: 128 }, body: Some(earth), system: None, focus: None },
        )
        .unwrap();
        let (p, _) = pipe::generate(&mut w, &GenerateOptions { until: Some(LayerId::Geology), ..GenerateOptions::default() }, |_| {}).unwrap();
        let grid = w.manifest.grid.build().unwrap();
        let segments = geoforge_surface::tectonics::boundaries::analyze(&grid, p.plates().unwrap(), 6371.0, &geoforge_surface::tectonics::TectonicsParams::default().boundaries);
        let (provinces, crust) = (p.provinces().unwrap(), p.crust_type().unwrap());
        for (key, cells) in provinces.cells_by_owner() {
            let prov = provinces.table.get(key).unwrap();
            if prov.kind != GeologicProvince::BasementUplift {
                continue;
            }
            found += 1;
            let seg = &segments[prov.boundary.expect("built from a margin")];
            assert_eq!(seg.kind, BoundaryKind::Convergent);
            let margin: Vec<_> = seg.cells.iter().map(|&c| grid.center(c)).collect();
            for c in cells {
                assert_eq!(*crust.get(c), CrustType::Continental, "seed {seed}");
                let km = margin.iter().map(|m| m.angle_to(grid.center(c)) * 6371.0).fold(f64::INFINITY, f64::min);
                assert!(km > 300.0, "seed {seed}: an uplift {km:.0} km from its margin");
            }
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
    assert!(found > 0, "no flat-slab uplifts on six Earth-like seeds");
}

/// The same seed makes the same world at any resolution: a coarse grid's land and heights
/// match a fine grid's averaged down to it.
#[test]
fn a_seed_is_the_same_world_at_any_resolution() {
    let make = |w: u32| {
        let dir = std::env::temp_dir().join(format!("geoforge-res-{w}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let mut world = pipe::create_world(
            &dir,
            &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(12345), grid: GridSpec::Equirect { width: w, height: w / 2 }, body: None, system: None, focus: None },
        )
        .unwrap();
        let (p, _) = pipe::generate(&mut world, &GenerateOptions { until: Some(LayerId::Elevation), ..GenerateOptions::default() }, |_| {}).unwrap();
        let e: Vec<f64> = p.rough_elevation().unwrap().values.iter().map(|&v| f64::from(v)).collect();
        let _ = std::fs::remove_dir_all(&dir);
        e
    };
    let (coarse, fine) = (make(128), make(512));
    // Average the fine grid's 4×4 blocks onto the coarse grid.
    let (w, k) = (128usize, 4usize);
    let down: Vec<f64> = (0..coarse.len())
        .map(|i| {
            let (x, y) = (i % w, i / w);
            (0..k).flat_map(|dy| (0..k).map(move |dx| (y * k + dy) * w * k + x * k + dx)).map(|j| fine[j]).sum::<f64>() / (k * k) as f64
        })
        .collect();
    let agree = coarse.iter().zip(&down).filter(|(a, b)| (**a >= 0.0) == (**b >= 0.0)).count() as f64 / coarse.len() as f64;
    let n = coarse.len() as f64;
    let (ma, mb) = (coarse.iter().sum::<f64>() / n, down.iter().sum::<f64>() / n);
    let cov: f64 = coarse.iter().zip(&down).map(|(a, b)| (a - ma) * (b - mb)).sum();
    let corr = cov / (coarse.iter().map(|a| (a - ma).powi(2)).sum::<f64>().sqrt() * down.iter().map(|b| (b - mb).powi(2)).sum::<f64>().sqrt());
    assert!(agree > 0.9, "land and sea agree on {:.0}% of cells", 100.0 * agree);
    assert!(corr > 0.85, "heights correlate {corr:.2}");
}

/// A pulsar's planet, lit only by its neutron star's glow (PSR B1257+12's are): every layer
/// either runs or is skipped with a reason, and the world comes out frozen and airless or
/// nearly so; the neutron star is a radio pulsar or quiet, and recycled by holding planets.
#[test]
fn a_pulsar_planet_runs_end_to_end() {
    let galaxy = geoforge_cosmos::galaxy::Galaxy::milky_way();
    let found = (0..5000u64).find_map(|s| {
        let mut ph = geoforge_cosmos::population::draw(&galaxy, Seed(s));
        (ph.primary_mass, ph.age_gyr) = (12.0, 3.0);
        ph.companions.clear();
        let g = geoforge_cosmos::galactic::build(Seed(s), &ph).ok()?;
        let ns = g.system.nodes.iter().find(|n| matches!(&n.node, geoforge_core::body::NodeKind::Star(st) if st.stage == geoforge_core::body::StellarStage::NeutronStar))?;
        Some((g.system.clone(), g.focus, ns.key))
    });
    let (system, focus, ns) = found.expect("a pulsar with planets within 5000 draws");
    let notes = geoforge_cosmos::remnants::notes(&system, ns);
    assert!(notes.iter().any(|n| n.starts_with("millisecond pulsar")), "{notes:?}");
    let body = geoforge_cosmos::any_body(&system, focus).unwrap();
    assert!(body.irradiance.flux_w_m2 < 1.0, "a pulsar barely lights it: {} W/m²", body.irradiance.flux_w_m2);
    let mut w = world("pulsar", Some(body), None);
    let (products, _) = pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
    assert!(products.elevation().is_some(), "it has a surface");
    let plan = pipe::plan(&w).unwrap();
    assert!(plan.layers.iter().all(|lp| lp.implementation.is_some() || lp.not_applicable.is_some()), "every layer runs or says why not");
    if let Some(t) = products.temperature() {
        let mean = t.values.iter().map(|&v| f64::from(v)).sum::<f64>() / t.values.len() as f64;
        assert!(mean < -150.0, "frozen: {mean} °C");
    }
}
