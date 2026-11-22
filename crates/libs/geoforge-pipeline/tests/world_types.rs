//! Every kind of body runs end to end: each layer is filled by the implementation its regime
//! calls for, or skipped with a reason, and the same seed gives the same world.

use geoforge_core::grid::GridSpec;
use geoforge_core::product::{LayerId, ProductId};
use geoforge_core::{Archetype, Key, Seed};
use geoforge_pipeline::{self as pipe, GenerateOptions, Implementation as I, NewWorld};

fn world(name: &str, body: Option<geoforge_core::BodyContext>, focus: Option<&str>) -> geoforge_store::World {
    let dir = std::env::temp_dir().join(format!("geoforge-types-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    pipe::create_world(&dir, &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(11), grid: GridSpec::Equirect { width: 128, height: 64 }, body, focus: focus.map(str::to_string) })
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
        (Archetype::Europa, vec![None; 10]),
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
        &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(11), grid: GridSpec::Equirect { width: 64, height: 32 }, body: None, focus: Some("Nowhere".into()) },
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
