//! Every kind of body runs end to end: each layer is filled by the implementation its regime
//! calls for, or skipped with a reason, and the same seed gives the same world.

use geoforge_core::grid::GridSpec;
use geoforge_core::product::{LayerId, ProductId};
use geoforge_core::{Archetype, Key, Seed};
use geoforge_pipeline::{self as pipe, GenerateOptions, Implementation as I, NewWorld};

fn world(name: &str, body: Option<geoforge_core::BodyContext>, focus: Option<&str>) -> geoforge_store::World {
    let dir = std::env::temp_dir().join(format!("geoforge-types-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    pipe::create_world(&dir, &NewWorld { seed: Seed(11), grid: GridSpec::Equirect { width: 128, height: 64 }, body, focus: focus.map(str::to_string) }).unwrap()
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
        (Archetype::Earth, vec![Some(Tectonics), Some(Geology), Some(Elevation), Some(Circulation), Some(Climate), Some(Erosion), Some(Hydrology), Some(Biomes)]),
        (Archetype::Mars, vec![Some(LidTectonics), Some(LidGeology), Some(LidElevation), Some(Circulation), Some(Climate), Some(Erosion), None, None]),
        (Archetype::Venus, vec![Some(LidTectonics), Some(LidGeology), Some(LidElevation), Some(Circulation), Some(Climate), Some(Erosion), None, None]),
        (Archetype::Moon, vec![Some(LidTectonics), Some(LidGeology), Some(LidElevation), None, None, Some(Erosion), None, None]),
        (Archetype::Europa, vec![None; 8]),
        (Archetype::GasGiant, vec![None; 8]),
    ];
    assert_eq!(LayerId::ALL.len(), 8);
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
    let err = pipe::create_world(&dir, &NewWorld { seed: Seed(11), grid: GridSpec::Equirect { width: 64, height: 32 }, body: None, focus: Some("Nowhere".into()) });
    assert!(err.is_err());
}
