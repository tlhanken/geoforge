//! A universe directory: its default body made from the seed, another star's system made in
//! the same universe sharing its galaxy (and the edits over its stars), the focus switched
//! between them, and nothing from another universe let in.

use geoforge_core::address::{Address, Level, StarRef};
use geoforge_core::galactic::StarEdits;
use geoforge_core::grid::GridSpec;
use geoforge_core::Seed;
use geoforge_cosmos::resolve;
use geoforge_store::{UniverseDir, World};

#[test]
fn bodies_of_one_universe_share_its_galaxy() {
    let dir = std::env::temp_dir().join(format!("geoforge-universe-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let grid = GridSpec::Equirect { width: 64, height: 32 };
    let mut u = UniverseDir::create(&dir, Address::universe(Seed(7))).unwrap();
    let first = geoforge_pipeline::create_at(&mut u, &Address::universe(Seed(7)), None, grid).unwrap();
    let a = first.manifest.address.unwrap();
    assert_eq!(a, resolve::default_address(Seed(7)));
    assert_eq!(World::open(&dir).unwrap().root, first.root, "the directory opens at its focus");

    // Another star of the same galaxy.
    let galaxy = resolve::galaxy(&a).unwrap();
    let field = galaxy.starfield(StarEdits::default());
    let here = match a.star {
        Some(StarRef::Star(id)) => field.star(&id).unwrap().position_pc,
        _ => unreachable!(),
    };
    let other = field.near(here, 15.0, 14.0).into_iter().map(|s| a.up_to(Level::Galaxy).star(StarRef::Star(s.id))).find(|b| resolve::system(b, &field).is_ok() && b.star != a.star).unwrap();
    let second = geoforge_pipeline::create_at(&mut u, &other, None, grid).unwrap();
    assert_eq!(second.manifest.address.unwrap().up_to(Level::System), other);
    assert_eq!(u.doc.focus, second.manifest.address, "the newest is the focus");
    assert_eq!(first.galaxy_dir(), second.galaxy_dir(), "one galaxy");
    assert_ne!(first.load_system().unwrap().unwrap().root.seed, second.load_system().unwrap().unwrap().root.seed);

    // An edit made from one body is the other's too.
    let mut edits = first.load_star_edits().unwrap();
    edits.remove.push(edits.next_added());
    first.save_star_edits(&edits).unwrap();
    assert_eq!(second.load_star_edits().unwrap(), edits);

    // Asking again opens what is made; another universe's address is refused.
    let again = geoforge_pipeline::create_at(&mut u, &a, None, grid).unwrap();
    assert_eq!(again.root, first.root);
    assert_eq!(u.doc.focus, Some(a));
    assert!(geoforge_pipeline::create_at(&mut u, &Address::universe(Seed(8)), None, grid).is_err());
    assert_eq!(u.worlds().unwrap().len(), 2);
    let _ = std::fs::remove_dir_all(&dir);
}
