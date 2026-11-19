//! Zonal annual temperatures and seasonal ranges of every preset with an atmosphere.
use geoforge_core::{Archetype, Key};
use geoforge_surface::atmosphere::Radiation;

fn main() {
    let lats: Vec<f64> = (0..18).map(|i| 85.0 - 10.0 * f64::from(i)).collect();
    for a in Archetype::ALL {
        let b = a.context(Key(1), a.slug());
        if !b.interior.has_solid_surface() || !b.interior.has_atmosphere() {
            continue;
        }
        let r = Radiation::for_body(&b, 3.2, 0.0);
        let (t, mean) = r.zonal_temperature(&lats);
        let row: Vec<String> = t.iter().map(|v| format!("{v:.0}")).collect();
        println!("{:<6} mean {mean:>6.1}  85N..85S [{}]  range@45N {:.1}", a.slug(), row.join(" "), r.seasonal_range(45.0));
    }
}
