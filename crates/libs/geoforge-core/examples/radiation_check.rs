//! Balance temperature, optical depth and albedo of every preset.
use geoforge_core::radiation::{balance_temperature_c, open_albedo, Greenhouse};
use geoforge_core::{Archetype, Key};

fn main() {
    for a in Archetype::ALL {
        let b = a.context(Key(1), a.slug());
        let g = Greenhouse::of(&b);
        let t = balance_temperature_c(&b);
        println!("{:<10} T {:>7.1} °C  τ {:>8.3}  albedo {:.3}", a.slug(), t, g.optical_depth(t), open_albedo(&b.interior));
    }
}
