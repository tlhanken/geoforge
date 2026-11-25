use geoforge_cosmos::nebulae::{near, NebulaKind};
fn main() {
    let s = geoforge_cosmos::sectors::Sectors::new(&geoforge_cosmos::galaxy::Galaxy::milky_way(), geoforge_core::seed::Seed(5));
    let sun = [s.galaxy().reference_kpc * 1000.0, 0.0, 20.0];
    let t = std::time::Instant::now();
    let all = near(&s, sun, 1500.0);
    println!("{:?} for {}", t.elapsed(), all.len());
    for k in NebulaKind::ALL {
        let v: Vec<_> = all.iter().filter(|n| n.kind == k).collect();
        let within1 = v.iter().filter(|n| n.distance_pc(sun) < 1000.0).count();
        println!("{:?}: {} (within 1 kpc {})", k, v.len(), within1);
        for n in v.iter().take(4) {
            println!(
                "   {} d={:.0} pc r={:.2} pc age={:?} Myr m={:?} Q={:?} V={:?} ang={:.3}°",
                n.name,
                n.distance_pc(sun),
                n.radius_pc,
                n.age_myr.map(|a| (a * 1000.0).round() / 1000.0),
                n.mass_msun.map(|m| m.round()),
                n.ionising_per_s.map(|q| q.log10()),
                n.apparent_v(sun).map(|v| (v * 10.0).round() / 10.0),
                n.angular_radius_deg(sun)
            );
        }
    }
    let t = std::time::Instant::now();
    let far = near(&s, sun, 3000.0);
    println!("3 kpc: {:?} for {}", t.elapsed(), far.len());
    let bright: Vec<_> = far.iter().filter(|n| n.apparent_v(sun).is_some_and(|v| v < 7.0)).collect();
    println!("brighter than 7: {}", bright.len());
    for n in bright.iter().take(10) {
        println!("   {:?} {} V={:.1} ang={:.2}°", n.kind, n.name, n.apparent_v(sun).unwrap(), n.angular_radius_deg(sun));
    }
}
