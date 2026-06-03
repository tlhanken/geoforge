#[cfg(test)]
mod tests {
    use crate::CosmologyContext;
    use geoforge_types::cosmology::CosmologyPreset;
    use geoforge_types::Seed;

    #[test]
    fn minimal_is_deterministic() {
        let a = CosmologyContext::new(Seed::new(42), CosmologyPreset::Minimal);
        let b = CosmologyContext::new(Seed::new(42), CosmologyPreset::Minimal);
        let ma = a.primary_marker();
        let mb = b.primary_marker();
        assert_eq!(ma.position_ly.x_ly, mb.position_ly.x_ly);
        assert_eq!(ma.phenotype.primary_class, mb.phenotype.primary_class);
    }

    #[test]
    fn solar_system_matches_marker_position() {
        let ctx = CosmologyContext::new(Seed::new(99), CosmologyPreset::Minimal);
        let m = ctx.primary_marker();
        let sys = ctx.solar_system(ctx.primary_ref());
        assert_eq!(m.position_ly.x_ly, sys.barycenter_ly.x_ly);
        assert_eq!(m.phenotype.multiplicity, sys.phenotype.multiplicity);
    }

    #[test]
    fn rich_loads_more_markers_than_minimal() {
        let min = CosmologyContext::new(Seed::new(1), CosmologyPreset::Minimal);
        let rich = CosmologyContext::new(Seed::new(1), CosmologyPreset::Rich);
        assert!(rich.markers_near_focus(0).len() >= min.markers_near_focus(0).len());
    }
}
