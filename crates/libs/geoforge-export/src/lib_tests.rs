#[cfg(test)]
mod tests {
    use crate::run_pre_tectonic_layers;
    use crate::{markers_near_focus, primary_marker, primary_solar_system};
    use geoforge_types::cosmology::CosmologyPreset;
    use geoforge_types::{PipelineLayerId, Seed};

    #[test]
    fn minimal_is_deterministic() {
        let a = run_pre_tectonic_layers(
            Seed::new(42),
            CosmologyPreset::Minimal,
            PipelineLayerId::PlanetaryBody,
        );
        let b = run_pre_tectonic_layers(
            Seed::new(42),
            CosmologyPreset::Minimal,
            PipelineLayerId::PlanetaryBody,
        );
        let ma = primary_marker(&a);
        let mb = primary_marker(&b);
        assert_eq!(ma.position_ly.x_ly, mb.position_ly.x_ly);
        assert_eq!(ma.phenotype.primary_class, mb.phenotype.primary_class);
        assert_eq!(
            a.planetary.as_ref().unwrap().focus_slot_index,
            b.planetary.as_ref().unwrap().focus_slot_index
        );
    }

    #[test]
    fn solar_system_matches_marker_position() {
        let layers = run_pre_tectonic_layers(
            Seed::new(99),
            CosmologyPreset::Minimal,
            PipelineLayerId::SolarSystem,
        );
        let m = primary_marker(&layers);
        let sys = primary_solar_system(&layers);
        assert_eq!(m.position_ly.x_ly, sys.barycenter_ly.x_ly);
        assert_eq!(m.phenotype.multiplicity, sys.phenotype.multiplicity);
    }

    #[test]
    fn rich_loads_more_markers_than_minimal() {
        let min = run_pre_tectonic_layers(Seed::new(1), CosmologyPreset::Minimal, PipelineLayerId::Galaxy);
        let rich = run_pre_tectonic_layers(Seed::new(1), CosmologyPreset::Rich, PipelineLayerId::Galaxy);
        assert!(markers_near_focus(&rich, 0).len() >= markers_near_focus(&min, 0).len());
    }
}
