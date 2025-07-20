#[cfg(test)]
mod tests {
    use crate::map::spherical::PlanetaryParams;

    #[test]
    fn test_earth_parameters() {
        let earth = PlanetaryParams::earth();
        
        assert_eq!(earth.radius_km, 6371.0);
        assert_eq!(earth.gravity_ms2, 9.81);
        assert_eq!(earth.axial_tilt_degrees, 23.44);
        assert_eq!(earth.rotation_period_hours, 24.0);
        assert_eq!(earth.atmospheric_pressure_kpa, 101.325);
        assert_eq!(earth.greenhouse_factor, 1.0);
        assert_eq!(earth.stellar_luminosity, 1.0);
    }
    
    #[test]
    fn test_mars_parameters() {
        let mars = PlanetaryParams::mars();
        
        assert_eq!(mars.radius_km, 3389.5);
        assert_eq!(mars.gravity_ms2, 3.71);
        assert_eq!(mars.axial_tilt_degrees, 25.19);
        assert_eq!(mars.rotation_period_hours, 24.62);
        assert_eq!(mars.atmospheric_pressure_kpa, 0.636);
        assert_eq!(mars.greenhouse_factor, 0.18);
        assert_eq!(mars.stellar_luminosity, 1.0);
    }
    
    #[test]
    fn test_venus_parameters() {
        let venus = PlanetaryParams::venus();
        
        assert_eq!(venus.radius_km, 6051.8);
        assert_eq!(venus.gravity_ms2, 8.87);
        assert_eq!(venus.axial_tilt_degrees, 177.4);
        assert_eq!(venus.rotation_period_hours, 5832.5);
        assert_eq!(venus.atmospheric_pressure_kpa, 9200.0);
        assert_eq!(venus.greenhouse_factor, 15.15);
        assert_eq!(venus.stellar_luminosity, 1.0);
    }
    
    #[test]
    fn test_custom_from_radius() {
        let custom = PlanetaryParams::from_radius(5000.0);
        
        assert_eq!(custom.radius_km, 5000.0);
        assert!((custom.surface_area_km2 - (4.0 * std::f64::consts::PI * 5000.0 * 5000.0)).abs() < 1.0);
        // Should have Earth-like density but scaled
        assert_eq!(custom.density_kgm3, 5515.0);
    }
    
    #[test]
    fn test_distance_conversion() {
        let earth = PlanetaryParams::earth();
        let mars = PlanetaryParams::mars();
        
        let radians = std::f64::consts::PI / 4.0; // 45 degrees
        
        let earth_km = earth.radians_to_km(radians);
        let mars_km = mars.radians_to_km(radians);
        
        // Mars should have smaller distance for same angle due to smaller radius
        assert!(mars_km < earth_km);
        
        // Test round trip
        let back_to_radians = earth.km_to_radians(earth_km);
        assert!((back_to_radians - radians).abs() < 1e-10);
    }
    
    #[test]
    fn test_escape_velocity() {
        let earth = PlanetaryParams::earth();
        let mars = PlanetaryParams::mars();
        
        let earth_escape = earth.escape_velocity_kms();
        let mars_escape = mars.escape_velocity_kms();
        
        // Earth should have higher escape velocity
        assert!(earth_escape > mars_escape);
        
        // Earth's escape velocity should be around 11.2 km/s
        assert!((earth_escape - 11.2).abs() < 0.5);
    }
    
    #[test]
    fn test_seasonal_variation() {
        let earth = PlanetaryParams::earth();
        let venus = PlanetaryParams::venus();
        
        let earth_seasonal = earth.seasonal_variation_factor();
        let venus_seasonal = venus.seasonal_variation_factor();
        
        // Venus has extreme axial tilt (retrograde) so should have max variation
        assert_eq!(venus_seasonal, 1.0);
        
        // Earth should have moderate seasonal variation
        assert!(earth_seasonal > 0.0 && earth_seasonal < 1.0);
    }
    
    #[test]
    fn test_diurnal_variation() {
        let earth = PlanetaryParams::earth();
        let venus = PlanetaryParams::venus();
        
        let earth_diurnal = earth.diurnal_variation_factor();
        let venus_diurnal = venus.diurnal_variation_factor();
        
        // Venus with its very long day should have high diurnal variation
        assert!(venus_diurnal > earth_diurnal);
        
        // Fast-rotating planet should have low variation
        let fast_planet = PlanetaryParams::from_radius_and_density(6371.0, 5515.0);
        let mut fast_planet = fast_planet;
        fast_planet.rotation_period_hours = 8.0; // 8-hour day
        let fast_diurnal = fast_planet.diurnal_variation_factor();
        assert!(fast_diurnal < earth_diurnal);
    }
    
    #[test]
    fn test_inverse_square_law() {
        // Test fundamental radiation physics: flux ∝ 1/distance²
        let mut planet = PlanetaryParams::earth();
        
        // Baseline at 1 AU
        planet.semi_major_axis_au = 1.0;
        planet.perihelion_au = 1.0;
        planet.aphelion_au = 1.0;
        let baseline_insolation = planet.average_insolation_wm2();
        
        // Double the distance -> should get 1/4 the insolation
        planet.semi_major_axis_au = 2.0;
        planet.perihelion_au = 2.0;
        planet.aphelion_au = 2.0;
        let double_distance_insolation = planet.average_insolation_wm2();
        
        let expected_quarter = baseline_insolation / 4.0;
        assert!((double_distance_insolation - expected_quarter).abs() < 1.0);
        
        // Triple the distance -> should get 1/9 the insolation
        planet.semi_major_axis_au = 3.0;
        planet.perihelion_au = 3.0;
        planet.aphelion_au = 3.0;
        let triple_distance_insolation = planet.average_insolation_wm2();
        
        let expected_ninth = baseline_insolation / 9.0;
        assert!((triple_distance_insolation - expected_ninth).abs() < 1.0);
    }
    
    #[test]
    fn test_stellar_luminosity_scaling() {
        // Test luminosity scaling: flux ∝ luminosity
        let mut planet = PlanetaryParams::earth();
        planet.semi_major_axis_au = 1.0;
        planet.perihelion_au = 1.0;
        planet.aphelion_au = 1.0;
        
        // Baseline with solar luminosity
        planet.stellar_luminosity = 1.0;
        let baseline_insolation = planet.average_insolation_wm2();
        
        // Double stellar luminosity -> should get double insolation
        planet.stellar_luminosity = 2.0;
        let double_luminosity_insolation = planet.average_insolation_wm2();
        assert!((double_luminosity_insolation - 2.0 * baseline_insolation).abs() < 1.0);
        
        // Half stellar luminosity -> should get half insolation  
        planet.stellar_luminosity = 0.5;
        let half_luminosity_insolation = planet.average_insolation_wm2();
        assert!((half_luminosity_insolation - 0.5 * baseline_insolation).abs() < 1.0);
        
        // Ten times stellar luminosity -> should get 10x insolation
        planet.stellar_luminosity = 10.0;
        let ten_times_luminosity_insolation = planet.average_insolation_wm2();
        assert!((ten_times_luminosity_insolation - 10.0 * baseline_insolation).abs() < 1.0);
    }
    
    #[test]
    fn test_combined_distance_luminosity_effects() {
        // Test combined effects: flux ∝ luminosity/distance²
        let mut planet = PlanetaryParams::earth();
        
        // Scenario: 2x luminosity at 2x distance should give 1/2 the baseline insolation
        planet.stellar_luminosity = 1.0;
        planet.semi_major_axis_au = 1.0;
        planet.perihelion_au = 1.0;
        planet.aphelion_au = 1.0;
        let baseline = planet.average_insolation_wm2();
        
        planet.stellar_luminosity = 2.0;  // 2x brighter star
        planet.semi_major_axis_au = 2.0;  // 2x farther away
        planet.perihelion_au = 2.0;
        planet.aphelion_au = 2.0;
        let combined_effect = planet.average_insolation_wm2();
        
        // 2x luminosity / (2x distance)² = 2x / 4x = 0.5x
        let expected = baseline * 0.5;
        assert!((combined_effect - expected).abs() < 1.0);
    }
}