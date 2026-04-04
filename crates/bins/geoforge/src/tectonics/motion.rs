//! Plate motion assignment for tectonic plates
//!
//! This module handles assigning realistic motion vectors to tectonic plates,
//! including direction (azimuth) and velocity (cm/year). Motion can be assigned
//! randomly or using physics-based angular velocity models.

use crate::tectonics::plates::PlateSeed;
use crate::map::spherical::SphericalPoint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Configuration for plate motion assignment
#[derive(Debug, Clone)]
pub struct PlateMotionConfig {
    /// Random seed for deterministic generation
    pub seed: u64,

    /// Minimum plate speed in cm/year (default: 1.0)
    pub min_speed_cm_year: f64,

    /// Maximum plate speed in cm/year (default: 10.0)
    pub max_speed_cm_year: f64,

    /// Use physics-based angular velocity model (default: false)
    /// If false, uses random assignment
    pub use_physics_model: bool,

    /// Global rotation rate modifier (default: 1.0)
    /// Scales all velocities by this factor
    pub global_speed_scale: f64,
}

impl Default for PlateMotionConfig {
    fn default() -> Self {
        Self {
            seed: 0,
            min_speed_cm_year: 1.0,
            max_speed_cm_year: 10.0,
            use_physics_model: false,
            global_speed_scale: 1.0,
        }
    }
}

impl PlateMotionConfig {
    /// Create a new configuration with a specific seed
    pub fn with_seed(seed: u64) -> Self {
        Self {
            seed,
            ..Default::default()
        }
    }

    /// Set the speed range for plates
    pub fn with_speed_range(mut self, min: f64, max: f64) -> Self {
        self.min_speed_cm_year = min;
        self.max_speed_cm_year = max;
        self
    }

    /// Enable physics-based angular velocity model
    pub fn with_physics_model(mut self, enabled: bool) -> Self {
        self.use_physics_model = enabled;
        self
    }

    /// Set global speed scale factor
    pub fn with_speed_scale(mut self, scale: f64) -> Self {
        self.global_speed_scale = scale;
        self
    }
}

/// Assigns motion vectors to tectonic plates
pub struct PlateMotionAssigner {
    config: PlateMotionConfig,
    rng: StdRng,
}

impl PlateMotionAssigner {
    /// Create a new motion assigner with default configuration
    pub fn new() -> Self {
        Self::with_config(PlateMotionConfig::default())
    }

    /// Create a new motion assigner with custom configuration
    pub fn with_config(config: PlateMotionConfig) -> Self {
        let rng = StdRng::seed_from_u64(config.seed);
        Self { config, rng }
    }

    /// Assign motion vectors to all plates
    ///
    /// Modifies the plate seeds in-place, assigning realistic motion_direction
    /// and motion_speed values based on the configuration.
    ///
    /// # Arguments
    /// * `seeds` - Mutable slice of plate seeds to assign motion to
    pub fn assign_motion(&mut self, seeds: &mut [PlateSeed]) {
        if self.config.use_physics_model {
            self.assign_physics_based_motion(seeds);
        } else {
            self.assign_random_motion(seeds);
        }
    }

    /// Assign random motion vectors to plates
    ///
    /// Each plate gets a random direction (0-360°) and speed within the configured range.
    /// This is simpler but still produces reasonable tectonic behavior.
    fn assign_random_motion(&mut self, seeds: &mut [PlateSeed]) {
        for seed in seeds.iter_mut() {
            // Random direction (0-360 degrees)
            seed.motion_direction = self.rng.gen_range(0.0..360.0);

            // Random speed within range
            let speed = self.rng.gen_range(
                self.config.min_speed_cm_year..=self.config.max_speed_cm_year
            );
            seed.motion_speed = speed * self.config.global_speed_scale;
        }
    }

    /// Assign physics-based motion using angular velocity model
    ///
    /// Models plate motion as rotation about a pole on the sphere.
    /// More realistic but also more complex.
    fn assign_physics_based_motion(&mut self, seeds: &mut [PlateSeed]) {
        // For each plate, assign a random rotation pole and angular velocity
        for seed in seeds.iter_mut() {
            // Random rotation pole (latitude and longitude)
            let pole_lat = self.rng.gen_range(-90.0..=90.0);
            let pole_lon = self.rng.gen_range(-180.0..=180.0);
            let pole = SphericalPoint::from_lat_lon(pole_lat, pole_lon);

            // Random angular velocity (degrees per million years)
            // Converted to linear velocity at plate location
            let angular_vel_deg_per_my = self.rng.gen_range(0.1..=1.0);

            // Calculate motion direction and speed at this plate's location
            let (direction, speed) = self.calculate_motion_from_pole(
                seed.spherical_point(),
                &pole,
                angular_vel_deg_per_my,
            );

            seed.motion_direction = direction;
            seed.motion_speed = speed * self.config.global_speed_scale;
        }
    }

    /// Calculate motion direction and speed at a point given a rotation pole
    ///
    /// Returns (direction_degrees, speed_cm_per_year)
    fn calculate_motion_from_pole(
        &self,
        point: &SphericalPoint,
        pole: &SphericalPoint,
        angular_vel_deg_per_my: f64,
    ) -> (f64, f64) {
        use crate::map::spherical::bearing;

        // Angular distance from pole (in radians)
        let angular_distance = point.distance_to(pole);

        // Motion is perpendicular to the great circle connecting point and pole
        // Direction can be calculated using spherical trigonometry
        let point_to_pole_bearing = bearing(point, pole);

        // Motion direction is perpendicular to bearing to pole
        let motion_direction = (point_to_pole_bearing + 90.0) % 360.0;

        // Speed depends on distance from pole: v = ω × r × sin(θ)
        // Where θ is the angular distance from pole
        const EARTH_RADIUS_KM: f64 = 6371.0;
        let speed_km_per_my = angular_vel_deg_per_my.to_radians()
            * EARTH_RADIUS_KM
            * angular_distance.sin();

        // Convert km/million years to cm/year
        let speed_cm_per_year = speed_km_per_my * 100_000.0 / 1_000_000.0;

        // Clamp to reasonable range
        let speed_clamped = speed_cm_per_year.clamp(
            self.config.min_speed_cm_year,
            self.config.max_speed_cm_year,
        );

        (motion_direction, speed_clamped)
    }

    /// Calculate relative motion between two plates at a specific point
    ///
    /// Returns (relative_velocity_cm_per_year, angle_to_boundary_degrees)
    ///
    /// The angle indicates the direction of relative motion:
    /// - Near 0° or 180°: Motion perpendicular to boundary (convergent or divergent)
    /// - Near 90°: Motion parallel to boundary (transform)
    pub fn calculate_relative_motion(
        seed_a: &PlateSeed,
        seed_b: &PlateSeed,
        _boundary_point: &SphericalPoint,
    ) -> (f64, f64) {
        // Calculate velocity vectors in local tangent plane
        let vel_a = velocity_vector(seed_a);
        let vel_b = velocity_vector(seed_b);

        // Relative velocity: v_rel = v_a - v_b
        let rel_vel_x = vel_a.0 - vel_b.0;
        let rel_vel_y = vel_a.1 - vel_b.1;

        // Magnitude of relative velocity
        let rel_speed = (rel_vel_x * rel_vel_x + rel_vel_y * rel_vel_y).sqrt();

        // Direction of relative velocity
        let rel_direction = rel_vel_y.atan2(rel_vel_x).to_degrees();

        use crate::map::spherical::bearing;

        // Approximate boundary direction (perpendicular to line connecting seeds)
        let boundary_bearing = bearing(seed_a.spherical_point(), seed_b.spherical_point());
        let boundary_direction = (boundary_bearing + 90.0) % 360.0;

        // Angle between relative velocity and boundary
        let mut angle = (rel_direction - boundary_direction).abs();
        if angle > 180.0 {
            angle = 360.0 - angle;
        }
        // Normalize to 0-90 range (acute angle)
        if angle > 90.0 {
            angle = 180.0 - angle;
        }

        (rel_speed, angle)
    }

    /// Classify boundary type based on relative motion
    ///
    /// Returns the boundary interaction type based on the angle of relative motion
    /// and whether plates are moving toward or away from each other.
    pub fn classify_boundary_interaction(
        seed_a: &PlateSeed,
        seed_b: &PlateSeed,
        angle_threshold: f64,
    ) -> crate::tectonics::plates::PlateInteraction {
        use crate::tectonics::plates::PlateInteraction;

        // Get velocity vectors
        let vel_a = velocity_vector(seed_a);
        let vel_b = velocity_vector(seed_b);

        // Vector from A to B (approximate)
        let dx = seed_b.lon - seed_a.lon;
        let dy = seed_b.lat - seed_a.lat;
        let dist = (dx * dx + dy * dy).sqrt();

        if dist < 1e-10 {
            return PlateInteraction::Transform; // Degenerate case
        }

        // Normalized direction vector
        let nx = dx / dist;
        let ny = dy / dist;

        // Relative velocity
        let rel_vx = vel_a.0 - vel_b.0;
        let rel_vy = vel_a.1 - vel_b.1;

        // Dot product: positive = converging, negative = diverging
        let dot_product = rel_vx * nx + rel_vy * ny;

        // Cross product magnitude (for angle calculation)
        let cross_product = (rel_vx * ny - rel_vy * nx).abs();
        let parallel_component = dot_product.abs();

        // Angle from perpendicular motion
        let angle = cross_product.atan2(parallel_component).to_degrees();

        // Classify based on angle and dot product
        if angle < angle_threshold {
            // Nearly perpendicular motion
            if dot_product > 0.0 {
                PlateInteraction::Convergent // Moving toward each other
            } else {
                PlateInteraction::Divergent // Moving apart
            }
        } else {
            PlateInteraction::Transform // Sliding past
        }
    }
}

impl Default for PlateMotionAssigner {
    fn default() -> Self {
        Self::new()
    }
}

/// Helper function to convert plate seed to velocity vector
fn velocity_vector(seed: &PlateSeed) -> (f64, f64) {
    let dir_rad = seed.motion_direction.to_radians();
    let vx = seed.motion_speed * dir_rad.cos();
    let vy = seed.motion_speed * dir_rad.sin();
    (vx, vy)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_creation() {
        let config = PlateMotionConfig::default();
        assert_eq!(config.min_speed_cm_year, 1.0);
        assert_eq!(config.max_speed_cm_year, 10.0);
        assert!(!config.use_physics_model);

        let config = PlateMotionConfig::with_seed(42)
            .with_speed_range(2.0, 8.0)
            .with_physics_model(true)
            .with_speed_scale(1.5);

        assert_eq!(config.seed, 42);
        assert_eq!(config.min_speed_cm_year, 2.0);
        assert_eq!(config.max_speed_cm_year, 8.0);
        assert!(config.use_physics_model);
        assert_eq!(config.global_speed_scale, 1.5);
    }

    #[test]
    fn test_random_motion_assignment() {
        let config = PlateMotionConfig::with_seed(123);
        let mut assigner = PlateMotionAssigner::with_config(config);

        let mut seeds = vec![
            PlateSeed::new(1, 0, 0, 0.0, 0.0, 0.0, 0.0),
            PlateSeed::new(2, 100, 50, 45.0, -120.0, 0.0, 0.0),
            PlateSeed::new(3, 200, 100, -30.0, 60.0, 0.0, 0.0),
        ];

        assigner.assign_motion(&mut seeds);

        // All plates should have motion assigned
        for seed in &seeds {
            assert!(seed.motion_direction >= 0.0 && seed.motion_direction < 360.0);
            assert!(seed.motion_speed >= 1.0 && seed.motion_speed <= 10.0);
        }

        // Different plates should have different motion (with high probability)
        assert_ne!(seeds[0].motion_direction, seeds[1].motion_direction);
        assert_ne!(seeds[1].motion_speed, seeds[2].motion_speed);
    }

    #[test]
    fn test_motion_determinism() {
        let config = PlateMotionConfig::with_seed(456);

        let mut seeds1 = vec![
            PlateSeed::new(1, 0, 0, 0.0, 0.0, 0.0, 0.0),
            PlateSeed::new(2, 100, 50, 45.0, -120.0, 0.0, 0.0),
        ];
        let mut seeds2 = seeds1.clone();

        let mut assigner1 = PlateMotionAssigner::with_config(config.clone());
        let mut assigner2 = PlateMotionAssigner::with_config(config);

        assigner1.assign_motion(&mut seeds1);
        assigner2.assign_motion(&mut seeds2);

        // Same seed should produce identical results
        for (s1, s2) in seeds1.iter().zip(seeds2.iter()) {
            assert_eq!(s1.motion_direction, s2.motion_direction);
            assert_eq!(s1.motion_speed, s2.motion_speed);
        }
    }

    #[test]
    fn test_boundary_classification() {
        use crate::tectonics::plates::PlateInteraction;

        // Converging plates
        let seed_a = PlateSeed::new(1, 0, 0, 0.0, 0.0, 0.0, 5.0); // Moving east
        let seed_b = PlateSeed::new(2, 100, 0, 0.0, 10.0, 180.0, 5.0); // Moving west

        let interaction = PlateMotionAssigner::classify_boundary_interaction(&seed_a, &seed_b, 45.0);
        assert_eq!(interaction, PlateInteraction::Convergent);

        // Diverging plates
        let seed_c = PlateSeed::new(3, 0, 0, 0.0, 0.0, 180.0, 5.0); // Moving west
        let seed_d = PlateSeed::new(4, 100, 0, 0.0, 10.0, 0.0, 5.0); // Moving east

        let interaction = PlateMotionAssigner::classify_boundary_interaction(&seed_c, &seed_d, 45.0);
        assert_eq!(interaction, PlateInteraction::Divergent);

        // Transform boundary
        let seed_e = PlateSeed::new(5, 0, 0, 0.0, 0.0, 90.0, 5.0); // Moving north
        let seed_f = PlateSeed::new(6, 100, 0, 0.0, 10.0, 270.0, 5.0); // Moving south

        let interaction = PlateMotionAssigner::classify_boundary_interaction(&seed_e, &seed_f, 45.0);
        assert_eq!(interaction, PlateInteraction::Transform);
    }

    #[test]
    fn test_speed_range_respected() {
        let config = PlateMotionConfig::with_seed(789)
            .with_speed_range(3.0, 6.0);

        let mut assigner = PlateMotionAssigner::with_config(config);

        let mut seeds = vec![
            PlateSeed::new(1, 0, 0, 0.0, 0.0, 0.0, 0.0),
            PlateSeed::new(2, 100, 50, 45.0, -120.0, 0.0, 0.0),
            PlateSeed::new(3, 200, 100, -30.0, 60.0, 0.0, 0.0),
            PlateSeed::new(4, 300, 150, 60.0, 120.0, 0.0, 0.0),
        ];

        assigner.assign_motion(&mut seeds);

        for seed in &seeds {
            assert!(seed.motion_speed >= 3.0, "Speed {} below minimum", seed.motion_speed);
            assert!(seed.motion_speed <= 6.0, "Speed {} above maximum", seed.motion_speed);
        }
    }
}
