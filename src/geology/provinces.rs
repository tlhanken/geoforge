//! Geological province type definitions and characteristics
//!
//! Defines the types of geological provinces that can be generated from tectonic data.

/// Geological province types
///
/// These represent major geological domains created by tectonic processes.
/// Stage 2.1 focuses on orogenic belts; future stages will add other province types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GeologicProvince {
    // ========== Stage 2.1: Orogenic Belts ==========

    /// Collision orogen - Continental-continental convergence
    ///
    /// Forms when two continental plates collide, creating massive mountain ranges.
    /// Example: Himalayas (India-Eurasia collision)
    ///
    /// Characteristics:
    /// - Highest elevation potential
    /// - Broad, thick crustal zone
    /// - Significant crustal shortening
    CollisionOrogen,

    /// Subduction orogen - Oceanic-continental convergence
    ///
    /// Forms where oceanic crust subducts beneath continental crust,
    /// creating volcanic mountain ranges and deep ocean trenches.
    /// Example: Andes (Nazca-South America), Cascades (Juan de Fuca-North America)
    ///
    /// Characteristics:
    /// - Volcanic arc parallel to trench
    /// - Moderate to high elevation
    /// - Associated with deep ocean trench
    ///
    /// Note: Oceanic-oceanic convergence (island arcs like Marianas) will be
    /// handled in a future stage as part of Arc and Basin Systems.
    SubductionOrogen,

    // ========== Future Stages ==========
    // These will be implemented in later stages

    // Stage 2.2: Large Igneous Provinces
    // ContinentalFloodBasalt,
    // OceanicPlateau,
    // HotspotTrack,

    // Stage 2.3: Arc and Basin Systems
    // VolcanicArc,
    // ForearcBasin,
    // BackarcBasin,

    // Stage 2.4: Stable Continental Regions
    // Craton,
    // Platform,
    // IntracratonicBasin,

    // Stage 2.5: Extensional Zones
    // ContinentalRift,
    // ExtendedCrust,

    // Stage 2.6: Oceanic Domains
    // MidOceanRidge,
    // AbyssalPlain,
    // OceanTrench,
}

impl GeologicProvince {
    /// Get a human-readable name for this province type
    pub fn name(&self) -> &'static str {
        match self {
            GeologicProvince::CollisionOrogen => "Collision Orogen",
            GeologicProvince::SubductionOrogen => "Subduction Orogen",
        }
    }

    /// Get a short description of this province type
    pub fn description(&self) -> &'static str {
        match self {
            GeologicProvince::CollisionOrogen =>
                "Continental-continental collision zone (e.g., Himalayas)",
            GeologicProvince::SubductionOrogen =>
                "Oceanic-continental subduction zone (e.g., Andes)",
        }
    }
}

/// Characteristics and parameters of a geological province
///
/// These define how a province affects terrain generation, elevation, and other properties.
#[derive(Debug, Clone)]
pub struct ProvinceCharacteristics {
    /// Type of geological province
    pub province_type: GeologicProvince,

    /// Elevation intensity factor (0.0 to 1.0+)
    ///
    /// Multiplier for base elevation in this province:
    /// - 1.0 = maximum elevation potential (collision orogens)
    /// - 0.7 = moderate elevation (subduction orogens)
    /// - 0.5 = lower elevation (accretionary orogens)
    pub elevation_intensity: f64,

    /// Terrain roughness factor (0.0 to 1.0)
    ///
    /// How irregular/rugged the terrain is:
    /// - High roughness: sharp peaks, deep valleys
    /// - Low roughness: smoother, more uniform
    pub roughness: f64,

    /// Typical width of the province in kilometers
    ///
    /// Calculated dynamically from convergence rate and other factors.
    /// - Collision: 500-2000 km
    /// - Subduction: 200-800 km
    pub width_km: f64,

    /// Activity/intensity level (0.0 to 1.0)
    ///
    /// How active or strong the tectonic processes are:
    /// - Derived from convergence rate
    /// - Higher = more active mountain building
    /// - Lower = older, less active
    pub intensity: f64,

    /// Average convergence rate in cm/year
    ///
    /// Used to calculate dynamic width and intensity
    pub convergence_rate: f64,
}

impl ProvinceCharacteristics {
    /// Create characteristics for a collision orogen
    ///
    /// # Arguments
    /// * `convergence_rate` - Convergence rate in cm/year
    /// * `width_km` - Calculated width in kilometers
    pub fn collision_orogen(convergence_rate: f64, width_km: f64) -> Self {
        let intensity = (convergence_rate / 10.0).min(1.0); // Normalize to 0-1

        Self {
            province_type: GeologicProvince::CollisionOrogen,
            elevation_intensity: 1.0, // Highest elevation potential
            roughness: 0.9,           // Very rugged
            width_km,
            intensity,
            convergence_rate,
        }
    }

    /// Create characteristics for a subduction orogen
    ///
    /// # Arguments
    /// * `convergence_rate` - Convergence rate in cm/year
    /// * `width_km` - Calculated width in kilometers
    pub fn subduction_orogen(convergence_rate: f64, width_km: f64) -> Self {
        let intensity = (convergence_rate / 10.0).min(1.0);

        Self {
            province_type: GeologicProvince::SubductionOrogen,
            elevation_intensity: 0.7, // Moderate to high elevation
            roughness: 0.8,           // Quite rugged
            width_km,
            intensity,
            convergence_rate,
        }
    }
}

/// A region of specific geological province
///
/// Represents a contiguous area with uniform geological characteristics,
/// typically generated from a tectonic boundary.
#[derive(Debug, Clone)]
pub struct ProvinceRegion {
    /// Pixel coordinates that belong to this province
    pub pixels: Vec<(usize, usize)>,

    /// Geological characteristics of this region
    pub characteristics: ProvinceCharacteristics,

    /// Index of the boundary that created this province (if applicable)
    ///
    /// Links back to the BoundarySegment that generated this orogenic belt.
    /// None for provinces not generated from boundaries.
    pub source_boundary_index: Option<usize>,
}

impl ProvinceRegion {
    /// Create a new province region
    pub fn new(
        pixels: Vec<(usize, usize)>,
        characteristics: ProvinceCharacteristics,
        source_boundary_index: Option<usize>,
    ) -> Self {
        Self {
            pixels,
            characteristics,
            source_boundary_index,
        }
    }

    /// Get the number of pixels in this region
    pub fn pixel_count(&self) -> usize {
        self.pixels.len()
    }

    /// Get the approximate area in square kilometers
    ///
    /// # Arguments
    /// * `km_per_pixel` - The average kilometers per pixel for the map
    pub fn area_km2(&self, km_per_pixel: f64) -> f64 {
        self.pixels.len() as f64 * km_per_pixel * km_per_pixel
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_province_names() {
        assert_eq!(GeologicProvince::CollisionOrogen.name(), "Collision Orogen");
        assert_eq!(GeologicProvince::SubductionOrogen.name(), "Subduction Orogen");
    }

    #[test]
    fn test_collision_characteristics() {
        let chars = ProvinceCharacteristics::collision_orogen(5.0, 1000.0);
        assert_eq!(chars.province_type, GeologicProvince::CollisionOrogen);
        assert_eq!(chars.elevation_intensity, 1.0);
        assert_eq!(chars.width_km, 1000.0);
        assert_eq!(chars.convergence_rate, 5.0);
        assert!(chars.intensity > 0.0 && chars.intensity <= 1.0);
    }

    #[test]
    fn test_subduction_characteristics() {
        let chars = ProvinceCharacteristics::subduction_orogen(8.0, 500.0);
        assert_eq!(chars.province_type, GeologicProvince::SubductionOrogen);
        assert_eq!(chars.elevation_intensity, 0.7);
        assert_eq!(chars.width_km, 500.0);
        assert!(chars.intensity > 0.0 && chars.intensity <= 1.0);
    }

    #[test]
    fn test_province_region() {
        let pixels = vec![(0, 0), (1, 0), (2, 0)];
        let chars = ProvinceCharacteristics::collision_orogen(5.0, 1000.0);
        let region = ProvinceRegion::new(pixels, chars, Some(0));

        assert_eq!(region.pixel_count(), 3);
        assert_eq!(region.source_boundary_index, Some(0));

        // Test area calculation (assuming 10 km/pixel)
        let area = region.area_km2(10.0);
        assert_eq!(area, 300.0); // 3 pixels * 10^2
    }
}
