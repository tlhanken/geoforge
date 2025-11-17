//! Geological province type definitions and characteristics
//!
//! Defines the types of geological provinces that can be generated from tectonic data.

/// Geological province types
///
/// These represent major geological domains created by tectonic processes.
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
    SubductionOrogen,

    /// Accretionary orogen - Terrane accretion
    ///
    /// Forms when exotic terranes (island arcs, oceanic plateaus, microcontinents)
    /// are scraped off and accreted onto continental margins during subduction.
    /// Example: Western North America (Cordilleran orogeny)
    ///
    /// Characteristics:
    /// - Complex, deformed geology with mixed origins
    /// - Moderate elevation
    /// - Irregular structure
    AccretionaryOrogen,

    /// Extensional orogen - Core complex formation
    ///
    /// Forms in zones of continental extension where metamorphic core complexes
    /// are exhumed. Often associated with rifting but different from active rifts.
    /// Example: Basin and Range Province (western USA)
    ///
    /// Characteristics:
    /// - Moderate relief with uplifted metamorphic cores
    /// - Normal faulting and extension
    /// - Often follows earlier compression
    ExtensionalOrogen,

    // ========== Stage 2.2: Large Igneous Provinces ==========

    /// Continental flood basalt - Massive volcanic province on continental crust
    ///
    /// Formed by mantle plume eruptions creating vast lava plains.
    /// Example: Deccan Traps (India), Siberian Traps (Russia)
    ContinentalFloodBasalt,

    /// Oceanic plateau - Massive volcanic province on oceanic crust
    ///
    /// Formed by underwater mantle plume eruptions.
    /// Example: Ontong Java Plateau (Pacific Ocean)
    OceanicPlateau,

    /// Hotspot track - Volcanic island/seamount chain
    ///
    /// Formed as a tectonic plate moves over a stationary mantle plume.
    /// Example: Hawaiian Islands, Yellowstone track
    HotspotTrack,

    // ========== Stage 2.3: Arc and Basin Systems ==========

    /// Volcanic arc - Active subduction zone volcanism
    ///
    /// Chain of volcanoes above a subduction zone.
    /// Example: Aleutian Islands, Java (Indonesia)
    VolcanicArc,

    /// Forearc basin - Sedimentary basin between trench and arc
    ///
    /// Sediment accumulation zone on the oceanic side of a volcanic arc.
    /// Example: Great Valley (California)
    ForearcBasin,

    /// Backarc basin - Extensional basin behind volcanic arc
    ///
    /// Forms behind volcanic arcs due to slab rollback or mantle wedge dynamics.
    /// Example: Sea of Japan, Lau Basin
    BackarcBasin,

    // ========== Stage 2.4: Stable Continental Regions ==========

    /// Craton/Shield - Ancient stable continental core
    ///
    /// Very old (>1.5 Ga), tectonically stable continental crust with low relief.
    /// Example: Canadian Shield, Kaapvaal Craton (South Africa)
    Craton,

    /// Platform - Stable craton with sedimentary cover
    ///
    /// Cratonic regions covered by thin sedimentary layers.
    /// Example: Russian Platform, North American Platform
    Platform,

    /// Intracratonic basin - Subsided region within stable craton
    ///
    /// Basin formed by slow subsidence within an otherwise stable craton.
    /// Example: Michigan Basin, Illinois Basin
    IntracratonicBasin,

    // ========== Stage 2.5: Extensional Zones ==========

    /// Continental rift - Active extensional zone
    ///
    /// Zone where continental crust is actively being pulled apart.
    /// Example: East African Rift, Rio Grande Rift
    ContinentalRift,

    /// Extended crust - Thinned continental crust from extension
    ///
    /// Areas where continental crust has been stretched and thinned.
    /// Example: Basin and Range (western USA), Aegean Sea
    ExtendedCrust,

    // ========== Stage 2.6: Oceanic Domains ==========

    /// Mid-ocean ridge - Active seafloor spreading center
    ///
    /// Divergent boundary where new oceanic crust is created.
    /// Example: Mid-Atlantic Ridge, East Pacific Rise
    MidOceanRidge,

    /// Abyssal plain - Deep oceanic basin
    ///
    /// Flat, sediment-covered regions of the deep ocean floor.
    /// Example: Most of the deep Pacific, Atlantic basins
    AbyssalPlain,

    /// Ocean trench - Deep subduction zone
    ///
    /// Deepest parts of the ocean where oceanic crust subducts.
    /// Example: Mariana Trench, Peru-Chile Trench
    OceanTrench,

    /// Oceanic fracture zone - Transform fault system
    ///
    /// Large transform faults that offset mid-ocean ridge segments.
    /// Example: Romanche Fracture Zone, Charlie-Gibbs Fracture Zone
    OceanicFractureZone,

    /// Seamount field - Scattered underwater volcanoes
    ///
    /// Individual or clustered submarine volcanoes on oceanic crust.
    /// Different from hotspot tracks (which are linear chains).
    /// Example: Musicians Seamounts, Corner Rise Seamounts
    SeamountField,
}

impl GeologicProvince {
    /// Get a human-readable name for this province type
    pub fn name(&self) -> &'static str {
        match self {
            // Stage 2.1: Orogenic Belts
            GeologicProvince::CollisionOrogen => "Collision Orogen",
            GeologicProvince::SubductionOrogen => "Subduction Orogen",
            GeologicProvince::AccretionaryOrogen => "Accretionary Orogen",
            GeologicProvince::ExtensionalOrogen => "Extensional Orogen",

            // Stage 2.2: Large Igneous Provinces
            GeologicProvince::ContinentalFloodBasalt => "Continental Flood Basalt",
            GeologicProvince::OceanicPlateau => "Oceanic Plateau",
            GeologicProvince::HotspotTrack => "Hotspot Track",

            // Stage 2.3: Arc and Basin Systems
            GeologicProvince::VolcanicArc => "Volcanic Arc",
            GeologicProvince::ForearcBasin => "Forearc Basin",
            GeologicProvince::BackarcBasin => "Backarc Basin",

            // Stage 2.4: Stable Continental Regions
            GeologicProvince::Craton => "Craton/Shield",
            GeologicProvince::Platform => "Platform",
            GeologicProvince::IntracratonicBasin => "Intracratonic Basin",

            // Stage 2.5: Extensional Zones
            GeologicProvince::ContinentalRift => "Continental Rift",
            GeologicProvince::ExtendedCrust => "Extended Crust",

            // Stage 2.6: Oceanic Domains
            GeologicProvince::MidOceanRidge => "Mid-Ocean Ridge",
            GeologicProvince::AbyssalPlain => "Abyssal Plain",
            GeologicProvince::OceanTrench => "Ocean Trench",
            GeologicProvince::OceanicFractureZone => "Oceanic Fracture Zone",
            GeologicProvince::SeamountField => "Seamount Field",
        }
    }

    /// Get a short description of this province type
    pub fn description(&self) -> &'static str {
        match self {
            // Stage 2.1: Orogenic Belts
            GeologicProvince::CollisionOrogen =>
                "Continental-continental collision zone (e.g., Himalayas)",
            GeologicProvince::SubductionOrogen =>
                "Oceanic-continental subduction zone (e.g., Andes)",
            GeologicProvince::AccretionaryOrogen =>
                "Terrane accretion zone (e.g., Western North America)",
            GeologicProvince::ExtensionalOrogen =>
                "Core complex formation (e.g., Basin and Range)",

            // Stage 2.2: Large Igneous Provinces
            GeologicProvince::ContinentalFloodBasalt =>
                "Massive volcanic province (e.g., Deccan Traps)",
            GeologicProvince::OceanicPlateau =>
                "Underwater volcanic plateau (e.g., Ontong Java)",
            GeologicProvince::HotspotTrack =>
                "Volcanic island chain (e.g., Hawaiian Islands)",

            // Stage 2.3: Arc and Basin Systems
            GeologicProvince::VolcanicArc =>
                "Active volcanic arc above subduction (e.g., Aleutians)",
            GeologicProvince::ForearcBasin =>
                "Sedimentary basin between trench and arc (e.g., Great Valley)",
            GeologicProvince::BackarcBasin =>
                "Extensional basin behind arc (e.g., Sea of Japan)",

            // Stage 2.4: Stable Continental Regions
            GeologicProvince::Craton =>
                "Ancient stable core (e.g., Canadian Shield)",
            GeologicProvince::Platform =>
                "Stable region with sedimentary cover (e.g., Russian Platform)",
            GeologicProvince::IntracratonicBasin =>
                "Basin within stable craton (e.g., Michigan Basin)",

            // Stage 2.5: Extensional Zones
            GeologicProvince::ContinentalRift =>
                "Active extensional zone (e.g., East African Rift)",
            GeologicProvince::ExtendedCrust =>
                "Thinned continental crust (e.g., Basin and Range)",

            // Stage 2.6: Oceanic Domains
            GeologicProvince::MidOceanRidge =>
                "Seafloor spreading center (e.g., Mid-Atlantic Ridge)",
            GeologicProvince::AbyssalPlain =>
                "Deep ocean floor (e.g., Pacific abyssal plains)",
            GeologicProvince::OceanTrench =>
                "Deep subduction zone (e.g., Mariana Trench)",
            GeologicProvince::OceanicFractureZone =>
                "Transform fault system (e.g., Romanche Fracture Zone)",
            GeologicProvince::SeamountField =>
                "Underwater volcanic seamounts (e.g., Musicians Seamounts)",
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

    /// Create characteristics for an accretionary orogen
    pub fn accretionary_orogen(convergence_rate: f64, width_km: f64) -> Self {
        let intensity = (convergence_rate / 10.0).min(1.0);

        Self {
            province_type: GeologicProvince::AccretionaryOrogen,
            elevation_intensity: 0.5, // Moderate elevation
            roughness: 0.7,           // Moderately rugged
            width_km,
            intensity,
            convergence_rate,
        }
    }

    /// Create characteristics for an extensional orogen
    pub fn extensional_orogen(extension_rate: f64, width_km: f64) -> Self {
        let intensity = (extension_rate / 10.0).min(1.0);

        Self {
            province_type: GeologicProvince::ExtensionalOrogen,
            elevation_intensity: 0.4, // Moderate with basins
            roughness: 0.6,           // Moderately rough
            width_km,
            intensity,
            convergence_rate: extension_rate, // Reuse field for extension rate
        }
    }

    /// Create characteristics for continental flood basalts
    pub fn continental_flood_basalt(area_km2: f64) -> Self {
        Self {
            province_type: GeologicProvince::ContinentalFloodBasalt,
            elevation_intensity: 0.3, // Relatively flat lava plains
            roughness: 0.2,           // Quite smooth (flood basalts)
            width_km: (area_km2.sqrt()), // Approximate linear dimension
            intensity: 0.8,           // High volcanic activity
            convergence_rate: 0.0,    // Not applicable
        }
    }

    /// Create characteristics for oceanic plateaus
    pub fn oceanic_plateau(area_km2: f64) -> Self {
        Self {
            province_type: GeologicProvince::OceanicPlateau,
            elevation_intensity: -0.3, // Submarine but elevated above abyssal plain
            roughness: 0.3,
            width_km: (area_km2.sqrt()),
            intensity: 0.7,
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for hotspot tracks
    pub fn hotspot_track(length_km: f64) -> Self {
        Self {
            province_type: GeologicProvince::HotspotTrack,
            elevation_intensity: 0.6, // Island chains
            roughness: 0.5,
            width_km: length_km,
            intensity: 0.6,
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for volcanic arcs
    pub fn volcanic_arc(convergence_rate: f64, length_km: f64) -> Self {
        let intensity = (convergence_rate / 10.0).min(1.0);

        Self {
            province_type: GeologicProvince::VolcanicArc,
            elevation_intensity: 0.8, // High volcanic peaks
            roughness: 0.9,           // Very rough volcanic terrain
            width_km: length_km,
            intensity,
            convergence_rate,
        }
    }

    /// Create characteristics for forearc basins
    pub fn forearc_basin(width_km: f64) -> Self {
        Self {
            province_type: GeologicProvince::ForearcBasin,
            elevation_intensity: -0.2, // Below sea level typically
            roughness: 0.1,            // Smooth sedimentary basin
            width_km,
            intensity: 0.3,
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for backarc basins
    pub fn backarc_basin(width_km: f64) -> Self {
        Self {
            province_type: GeologicProvince::BackarcBasin,
            elevation_intensity: -0.4, // Oceanic or near-oceanic depths
            roughness: 0.2,
            width_km,
            intensity: 0.5, // Some spreading/extension
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for cratons
    pub fn craton(area_km2: f64) -> Self {
        Self {
            province_type: GeologicProvince::Craton,
            elevation_intensity: 0.2, // Low, stable elevation
            roughness: 0.1,           // Very smooth, eroded
            width_km: (area_km2.sqrt()),
            intensity: 0.0, // Tectonically inactive
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for platforms
    pub fn platform(area_km2: f64) -> Self {
        Self {
            province_type: GeologicProvince::Platform,
            elevation_intensity: 0.15, // Slightly lower than cratons
            roughness: 0.05,           // Very flat (sedimentary cover)
            width_km: (area_km2.sqrt()),
            intensity: 0.0,
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for intracratonic basins
    pub fn intracratonic_basin(area_km2: f64) -> Self {
        Self {
            province_type: GeologicProvince::IntracratonicBasin,
            elevation_intensity: 0.0, // At or below sea level
            roughness: 0.05,          // Very flat
            width_km: (area_km2.sqrt()),
            intensity: 0.0,
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for continental rifts
    pub fn continental_rift(extension_rate: f64, length_km: f64) -> Self {
        let intensity = (extension_rate / 10.0).min(1.0);

        Self {
            province_type: GeologicProvince::ContinentalRift,
            elevation_intensity: 0.3, // Uplifted rift shoulders, down-dropped basin
            roughness: 0.7,           // Rough from faulting
            width_km: length_km,
            intensity,
            convergence_rate: -extension_rate, // Negative for extension
        }
    }

    /// Create characteristics for extended crust
    pub fn extended_crust(area_km2: f64) -> Self {
        Self {
            province_type: GeologicProvince::ExtendedCrust,
            elevation_intensity: 0.1, // Lowered by extension
            roughness: 0.5,           // Moderately rough from faulting
            width_km: (area_km2.sqrt()),
            intensity: 0.3,
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for mid-ocean ridges
    pub fn mid_ocean_ridge(spreading_rate: f64, length_km: f64) -> Self {
        let intensity = (spreading_rate / 10.0).min(1.0);

        Self {
            province_type: GeologicProvince::MidOceanRidge,
            elevation_intensity: -0.5, // Submarine but elevated
            roughness: 0.7,            // Rugged from volcanic activity
            width_km: length_km,
            intensity,
            convergence_rate: -spreading_rate, // Negative for divergence
        }
    }

    /// Create characteristics for abyssal plains
    pub fn abyssal_plain(area_km2: f64) -> Self {
        Self {
            province_type: GeologicProvince::AbyssalPlain,
            elevation_intensity: -1.0, // Deepest ocean floors
            roughness: 0.05,           // Very flat (sediment covered)
            width_km: (area_km2.sqrt()),
            intensity: 0.0,
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for ocean trenches
    pub fn ocean_trench(convergence_rate: f64, length_km: f64) -> Self {
        let intensity = (convergence_rate / 10.0).min(1.0);

        Self {
            province_type: GeologicProvince::OceanTrench,
            elevation_intensity: -1.2, // Deepest parts of ocean
            roughness: 0.4,            // Steep walls
            width_km: length_km,
            intensity,
            convergence_rate,
        }
    }

    /// Create characteristics for oceanic fracture zones
    pub fn oceanic_fracture_zone(length_km: f64) -> Self {
        Self {
            province_type: GeologicProvince::OceanicFractureZone,
            elevation_intensity: -0.7, // Elevated above abyssal plain
            roughness: 0.6,            // Rough from faulting
            width_km: length_km,
            intensity: 0.3,
            convergence_rate: 0.0,
        }
    }

    /// Create characteristics for seamount fields
    pub fn seamount_field(area_km2: f64) -> Self {
        Self {
            province_type: GeologicProvince::SeamountField,
            elevation_intensity: -0.2, // Rises from ocean floor but submerged
            roughness: 0.8,            // Very rough volcanic terrain
            width_km: (area_km2.sqrt()),
            intensity: 0.5, // Moderate volcanic activity
            convergence_rate: 0.0,
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
