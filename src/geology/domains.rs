//! Geologic domain classification and modeling
//! 
//! This module implements the Stage 2 geological provinces system that builds on tectonic plates
//! to create realistic geological domains including orogenic belts, large igneous provinces,
//! arc systems, stable continental regions, extensional zones, and oceanic domains.

use crate::map::spherical::SphericalPoint;
use crate::tectonics::{PlateSeed, PlateInteraction};

/// Major geologic domain types based on tectonic setting
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeologicDomain {
    /// Orogenic Belts - Mountain-building zones
    Orogenic(OrogenicType),
    /// Large Igneous Provinces - Major volcanic regions
    IgneousProvince(IgneousType),
    /// Arc and Basin Systems - Subduction-related features
    ArcSystem(ArcType),
    /// Stable Continental Regions - Ancient, stable areas
    StableContinental(StableType),
    /// Extensional Zones - Areas of crustal extension
    Extensional(ExtensionalType),
    /// Oceanic Domains - Deep ocean features
    Oceanic(OceanicType),
}

/// Types of orogenic (mountain-building) belts
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrogenicType {
    /// Continental plate collision (Himalayas-style)
    CollisionOrogeny,
    /// Oceanic-continental convergence (Andes-style)
    SubductionOrogeny,
    /// Terrane accretion and exotic block collision
    AccretionaryOrogeny,
    /// Core complex formation in rifting zones
    ExtensionalOrogeny,
}

/// Types of large igneous provinces
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IgneousType {
    /// Massive continental volcanic provinces (Deccan Traps-style)
    ContinentalFloodBasalt,
    /// Underwater volcanic provinces
    OceanicPlateau,
    /// Volcanic island chains and seamount trails
    HotspotTrack,
    /// Radiating intrusive networks
    DykeSwarm,
}

/// Types of arc and basin systems
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArcType {
    /// Active subduction zone volcanism
    VolcanicArc,
    /// Sedimentary basins between trench and arc
    ForearcBasin,
    /// Extensional basins behind volcanic arcs
    BackarcBasin,
    /// Spreading centers in backarc regions
    BackarcRidge,
}

/// Types of stable continental regions
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StableType {
    /// Ancient, stable continental cores (>1.5 Ga)
    Craton,
    /// Stable cratonic areas with thin sedimentary cover
    Platform,
    /// Subsided regions within stable cratons
    IntracratonicBasin,
}

/// Types of extensional zones
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtensionalType {
    /// Active extension and normal faulting
    ContinentalRift,
    /// Thinned continental crust from extension
    ExtendedCrust,
    /// Continent-ocean boundary zones
    TransitionalCrust,
}

/// Types of oceanic domains
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OceanicType {
    /// Active seafloor spreading centers
    MidOceanRidge,
    /// Deep oceanic basins with sediment cover
    AbyssalPlain,
    /// Transform fault systems
    FractureZone,
    /// Subduction zone depocenters
    DeepTrench,
}

/// Information about a geologic domain at a specific location
#[derive(Debug, Clone)]
pub struct DomainInfo {
    pub domain_type: GeologicDomain,
    pub age_ma: f64,           // Age in millions of years
    pub elevation_m: f64,      // Elevation in meters (negative for below sea level)
    pub crustal_thickness_km: f64,  // Crustal thickness in kilometers
    pub associated_plates: Vec<u16>, // IDs of associated tectonic plates
}

impl DomainInfo {
    /// Create a new domain info structure
    pub fn new(
        domain_type: GeologicDomain,
        age_ma: f64,
        elevation_m: f64,
        crustal_thickness_km: f64,
        associated_plates: Vec<u16>,
    ) -> Self {
        Self {
            domain_type,
            age_ma,
            elevation_m,
            crustal_thickness_km,
            associated_plates,
        }
    }
    
    /// Check if this domain is associated with a specific plate
    pub fn involves_plate(&self, plate_id: u16) -> bool {
        self.associated_plates.contains(&plate_id)
    }
    
    /// Get typical crustal density for this domain type (kg/m³)
    pub fn typical_crustal_density(&self) -> f64 {
        match self.domain_type {
            GeologicDomain::Oceanic(_) => 2900.0,  // Oceanic basalt
            GeologicDomain::StableContinental(_) => 2700.0,  // Continental granite
            GeologicDomain::Orogenic(_) => 2750.0,  // Mixed igneous/metamorphic
            GeologicDomain::IgneousProvince(_) => 2850.0,  // Mafic volcanics
            GeologicDomain::ArcSystem(_) => 2650.0,  // Andesitic composition
            GeologicDomain::Extensional(_) => 2600.0,  // Thinned continental
        }
    }
}

/// Generator for geologic domains based on tectonic plate configuration
pub struct GeologicDomainGenerator {
    width: usize,
    height: usize,
    plate_seeds: Vec<PlateSeed>,
}

impl GeologicDomainGenerator {
    /// Create a new geologic domain generator
    pub fn new(width: usize, height: usize, plate_seeds: Vec<PlateSeed>) -> Self {
        Self {
            width,
            height,
            plate_seeds,
        }
    }
    
    /// Generate geologic domains based on plate configuration
    pub fn generate(&self, _plate_map: &[u16]) -> Result<Vec<DomainInfo>, GeologyError> {
        // TODO: Implement domain generation logic
        // This will analyze plate boundaries, interactions, and ages to determine domains
        Ok(vec![])
    }
    
    /// Create a simple geology map for testing purposes
    /// This assigns random but realistic domain types to demonstrate the system
    pub fn generate_simple_map(&self) -> crate::map::terrain::GeologyMap {
        use rand::prelude::*;
        
        let mut rng = StdRng::seed_from_u64(42); // Consistent results
        let mut geology_map = crate::map::terrain::TerrainMap::new(
            self.width, 
            self.height, 
            GeologicDomain::Oceanic(OceanicType::AbyssalPlain)
        );
        
        // Create a simple pattern for demonstration
        // In reality, this would be based on tectonic plate analysis
        for y in 0..self.height {
            for x in 0..self.width {
                let domain = match (x + y) % 22 {
                    0 => GeologicDomain::Orogenic(OrogenicType::CollisionOrogeny),
                    1 => GeologicDomain::Orogenic(OrogenicType::SubductionOrogeny),
                    2 => GeologicDomain::Orogenic(OrogenicType::AccretionaryOrogeny),
                    3 => GeologicDomain::Orogenic(OrogenicType::ExtensionalOrogeny),
                    4 => GeologicDomain::IgneousProvince(IgneousType::ContinentalFloodBasalt),
                    5 => GeologicDomain::IgneousProvince(IgneousType::OceanicPlateau),
                    6 => GeologicDomain::IgneousProvince(IgneousType::HotspotTrack),
                    7 => GeologicDomain::IgneousProvince(IgneousType::DykeSwarm),
                    8 => GeologicDomain::ArcSystem(ArcType::VolcanicArc),
                    9 => GeologicDomain::ArcSystem(ArcType::ForearcBasin),
                    10 => GeologicDomain::ArcSystem(ArcType::BackarcBasin),
                    11 => GeologicDomain::ArcSystem(ArcType::BackarcRidge),
                    12 => GeologicDomain::StableContinental(StableType::Craton),
                    13 => GeologicDomain::StableContinental(StableType::Platform),
                    14 => GeologicDomain::StableContinental(StableType::IntracratonicBasin),
                    15 => GeologicDomain::Extensional(ExtensionalType::ContinentalRift),
                    16 => GeologicDomain::Extensional(ExtensionalType::ExtendedCrust),
                    17 => GeologicDomain::Extensional(ExtensionalType::TransitionalCrust),
                    18 => GeologicDomain::Oceanic(OceanicType::MidOceanRidge),
                    19 => GeologicDomain::Oceanic(OceanicType::AbyssalPlain),
                    20 => GeologicDomain::Oceanic(OceanicType::FractureZone),
                    21 => GeologicDomain::Oceanic(OceanicType::DeepTrench),
                    _ => GeologicDomain::Oceanic(OceanicType::AbyssalPlain),
                };
                
                geology_map.set(x, y, domain);
            }
        }
        
        geology_map
    }
    
    /// Determine domain type at plate boundary based on interaction type
    fn classify_boundary_domain(
        &self,
        interaction: PlateInteraction,
        plate_a_age: f64,
        plate_b_age: f64,
    ) -> GeologicDomain {
        match interaction {
            PlateInteraction::Convergent => {
                if plate_a_age > 1500.0 && plate_b_age > 1500.0 {
                    // Old continental collision
                    GeologicDomain::Orogenic(OrogenicType::CollisionOrogeny)
                } else {
                    // Oceanic-continental subduction
                    GeologicDomain::Orogenic(OrogenicType::SubductionOrogeny)
                }
            }
            PlateInteraction::Divergent => {
                GeologicDomain::Oceanic(OceanicType::MidOceanRidge)
            }
            PlateInteraction::Transform => {
                GeologicDomain::Oceanic(OceanicType::FractureZone)
            }
        }
    }
}

/// Error types for geology operations
#[derive(Debug)]
pub enum GeologyError {
    InvalidPlateConfiguration,
    DomainGenerationFailed(String),
}

impl std::fmt::Display for GeologyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GeologyError::InvalidPlateConfiguration => {
                write!(f, "Invalid plate configuration for domain generation")
            }
            GeologyError::DomainGenerationFailed(msg) => {
                write!(f, "Domain generation failed: {}", msg)
            }
        }
    }
}

impl std::error::Error for GeologyError {}

impl GeologicDomain {
    /// Get the RGB color for this geologic domain (compatible with image::Rgb)
    pub fn color(&self) -> [u8; 3] {
        match self {
            // Orogenic Belts - Warm browns and earth tones for mountain-building zones
            GeologicDomain::Orogenic(OrogenicType::CollisionOrogeny) => [139, 69, 19],    // Saddle Brown - Continental collision
            GeologicDomain::Orogenic(OrogenicType::SubductionOrogeny) => [160, 82, 45],   // Sienna - Oceanic-continental
            GeologicDomain::Orogenic(OrogenicType::AccretionaryOrogeny) => [210, 180, 140], // Tan - Terrane accretion
            GeologicDomain::Orogenic(OrogenicType::ExtensionalOrogeny) => [222, 184, 135], // Burlywood - Core complexes
            
            // Large Igneous Provinces - Red/orange volcanic colors
            GeologicDomain::IgneousProvince(IgneousType::ContinentalFloodBasalt) => [178, 34, 34],  // Fire Brick - Flood basalts
            GeologicDomain::IgneousProvince(IgneousType::OceanicPlateau) => [220, 20, 60],         // Crimson - Oceanic plateaus
            GeologicDomain::IgneousProvince(IgneousType::HotspotTrack) => [255, 69, 0],            // Red Orange - Hotspot tracks
            GeologicDomain::IgneousProvince(IgneousType::DykeSwarm) => [255, 99, 71],              // Tomato - Dyke swarms
            
            // Arc Systems - Purple/magenta for subduction-related features
            GeologicDomain::ArcSystem(ArcType::VolcanicArc) => [148, 0, 211],          // Dark Violet - Volcanic arcs
            GeologicDomain::ArcSystem(ArcType::ForearcBasin) => [186, 85, 211],        // Medium Orchid - Forearc basins
            GeologicDomain::ArcSystem(ArcType::BackarcBasin) => [221, 160, 221],       // Plum - Backarc basins
            GeologicDomain::ArcSystem(ArcType::BackarcRidge) => [238, 130, 238],       // Violet - Backarc ridges
            
            // Stable Continental - Green tones for ancient, stable areas
            GeologicDomain::StableContinental(StableType::Craton) => [34, 139, 34],    // Forest Green - Ancient cratons
            GeologicDomain::StableContinental(StableType::Platform) => [50, 205, 50],  // Lime Green - Platforms
            GeologicDomain::StableContinental(StableType::IntracratonicBasin) => [144, 238, 144], // Light Green - Intracratonic basins
            
            // Extensional Zones - Yellow/orange for rifting and extension
            GeologicDomain::Extensional(ExtensionalType::ContinentalRift) => [255, 215, 0],     // Gold - Continental rifts
            GeologicDomain::Extensional(ExtensionalType::ExtendedCrust) => [255, 228, 181],    // Moccasin - Extended crust
            GeologicDomain::Extensional(ExtensionalType::TransitionalCrust) => [255, 248, 220], // Cornsilk - Transitional crust
            
            // Oceanic Domains - Blue tones for deep ocean features
            GeologicDomain::Oceanic(OceanicType::MidOceanRidge) => [0, 0, 139],        // Dark Blue - Mid-ocean ridges
            GeologicDomain::Oceanic(OceanicType::AbyssalPlain) => [25, 25, 112],       // Midnight Blue - Abyssal plains
            GeologicDomain::Oceanic(OceanicType::FractureZone) => [72, 61, 139],       // Dark Slate Blue - Fracture zones
            GeologicDomain::Oceanic(OceanicType::DeepTrench) => [0, 0, 0],             // Black - Deep trenches
        }
    }
    
    /// Get the base color for the domain type (before subtype variation)
    pub fn base_color(&self) -> [u8; 3] {
        match self {
            GeologicDomain::Orogenic(_) => [139, 69, 19],        // Brown base for orogenic
            GeologicDomain::IgneousProvince(_) => [178, 34, 34],  // Red base for igneous
            GeologicDomain::ArcSystem(_) => [148, 0, 211],        // Purple base for arcs
            GeologicDomain::StableContinental(_) => [34, 139, 34], // Green base for stable
            GeologicDomain::Extensional(_) => [255, 215, 0],      // Yellow base for extensional
            GeologicDomain::Oceanic(_) => [0, 0, 139],            // Blue base for oceanic
        }
    }
    
    /// Get a human-readable name for this domain
    pub fn name(&self) -> &'static str {
        match self {
            GeologicDomain::Orogenic(OrogenicType::CollisionOrogeny) => "Collision Orogeny",
            GeologicDomain::Orogenic(OrogenicType::SubductionOrogeny) => "Subduction Orogeny", 
            GeologicDomain::Orogenic(OrogenicType::AccretionaryOrogeny) => "Accretionary Orogeny",
            GeologicDomain::Orogenic(OrogenicType::ExtensionalOrogeny) => "Extensional Orogeny",
            
            GeologicDomain::IgneousProvince(IgneousType::ContinentalFloodBasalt) => "Continental Flood Basalt",
            GeologicDomain::IgneousProvince(IgneousType::OceanicPlateau) => "Oceanic Plateau",
            GeologicDomain::IgneousProvince(IgneousType::HotspotTrack) => "Hotspot Track",
            GeologicDomain::IgneousProvince(IgneousType::DykeSwarm) => "Dyke Swarm",
            
            GeologicDomain::ArcSystem(ArcType::VolcanicArc) => "Volcanic Arc",
            GeologicDomain::ArcSystem(ArcType::ForearcBasin) => "Forearc Basin", 
            GeologicDomain::ArcSystem(ArcType::BackarcBasin) => "Backarc Basin",
            GeologicDomain::ArcSystem(ArcType::BackarcRidge) => "Backarc Ridge",
            
            GeologicDomain::StableContinental(StableType::Craton) => "Craton",
            GeologicDomain::StableContinental(StableType::Platform) => "Platform",
            GeologicDomain::StableContinental(StableType::IntracratonicBasin) => "Intracratonic Basin",
            
            GeologicDomain::Extensional(ExtensionalType::ContinentalRift) => "Continental Rift",
            GeologicDomain::Extensional(ExtensionalType::ExtendedCrust) => "Extended Crust",
            GeologicDomain::Extensional(ExtensionalType::TransitionalCrust) => "Transitional Crust",
            
            GeologicDomain::Oceanic(OceanicType::MidOceanRidge) => "Mid-Ocean Ridge",
            GeologicDomain::Oceanic(OceanicType::AbyssalPlain) => "Abyssal Plain",
            GeologicDomain::Oceanic(OceanicType::FractureZone) => "Fracture Zone",
            GeologicDomain::Oceanic(OceanicType::DeepTrench) => "Deep Trench",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_domain_info_creation() {
        let domain = DomainInfo::new(
            GeologicDomain::Orogenic(OrogenicType::CollisionOrogeny),
            50.0,
            2000.0,
            45.0,
            vec![1, 2],
        );
        
        assert!(domain.involves_plate(1));
        assert!(domain.involves_plate(2));
        assert!(!domain.involves_plate(3));
        assert_eq!(domain.age_ma, 50.0);
    }
    
    #[test]
    fn test_crustal_density() {
        let oceanic = DomainInfo::new(
            GeologicDomain::Oceanic(OceanicType::AbyssalPlain),
            100.0, -4000.0, 7.0, vec![1]
        );
        let continental = DomainInfo::new(
            GeologicDomain::StableContinental(StableType::Craton),
            2500.0, 500.0, 35.0, vec![2]
        );
        
        assert!(oceanic.typical_crustal_density() > continental.typical_crustal_density());
    }
}