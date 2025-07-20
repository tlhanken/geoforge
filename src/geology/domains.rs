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
    pub fn generate(&self, plate_map: &[u16]) -> Result<Vec<DomainInfo>, GeologyError> {
        // TODO: Implement domain generation logic
        // This will analyze plate boundaries, interactions, and ages to determine domains
        Ok(vec![])
    }
    
    /// Generate realistic geology map based on tectonic plate analysis
    pub fn generate_geology_map(&self, plate_map: &crate::map::terrain::PlateMap) -> crate::map::terrain::GeologyMap {
        let mut geology_map = crate::map::terrain::TerrainMap::new(
            self.width, 
            self.height, 
            GeologicDomain::Oceanic(OceanicType::AbyssalPlain)
        );
        
        // First pass: Detect plate boundaries and classify pixels
        let boundary_info = self.detect_plate_boundaries(plate_map);
        
        // Generate domains with progress indication
        println!("Classifying geological domains...");
        let total_pixels = self.width * self.height;
        let progress_interval = total_pixels / 10; // 10% intervals
        let mut processed = 0;
        
        for y in 0..self.height {
            for x in 0..self.width {
                let domain = self.classify_pixel_domain(x, y, plate_map, &boundary_info);
                geology_map.set(x, y, domain);
                
                processed += 1;
                if processed % progress_interval == 0 {
                    let progress = (processed * 100) / total_pixels;
                    println!("  Progress: {}%", progress);
                }
            }
        }
        
        geology_map
    }
    
    /// Detect plate boundaries and classify them
    fn detect_plate_boundaries(&self, plate_map: &crate::map::terrain::PlateMap) -> BoundaryInfo {
        let mut boundary_pixels = std::collections::HashSet::new();
        let mut boundary_types = std::collections::HashMap::new();
        
        for y in 0..self.height {
            for x in 0..self.width {
                let current_plate = plate_map.get(x, y).copied().unwrap_or(0);
                let neighbors = plate_map.get_neighbors(x, y);
                
                // Check if this pixel is on a plate boundary
                for (nx, ny) in neighbors {
                    if let Some(neighbor_plate) = plate_map.get(nx, ny) {
                        if *neighbor_plate != current_plate {
                            boundary_pixels.insert((x, y));
                            
                            // Classify boundary type based on plate motion
                            let boundary_type = self.classify_boundary_interaction(
                                current_plate, *neighbor_plate
                            );
                            boundary_types.insert((x, y), boundary_type);
                            break;
                        }
                    }
                }
            }
        }
        
        BoundaryInfo {
            boundary_pixels,
            boundary_types,
        }
    }
    
    /// Classify the geological domain for a specific pixel
    fn classify_pixel_domain(
        &self, 
        x: usize, 
        y: usize, 
        plate_map: &crate::map::terrain::PlateMap,
        boundary_info: &BoundaryInfo
    ) -> GeologicDomain {
        let plate_id = plate_map.get(x, y).copied().unwrap_or(0);
        
        // Distance to nearest boundary (simplified)
        let distance_to_boundary = self.distance_to_nearest_boundary(x, y, &boundary_info.boundary_pixels);
        
        // Classify based on location relative to boundaries
        if boundary_info.boundary_pixels.contains(&(x, y)) {
            // This pixel is on a plate boundary
            if let Some(boundary_type) = boundary_info.boundary_types.get(&(x, y)) {
                self.classify_boundary_domain(*boundary_type, plate_id)
            } else {
                GeologicDomain::Oceanic(OceanicType::FractureZone)
            }
        } else if distance_to_boundary < 5 {
            // Near boundary - create transitional zones
            self.classify_near_boundary_domain(x, y, plate_id, distance_to_boundary, boundary_info)
        } else {
            // Interior of plate - stable regions
            self.classify_plate_interior_domain(plate_id, distance_to_boundary)
        }
    }
    
    /// Classify boundary interaction type between two plates
    fn classify_boundary_interaction(&self, plate_a: u16, plate_b: u16) -> PlateInteraction {
        // Find the seeds for these plates
        let seed_a = self.plate_seeds.iter().find(|s| s.id == plate_a);
        let seed_b = self.plate_seeds.iter().find(|s| s.id == plate_b);
        
        match (seed_a, seed_b) {
            (Some(a), Some(b)) => {
                // Calculate relative motion vectors
                let motion_diff = (a.motion_direction - b.motion_direction).abs();
                let motion_diff = if motion_diff > 180.0 { 360.0 - motion_diff } else { motion_diff };
                
                // Classify based on motion direction difference
                if motion_diff < 45.0 {
                    // Similar directions - transform boundary
                    PlateInteraction::Transform
                } else if motion_diff > 135.0 {
                    // Opposite directions - convergent boundary
                    PlateInteraction::Convergent
                } else {
                    // Perpendicular - divergent boundary
                    PlateInteraction::Divergent
                }
            }
            _ => PlateInteraction::Transform // Default fallback
        }
    }
    
    /// Classify domain at plate boundary based on interaction type
    fn classify_boundary_domain(&self, interaction: PlateInteraction, plate_id: u16) -> GeologicDomain {
        // Get plate characteristics
        let plate_age = self.estimate_plate_age(plate_id);
        let plate_size = self.estimate_plate_size(plate_id);
        
        match interaction {
            PlateInteraction::Convergent => {
                if plate_age > 1500.0 && plate_size > 0.05 {
                    // Large, old plates - continental collision
                    GeologicDomain::Orogenic(OrogenicType::CollisionOrogeny)
                } else if plate_age < 200.0 {
                    // Young oceanic plate - subduction
                    GeologicDomain::ArcSystem(ArcType::VolcanicArc)
                } else {
                    // Mixed - subduction orogeny
                    GeologicDomain::Orogenic(OrogenicType::SubductionOrogeny)
                }
            }
            PlateInteraction::Divergent => {
                if plate_age < 100.0 {
                    // Young spreading center
                    GeologicDomain::Oceanic(OceanicType::MidOceanRidge)
                } else {
                    // Continental rifting
                    GeologicDomain::Extensional(ExtensionalType::ContinentalRift)
                }
            }
            PlateInteraction::Transform => {
                if plate_age < 500.0 {
                    // Oceanic transform
                    GeologicDomain::Oceanic(OceanicType::FractureZone)
                } else {
                    // Continental transform
                    GeologicDomain::Extensional(ExtensionalType::ExtendedCrust)
                }
            }
        }
    }
    
    /// Classify domains near plate boundaries (simplified for performance)
    fn classify_near_boundary_domain(
        &self, 
        _x: usize, 
        _y: usize, 
        plate_id: u16, 
        distance: usize,
        _boundary_info: &BoundaryInfo
    ) -> GeologicDomain {
        let plate_age = self.estimate_plate_age(plate_id);
        
        // Simplified classification based on distance and plate age
        match distance {
            1..=2 => {
                if plate_age < 200.0 {
                    GeologicDomain::ArcSystem(ArcType::VolcanicArc)
                } else {
                    GeologicDomain::Orogenic(OrogenicType::SubductionOrogeny)
                }
            }
            3..=4 => {
                if plate_age < 500.0 {
                    GeologicDomain::ArcSystem(ArcType::BackarcBasin)
                } else {
                    GeologicDomain::Extensional(ExtensionalType::ExtendedCrust)
                }
            }
            _ => {
                GeologicDomain::Extensional(ExtensionalType::TransitionalCrust)
            }
        }
    }
    
    /// Classify stable plate interior domains
    fn classify_plate_interior_domain(&self, plate_id: u16, distance_to_boundary: usize) -> GeologicDomain {
        let plate_age = self.estimate_plate_age(plate_id);
        let plate_size = self.estimate_plate_size(plate_id);
        
        if plate_age > 2000.0 && plate_size > 0.03 {
            // Ancient, large plates - cratons
            if distance_to_boundary > 20 {
                GeologicDomain::StableContinental(StableType::Craton)
            } else {
                GeologicDomain::StableContinental(StableType::Platform)
            }
        } else if plate_age > 500.0 {
            // Intermediate age - platforms and basins
            if distance_to_boundary > 15 {
                GeologicDomain::StableContinental(StableType::IntracratonicBasin)
            } else {
                GeologicDomain::StableContinental(StableType::Platform)
            }
        } else if plate_age < 100.0 {
            // Young oceanic plate
            if distance_to_boundary > 10 {
                GeologicDomain::Oceanic(OceanicType::AbyssalPlain)
            } else {
                GeologicDomain::IgneousProvince(IgneousType::HotspotTrack)
            }
        } else {
            // Default oceanic
            GeologicDomain::Oceanic(OceanicType::AbyssalPlain)
        }
    }
    
    /// Estimate plate age based on its size and position
    fn estimate_plate_age(&self, plate_id: u16) -> f64 {
        // Larger plates tend to be older (simplified model)
        let plate_size = self.estimate_plate_size(plate_id);
        
        if plate_size > 0.1 {
            2500.0 + (plate_size * 1000.0) // Very old for large plates
        } else if plate_size > 0.05 {
            1000.0 + (plate_size * 500.0)  // Old for medium plates
        } else {
            100.0 + (plate_size * 200.0)   // Young for small plates
        }
    }
    
    /// Estimate relative plate size (0.0 to 1.0)
    fn estimate_plate_size(&self, plate_id: u16) -> f64 {
        // Use a pseudo-random but deterministic calculation
        // In reality, this would come from actual plate statistics
        let hash = ((plate_id as u64 * 2654435761u64) % 1000) as f64 / 1000.0;
        
        // Create realistic size distribution using power law
        if hash < 0.1 {
            0.1 + hash * 0.9  // 10% chance of large plates (0.1-1.0)
        } else if hash < 0.3 {
            0.05 + hash * 0.1 // 20% chance of medium plates (0.05-0.15)
        } else {
            0.001 + hash * 0.05 // 70% chance of small plates (0.001-0.051)
        }
    }
    
    /// Calculate distance to nearest boundary (optimized with early termination)
    fn distance_to_nearest_boundary(&self, x: usize, y: usize, boundaries: &std::collections::HashSet<(usize, usize)>) -> usize {
        let mut min_distance = 100; // Start with large distance
        
        // Only check nearby boundaries for efficiency
        for radius in 1..=10 {
            for dy in -(radius as i32)..=(radius as i32) {
                for dx in -(radius as i32)..=(radius as i32) {
                    let check_x = (x as i32 + dx) as usize;
                    let check_y = (y as i32 + dy) as usize;
                    
                    if boundaries.contains(&(check_x, check_y)) {
                        let distance = (dx.abs() + dy.abs()) as usize;
                        if distance < min_distance {
                            min_distance = distance;
                        }
                        if min_distance == 1 {
                            return 1; // Early exit for immediate neighbors
                        }
                    }
                }
            }
            if min_distance <= radius as usize {
                break; // Found a boundary within current radius
            }
        }
        
        min_distance
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

/// Information about plate boundaries and their characteristics
#[derive(Debug, Clone)]
struct BoundaryInfo {
    boundary_pixels: std::collections::HashSet<(usize, usize)>,
    boundary_types: std::collections::HashMap<(usize, usize), PlateInteraction>,
}

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