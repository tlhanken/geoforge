//! Geologic domain classification and modeling
//! 
//! This module implements the Stage 2 geological provinces system that builds on tectonic plates
//! to create realistic geological domains including orogenic belts, large igneous provinces,
//! arc systems, stable continental regions, extensional zones, and oceanic domains.

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
    
    /// Generate realistic geology map based on tectonic plate analysis (optimized)
    pub fn generate_geology_map(&self, plate_map: &crate::map::terrain::PlateMap) -> crate::map::terrain::GeologyMap {
        println!("Generating optimized geological domains...");
        
        // Step 1: Pre-compute plate statistics (O(plates) instead of O(pixels))
        let plate_stats = self.compute_plate_statistics(plate_map);
        println!("  Computed statistics for {} plates", plate_stats.len());
        
        // Step 2: Build efficient distance field for boundaries (O(pixels) single pass)
        let distance_field = self.build_distance_field(plate_map);
        println!("  Built boundary distance field");
        
        // Step 3: Pre-compute boundary interactions (O(boundaries) instead of O(pixels))
        let boundary_interactions = self.compute_boundary_interactions(plate_map, &plate_stats);
        println!("  Computed boundary interactions");
        
        // Step 4: Generate domains using cached data (O(pixels) single pass)
        let mut geology_map = crate::map::terrain::TerrainMap::new(
            self.width, 
            self.height, 
            GeologicDomain::Oceanic(OceanicType::AbyssalPlain)
        );
        
        for y in 0..self.height {
            for x in 0..self.width {
                let domain = self.classify_pixel_optimized(
                    x, y, plate_map, &distance_field, &plate_stats, &boundary_interactions
                );
                geology_map.set(x, y, domain);
            }
        }
        
        println!("  Geological domain generation complete");
        geology_map
    }
    
    
    /// Pre-compute plate statistics for efficient domain classification
    fn compute_plate_statistics(&self, plate_map: &crate::map::terrain::PlateMap) -> std::collections::HashMap<u16, PlateStatistics> {
        let mut stats = std::collections::HashMap::new();
        let mut plate_pixels: std::collections::HashMap<u16, Vec<(usize, usize)>> = std::collections::HashMap::new();
        
        // Count pixels per plate
        for y in 0..self.height {
            for x in 0..self.width {
                if let Some(&plate_id) = plate_map.get(x, y) {
                    plate_pixels.entry(plate_id).or_insert_with(Vec::new).push((x, y));
                }
            }
        }
        
        // Calculate statistics for each plate
        for (plate_id, pixels) in plate_pixels {
            let size_ratio = pixels.len() as f64 / (self.width * self.height) as f64;
            let age = self.calculate_realistic_plate_age(plate_id, size_ratio);
            let motion = self.get_plate_motion(plate_id);
            
            stats.insert(plate_id, PlateStatistics {
                size_ratio,
                age_ma: age,
                motion_direction: motion.0,
                motion_speed: motion.1,
                center: self.calculate_plate_center(&pixels),
                is_oceanic: age < 800.0 || size_ratio < 0.03, // More oceanic classification
            });
        }
        
        stats
    }
    
    /// Build distance field using efficient flood-fill algorithm with expanded ranges
    fn build_distance_field(&self, plate_map: &crate::map::terrain::PlateMap) -> Vec<Vec<u8>> {
        let mut distance_field = vec![vec![255u8; self.width]; self.height];
        let mut queue = std::collections::VecDeque::new();
        
        // Find all boundary pixels and set them to distance 0
        for y in 0..self.height {
            for x in 0..self.width {
                if self.is_boundary_pixel(x, y, plate_map) {
                    distance_field[y][x] = 0;
                    queue.push_back((x, y, 0u8));
                }
            }
        }
        
        // Flood-fill to compute distances (much larger range for realistic provinces)
        while let Some((x, y, dist)) = queue.pop_front() {
            if dist >= 80 { continue; } // Expanded from 20 to 80 for larger provinces
            
            for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let nx = (x as i32 + dx) as usize;
                let ny = (y as i32 + dy) as usize;
                
                if nx < self.width && ny < self.height {
                    let new_dist = dist + 1;
                    if distance_field[ny][nx] > new_dist {
                        distance_field[ny][nx] = new_dist;
                        queue.push_back((nx, ny, new_dist));
                    }
                }
            }
        }
        
        distance_field
    }
    
    /// Check if a pixel is on a plate boundary
    fn is_boundary_pixel(&self, x: usize, y: usize, plate_map: &crate::map::terrain::PlateMap) -> bool {
        let current_plate = plate_map.get(x, y).copied().unwrap_or(0);
        
        for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
            let nx = (x as i32 + dx) as usize;
            let ny = (y as i32 + dy) as usize;
            
            if nx < self.width && ny < self.height {
                if let Some(&neighbor_plate) = plate_map.get(nx, ny) {
                    if neighbor_plate != current_plate {
                        return true;
                    }
                }
            }
        }
        false
    }
    
    /// Pre-compute boundary interactions for efficient lookup
    fn compute_boundary_interactions(
        &self, 
        plate_map: &crate::map::terrain::PlateMap,
        plate_stats: &std::collections::HashMap<u16, PlateStatistics>
    ) -> std::collections::HashMap<(u16, u16), PlateInteraction> {
        let mut interactions = std::collections::HashMap::new();
        
        for y in 0..self.height {
            for x in 0..self.width {
                if self.is_boundary_pixel(x, y, plate_map) {
                    let current_plate = plate_map.get(x, y).copied().unwrap_or(0);
                    
                    for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                        let nx = (x as i32 + dx) as usize;
                        let ny = (y as i32 + dy) as usize;
                        
                        if nx < self.width && ny < self.height {
                            if let Some(&neighbor_plate) = plate_map.get(nx, ny) {
                                if neighbor_plate != current_plate {
                                    let key = if current_plate < neighbor_plate {
                                        (current_plate, neighbor_plate)
                                    } else {
                                        (neighbor_plate, current_plate)
                                    };
                                    
                                    if !interactions.contains_key(&key) {
                                        let interaction = self.determine_plate_interaction(
                                            current_plate, neighbor_plate, plate_stats
                                        );
                                        interactions.insert(key, interaction);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        interactions
    }
    
    /// Classify pixels into large coherent geological regions
    fn classify_pixel_optimized(
        &self,
        x: usize,
        y: usize,
        plate_map: &crate::map::terrain::PlateMap,
        distance_field: &[Vec<u8>],
        plate_stats: &std::collections::HashMap<u16, PlateStatistics>,
        boundary_interactions: &std::collections::HashMap<(u16, u16), PlateInteraction>,
    ) -> GeologicDomain {
        let plate_id = plate_map.get(x, y).copied().unwrap_or(0);
        let distance_to_boundary = distance_field[y][x];
        let plate_stat = plate_stats.get(&plate_id);
        
        match distance_to_boundary {
            0..=2 => {
                // Direct boundary zones
                self.classify_active_boundary_zone(x, y, plate_map, plate_stats, boundary_interactions)
            }
            3..=20 => {
                // Near boundary zones
                self.classify_near_boundary(plate_id, plate_stat)
            }
            _ => {
                // Plate interior - large coherent regions
                self.classify_plate_interior(plate_id, plate_stat)
            }
        }
    }
    
    /// Clean near-boundary classification
    fn classify_near_boundary(
        &self,
        _plate_id: u16,
        plate_stat: Option<&PlateStatistics>,
    ) -> GeologicDomain {
        match plate_stat {
            Some(stat) => {
                if stat.is_oceanic {
                    GeologicDomain::Oceanic(OceanicType::AbyssalPlain)
                } else if stat.age_ma > 1500.0 {
                    GeologicDomain::StableContinental(StableType::Platform)
                } else {
                    GeologicDomain::Extensional(ExtensionalType::ExtendedCrust)
                }
            }
            None => GeologicDomain::Oceanic(OceanicType::AbyssalPlain)
        }
    }
    
    /// Clean plate interior classification
    fn classify_plate_interior(
        &self,
        _plate_id: u16,
        plate_stat: Option<&PlateStatistics>,
    ) -> GeologicDomain {
        match plate_stat {
            Some(stat) => {
                if stat.is_oceanic {
                    GeologicDomain::Oceanic(OceanicType::AbyssalPlain)
                } else if stat.age_ma > 2500.0 && stat.size_ratio > 0.08 {
                    GeologicDomain::StableContinental(StableType::Craton)
                } else if stat.age_ma > 1000.0 {
                    GeologicDomain::StableContinental(StableType::Platform)
                } else {
                    GeologicDomain::StableContinental(StableType::Platform)
                }
            }
            None => GeologicDomain::Oceanic(OceanicType::AbyssalPlain)
        }
    }
    
    /// Classify active boundary zones with realistic geological features
    fn classify_active_boundary_zone(
        &self,
        x: usize,
        y: usize,
        plate_map: &crate::map::terrain::PlateMap,
        plate_stats: &std::collections::HashMap<u16, PlateStatistics>,
        boundary_interactions: &std::collections::HashMap<(u16, u16), PlateInteraction>,
    ) -> GeologicDomain {
        let plate_id = plate_map.get(x, y).copied().unwrap_or(0);
        
        // Find interaction with neighboring plates
        for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
            let nx = (x as i32 + dx) as usize;
            let ny = (y as i32 + dy) as usize;
            
            if nx < self.width && ny < self.height {
                if let Some(&neighbor_plate) = plate_map.get(nx, ny) {
                    if neighbor_plate != plate_id {
                        let key = if plate_id < neighbor_plate {
                            (plate_id, neighbor_plate)
                        } else {
                            (neighbor_plate, plate_id)
                        };
                        
                        if let Some(&interaction) = boundary_interactions.get(&key) {
                            return self.interaction_to_domain(interaction, plate_id, neighbor_plate, plate_stats);
                        }
                    }
                }
            }
        }
        
        // Fallback
        GeologicDomain::Oceanic(OceanicType::FractureZone)
    }
    
    /// Convert plate interaction to geological domain
    fn interaction_to_domain(
        &self,
        interaction: PlateInteraction,
        plate_a: u16,
        plate_b: u16,
        plate_stats: &std::collections::HashMap<u16, PlateStatistics>,
    ) -> GeologicDomain {
        let stat_a = plate_stats.get(&plate_a);
        let stat_b = plate_stats.get(&plate_b);
        
        match interaction {
            PlateInteraction::Convergent => {
                match (stat_a, stat_b) {
                    (Some(a), Some(b)) => {
                        if !a.is_oceanic && !b.is_oceanic {
                            // Continental-continental collision
                            GeologicDomain::Orogenic(OrogenicType::CollisionOrogeny)
                        } else if a.is_oceanic || b.is_oceanic {
                            // Oceanic-continental subduction
                            if a.age_ma > 100.0 || b.age_ma > 100.0 {
                                GeologicDomain::Orogenic(OrogenicType::SubductionOrogeny)
                            } else {
                                GeologicDomain::ArcSystem(ArcType::VolcanicArc)
                            }
                        } else {
                            GeologicDomain::Orogenic(OrogenicType::AccretionaryOrogeny)
                        }
                    }
                    _ => GeologicDomain::Orogenic(OrogenicType::CollisionOrogeny)
                }
            }
            PlateInteraction::Divergent => {
                match (stat_a, stat_b) {
                    (Some(a), Some(b)) if a.is_oceanic && b.is_oceanic => {
                        GeologicDomain::Oceanic(OceanicType::MidOceanRidge)
                    }
                    _ => GeologicDomain::Extensional(ExtensionalType::ContinentalRift)
                }
            }
            PlateInteraction::Transform => {
                match (stat_a, stat_b) {
                    (Some(a), Some(b)) if a.is_oceanic && b.is_oceanic => {
                        GeologicDomain::Oceanic(OceanicType::FractureZone)
                    }
                    _ => GeologicDomain::Extensional(ExtensionalType::ExtendedCrust)
                }
            }
        }
    }
    
    
    
    /// Calculate realistic plate age with more oceanic plates
    fn calculate_realistic_plate_age(&self, plate_id: u16, size_ratio: f64) -> f64 {
        let base_hash = ((plate_id as u64 * 2654435761u64) % 1000) as f64 / 1000.0;
        
        // More oceanic plates - shifted ages younger, favoring oceanic classification
        if size_ratio > 0.15 {
            // Only very large plates are continental (like supercontinents)
            1800.0 + base_hash * 1700.0 // Ancient continental cores (1.8-3.5 Ga)
        } else if size_ratio > 0.08 {
            // Large plates - mix of old continental and younger oceanic
            400.0 + base_hash * 1200.0  // Mixed (400-1600 Ma) 
        } else if size_ratio > 0.03 {
            // Medium plates - mostly oceanic
            50.0 + base_hash * 600.0    // Mostly oceanic (50-650 Ma)
        } else {
            // Small plates - almost entirely oceanic
            5.0 + base_hash * 200.0     // Young oceanic (5-205 Ma)
        }
    }
    
    /// Get plate motion information
    fn get_plate_motion(&self, plate_id: u16) -> (f64, f64) {
        if let Some(seed) = self.plate_seeds.iter().find(|s| s.id == plate_id) {
            (seed.motion_direction, seed.motion_speed)
        } else {
            (0.0, 0.0)
        }
    }
    
    /// Calculate plate center from pixel coordinates
    fn calculate_plate_center(&self, pixels: &[(usize, usize)]) -> (f64, f64) {
        if pixels.is_empty() {
            return (0.0, 0.0);
        }
        
        let sum_x: usize = pixels.iter().map(|(x, _)| x).sum();
        let sum_y: usize = pixels.iter().map(|(_, y)| y).sum();
        let count = pixels.len();
        
        (sum_x as f64 / count as f64, sum_y as f64 / count as f64)
    }
    
    /// Determine interaction between two plates
    fn determine_plate_interaction(
        &self,
        plate_a: u16,
        plate_b: u16,
        plate_stats: &std::collections::HashMap<u16, PlateStatistics>,
    ) -> PlateInteraction {
        match (plate_stats.get(&plate_a), plate_stats.get(&plate_b)) {
            (Some(stat_a), Some(stat_b)) => {
                let motion_diff = (stat_a.motion_direction - stat_b.motion_direction).abs();
                let motion_diff = if motion_diff > 180.0 { 360.0 - motion_diff } else { motion_diff };
                
                if motion_diff < 30.0 {
                    PlateInteraction::Transform
                } else if motion_diff > 150.0 {
                    PlateInteraction::Convergent
                } else {
                    PlateInteraction::Divergent
                }
            }
            _ => PlateInteraction::Transform
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


/// Pre-computed statistics for efficient plate analysis
#[derive(Debug, Clone)]
struct PlateStatistics {
    pub size_ratio: f64,           // Relative size (0.0 to 1.0)
    pub age_ma: f64,               // Age in millions of years
    pub motion_direction: f64,      // Motion direction in degrees
    #[allow(dead_code)]
    pub motion_speed: f64,          // Motion speed (future use)
    #[allow(dead_code)]
    pub center: (f64, f64),        // Plate center coordinates (future use)
    pub is_oceanic: bool,          // Whether this is an oceanic plate
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