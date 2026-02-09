//! Electrostatic simulation for tectonic plate generation
//! 
//! This module implements an electrostatic charge simulation on a sphere
//! to generate natural tectonic plate boundaries. Point charges are placed
//! randomly and allowed to reach equilibrium through mutual repulsion.

use crate::map::spherical::SphericalPoint;
use crate::tectonics::plates::PlateSeed;
use crate::tectonics::PlateError;
use rand::prelude::*;
use rand::rngs::StdRng;

/// A point charge on the sphere for electrostatic simulation
#[derive(Debug, Clone)]
pub struct PointCharge {
    pub position: SphericalPoint,
    pub charge: f64,
    pub id: u16,
}

impl PointCharge {
    pub fn new(position: SphericalPoint, charge: f64, id: u16) -> Self {
        Self { position, charge, id }
    }

    /// Calculate electrostatic force from another charge
    pub fn force_from(&self, other: &PointCharge) -> (f64, f64, f64) {
        // Convert to Cartesian for force calculation
        let pos1 = self.position.to_cartesian();
        let pos2 = other.position.to_cartesian();
        
        // Vector from other to self
        let dx = pos1.0 - pos2.0;
        let dy = pos1.1 - pos2.1;
        let dz = pos1.2 - pos2.2;
        
        // Distance between charges
        let distance_sq = dx * dx + dy * dy + dz * dz;
        let distance = distance_sq.sqrt();
        
        // Avoid division by zero for coincident points
        if distance < 1e-10 {
            return (0.0, 0.0, 0.0);
        }
        
        // Coulomb's law: F = k * q1 * q2 / r^2
        // We'll use a simplified constant and scale appropriately
        let force_magnitude = (self.charge * other.charge) / distance_sq;
        
        // Force direction (repulsive, pointing away from other charge)
        let force_x = force_magnitude * (dx / distance);
        let force_y = force_magnitude * (dy / distance);
        let force_z = force_magnitude * (dz / distance);
        
        (force_x, force_y, force_z)
    }
}

/// Electrostatic simulation parameters
pub struct ElectrostaticConfig {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub damping_factor: f64,
    pub time_step: f64,
    pub min_charge: f64,
    pub max_charge: f64,
}

impl Default for ElectrostaticConfig {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            convergence_threshold: 1e-6,
            damping_factor: 0.9,
            time_step: 0.01,
            min_charge: 0.01,  // Microscopic for micro-plates (Juan de Fuca level)
            max_charge: 100.0, // Colossal for superplates (Pacific level)
        }
    }
}

/// Generate initial random charges on the sphere
pub fn generate_random_charges(
    num_charges: usize,
    rng: &mut StdRng,
    config: &ElectrostaticConfig,
) -> Result<Vec<PointCharge>, PlateError> {
    let mut charges = Vec::with_capacity(num_charges);
    
    for i in 0..num_charges {
        // Generate random point on sphere with uniform distribution
        let u = rng.gen::<f64>();
        let v = rng.gen::<f64>();
        
        let theta = 2.0 * std::f64::consts::PI * u;  // longitude
        let phi = (2.0 * v - 1.0).acos();           // latitude from pole
        
        // Convert to lat/lon
        let lat = (std::f64::consts::PI / 2.0 - phi).to_degrees();
        let lon = theta.to_degrees();
        
        // Normalize longitude to [-180, 180]
        let lon = if lon > 180.0 { lon - 360.0 } else { lon };
        
        let position = SphericalPoint::from_lat_lon(lat, lon);
        
        // Power-law charge distribution for realistic size variety
        // Most plates are small, few are large (like Earth)
        let charge = generate_power_law_charge(i, num_charges, rng, config);
        
        charges.push(PointCharge::new(position, charge, (i + 1) as u16));
    }
    
    // Sort charges by magnitude for display
    let mut charge_values: Vec<f64> = charges.iter().map(|c| c.charge).collect();
    charge_values.sort_by(|a, b| b.partial_cmp(a).unwrap());
    
    println!("Generated {} charges with power-law distribution:", charges.len());
    println!("  Largest: {:.2}, Smallest: {:.2}, Ratio: {:.1}x", 
             charge_values[0], charge_values[charge_values.len()-1], 
             charge_values[0] / charge_values[charge_values.len()-1]);
    Ok(charges)
}

/// Generate power-law distributed charges for realistic plate size variety
/// Creates a few large plates and many small ones, like Earth
fn generate_power_law_charge(
    index: usize,
    total_plates: usize,
    rng: &mut StdRng,
    config: &ElectrostaticConfig,
) -> f64 {
    // Extreme Earth-like distribution:
    // - 1 superplate (Pacific-like: ~50% of surface)
    // - 1-2 major plates (Eurasian, African size: ~20% each)  
    // - 2-3 medium plates (North American, Antarctic size: ~8-12%)
    // - Rest are small/micro plates (Caribbean, Juan de Fuca size: <2%)
    
    let superplate_count = 1; // Only 1 massive plate
    let major_plate_count = (total_plates as f64 * 0.15).max(2.0) as usize; // ~15%
    let medium_plate_count = (total_plates as f64 * 0.25).max(3.0) as usize; // ~25%
    // Remaining ~60% are small/micro plates
    
    // Earth-level extreme charge distribution
    let charge = if index < superplate_count {
        // Superplate: 85-100% of max charge (COLOSSAL like Pacific - 50% of Earth's surface)
        rng.gen_range(config.max_charge * 0.85..=config.max_charge)
    } else if index < major_plate_count {
        // Major plates: 35-65% of max charge (Eurasian, African scale)
        rng.gen_range(config.max_charge * 0.35..=config.max_charge * 0.65)
    } else if index < medium_plate_count {
        // Medium plates: 15-35% of max charge (North American, Antarctic scale)
        rng.gen_range(config.max_charge * 0.15..=config.max_charge * 0.35)
    } else {
        // Small/micro plates: Extreme exponential decay for Earth-like variety
        let remaining_plates = total_plates - medium_plate_count;
        let position_in_tail = (index - medium_plate_count) as f64 / remaining_plates as f64;
        
        // More extreme exponential decay targeting 429x ratio
        let base_charge = config.max_charge * 0.15 * (-5.0 * position_in_tail).exp();
        
        // For the tiniest plates, go down to minimum charge
        if position_in_tail > 0.8 {
            rng.gen_range(config.min_charge..=config.min_charge * 2.0)
        } else {
            rng.gen_range(config.min_charge..=base_charge.max(config.min_charge * 5.0))
        }
    };
    
    // Less variation to maintain extreme differences
    let variation = rng.gen_range(0.9..1.1);
    (charge * variation).clamp(config.min_charge, config.max_charge)
}

/// Simulate electrostatic equilibrium using iterative force calculation
pub fn simulate_equilibrium(
    charges: &mut [PointCharge],
    config: &ElectrostaticConfig,
) -> Result<(), PlateError> {
    println!("Starting electrostatic simulation...");
    
    let mut velocities = vec![(0.0, 0.0, 0.0); charges.len()];
    
    for iteration in 0..config.max_iterations {
        let mut max_force: f64 = 0.0;
        let mut total_forces = vec![(0.0, 0.0, 0.0); charges.len()];
        
        // Calculate forces between all pairs of charges
        for i in 0..charges.len() {
            for j in 0..charges.len() {
                if i != j {
                    let force = charges[i].force_from(&charges[j]);
                    total_forces[i].0 += force.0;
                    total_forces[i].1 += force.1;
                    total_forces[i].2 += force.2;
                }
            }
            
            // Track maximum force for convergence
            let force_magnitude = (total_forces[i].0.powi(2) + 
                                  total_forces[i].1.powi(2) + 
                                  total_forces[i].2.powi(2)).sqrt();
            max_force = max_force.max(force_magnitude);
        }
        
        // Check for convergence
        if max_force < config.convergence_threshold {
            println!("Converged after {} iterations (max force: {:.2e})", iteration + 1, max_force);
            break;
        }
        
        // Update positions using velocity Verlet integration
        for i in 0..charges.len() {
            // Update velocity with damping
            velocities[i].0 = velocities[i].0 * config.damping_factor + total_forces[i].0 * config.time_step;
            velocities[i].1 = velocities[i].1 * config.damping_factor + total_forces[i].1 * config.time_step;
            velocities[i].2 = velocities[i].2 * config.damping_factor + total_forces[i].2 * config.time_step;
            
            // Update position
            let current_pos = charges[i].position.to_cartesian();
            let new_x = current_pos.0 + velocities[i].0 * config.time_step;
            let new_y = current_pos.1 + velocities[i].1 * config.time_step;
            let new_z = current_pos.2 + velocities[i].2 * config.time_step;
            
            // Project back to unit sphere
            let norm = (new_x * new_x + new_y * new_y + new_z * new_z).sqrt();
            let normalized_x = new_x / norm;
            let normalized_y = new_y / norm;
            let normalized_z = new_z / norm;
            
            // Convert back to spherical coordinates
            let lat = normalized_z.asin().to_degrees();
            let lon = normalized_y.atan2(normalized_x).to_degrees();
            
            charges[i].position = SphericalPoint::from_lat_lon(lat, lon);
        }
        
        // Progress reporting
        if iteration % 100 == 0 || iteration < 10 {
            println!("Iteration {}: max force = {:.2e}", iteration + 1, max_force);
        }
    }
    
    println!("Electrostatic simulation complete");
    Ok(())
}

/// Convert equilibrium charges to plate seeds
pub fn charges_to_seeds(
    charges: &[PointCharge],
    width: usize,
    height: usize,
    rng: &mut StdRng,
) -> Vec<PlateSeed> {
    charges
        .iter()
        .map(|charge| {
            let (lat, lon) = charge.position.to_lat_lon();
            
            // Convert to pixel coordinates
            let x = ((lon + 180.0) / 360.0 * width as f64) as usize;
            let y = ((90.0 - lat) / 180.0 * height as f64) as usize;
            
            // Clamp to valid range
            let x = x.min(width - 1);
            let y = y.min(height - 1);
            
            PlateSeed::new(
                charge.id,
                x,
                y,
                lat,
                lon,
                rng.gen_range(0.0..360.0),    // Random motion direction
                charge.charge,                 // Use charge as velocity magnitude
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_charge_generation() {
        let mut rng = StdRng::seed_from_u64(42);
        let config = ElectrostaticConfig::default();
        let charges = generate_random_charges(5, &mut rng, &config).unwrap();
        
        assert_eq!(charges.len(), 5);
        
        // Check that charges have valid coordinates
        for charge in &charges {
            let (lat, lon) = charge.position.to_lat_lon();
            assert!(lat >= -90.0 && lat <= 90.0);
            assert!(lon >= -180.0 && lon <= 180.0);
            assert!(charge.charge >= config.min_charge && charge.charge <= config.max_charge);
        }
    }
    
    #[test]
    fn test_force_calculation() {
        let pos1 = SphericalPoint::from_lat_lon(0.0, 0.0);
        let pos2 = SphericalPoint::from_lat_lon(0.0, 90.0);
        
        let charge1 = PointCharge::new(pos1, 1.0, 1);
        let charge2 = PointCharge::new(pos2, 1.0, 2);
        
        let force = charge1.force_from(&charge2);
        
        // Force should be non-zero and repulsive
        let force_magnitude = (force.0.powi(2) + force.1.powi(2) + force.2.powi(2)).sqrt();
        assert!(force_magnitude > 0.0);
    }
    
    #[test]
    fn test_coincident_charges() {
        let pos = SphericalPoint::from_lat_lon(45.0, -120.0);
        let charge1 = PointCharge::new(pos, 1.0, 1);
        let charge2 = PointCharge::new(pos, 1.0, 2);
        
        let force = charge1.force_from(&charge2);
        
        // Force should be zero for coincident charges
        let force_magnitude = (force.0.powi(2) + force.1.powi(2) + force.2.powi(2)).sqrt();
        assert!(force_magnitude < 1e-9);
    }
    
    #[test]
    fn test_power_law_distribution() {
        let mut rng = StdRng::seed_from_u64(42);
        let config = ElectrostaticConfig::default();
        let charges = generate_random_charges(20, &mut rng, &config).unwrap();
        
        let mut charge_values: Vec<f64> = charges.iter().map(|c| c.charge).collect();
        charge_values.sort_by(|a, b| b.partial_cmp(a).unwrap());
        
        // Should have significant size variety
        let largest = charge_values[0];
        let smallest = charge_values[charge_values.len() - 1];
        let ratio = largest / smallest;
        
        // Should achieve Earth-like ratios (hundreds to thousands)
        assert!(ratio > 100.0);
        println!("Charge ratio: {:.1}x", ratio);
    }
    
    #[test]
    fn test_seed_reproducibility() {
        let mut rng1 = StdRng::seed_from_u64(12345);
        let mut rng2 = StdRng::seed_from_u64(12345);
        let config = ElectrostaticConfig::default();
        
        let charges1 = generate_random_charges(10, &mut rng1, &config).unwrap();
        let charges2 = generate_random_charges(10, &mut rng2, &config).unwrap();
        
        // Should produce identical results with same seed
        assert_eq!(charges1.len(), charges2.len());
        for (c1, c2) in charges1.iter().zip(charges2.iter()) {
            assert_eq!(c1.id, c2.id);
            assert!((c1.charge - c2.charge).abs() < 1e-10);
            let (lat1, lon1) = c1.position.to_lat_lon();
            let (lat2, lon2) = c2.position.to_lat_lon();
            assert!((lat1 - lat2).abs() < 1e-10);
            assert!((lon1 - lon2).abs() < 1e-10);
        }
    }
    
    #[test]
    fn test_simulation_convergence() {
        let mut rng = StdRng::seed_from_u64(42);
        let config = ElectrostaticConfig {
            max_iterations: 50,
            convergence_threshold: 1e-3,
            ..Default::default()
        };
        
        let mut charges = generate_random_charges(5, &mut rng, &config).unwrap();
        let result = simulate_equilibrium(&mut charges, &config);
        
        // Should converge or reach max iterations without error
        assert!(result.is_ok());
        
        // Charges should still be on unit sphere
        for charge in charges {
            let (lat, lon) = charge.position.to_lat_lon();
            assert!(lat >= -90.0 && lat <= 90.0);
            assert!(lon >= -180.0 && lon <= 180.0);
        }
    }
}