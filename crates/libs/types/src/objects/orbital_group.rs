use uom::si::f64::{Length, Mass, MassDensity, Ratio};

/// Represents a collection of small bodies (asteroids, dust, ice) in an orbital shell.
/// Used for both planetary rings and solar asteroid belts.
#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalGroup {
    /// Width of the belt/ring
    pub width: Length,
    /// Total mass of all bodies in the group
    pub total_mass: Mass,
    /// Primary material composition
    pub density: MassDensity,
    /// How densely packed the group is
    pub number_density: MassDensity,
    /// Min and max size of individual particles/objects
    pub particle_size_range: (Length, Length),
    /// Visual opacity (0.0 = transparent, 1.0 = opaque)
    pub optical_depth: Ratio,
}
