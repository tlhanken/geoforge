use uom::si::f64::*;
use uom::typenum::{N1, N2, P3, Z0};

/// Inverse area: m⁻² — L⁻² (cosmological constant units)
type InverseArea = uom::si::Quantity<uom::si::ISQ<N2, Z0, Z0, Z0, Z0, Z0, Z0>, uom::si::SI<f64>, f64>;

#[derive(Debug, Clone, PartialEq)]
pub struct Universe {
    pub age: Time,
    pub constants: UniversalConstants,
}

/// Gravitational constant: m³·kg⁻¹·s⁻² — L³ M⁻¹ T⁻²
type GravitationalConstantUnit = uom::si::Quantity<
    uom::si::ISQ<P3, N1, N2, Z0, Z0, Z0, Z0>,
    uom::si::SI<f64>,
    f64,
>;

/// Fundamental constants — fixed by nature, not tunable
#[derive(Debug, Clone, PartialEq)]
pub struct UniversalConstants {
    /// Expansion rate of space itself. Higher = faster universe expansion,
    /// less time for structure to form before everything flies apart.
    pub hubble: Frequency,

    /// Maximum speed of causality. Sets the size of the observable universe
    /// and caps how fast information (and gravity) can propagate.
    pub speed_of_light: Velocity,

    /// Strength of gravity between masses. The primary driver of collapse —
    /// stars, structure, everything. Nudge this and stellar lifetimes change dramatically.
    pub gravitational: GravitationalConstantUnit,

    /// Smallest meaningful length scale. Below this, spacetime itself breaks down.
    /// Marks where quantum gravity effects dominate.
    pub planck_length: Length,

    /// Smallest meaningful time interval — one Planck length / c.
    /// The "tick rate" of the universe at quantum scales.
    pub planck_time: Time,

    /// Temperature of the CMB today. Encodes the universe's thermal history
    /// since recombination at 380,000 years old. Redshifts lower as universe expands.
    pub cmb_temperature: ThermodynamicTemperature,

    /// Current age of the universe. Constrains how long structure has had to form,
    /// and sets the upper bound on the particle horizon.
    pub age: Time,

    /// Dark energy density / spacetime curvature term. Drives accelerating expansion.
    /// If too large, overcomes gravity before anything can bind.
    pub cosmological_constant: InverseArea,
}