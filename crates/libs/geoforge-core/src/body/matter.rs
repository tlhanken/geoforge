//! What shells are made of.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Species {
    Iron,
    Nickel,
    Sulfur,
    Silicate,
    Basalt,
    Granite,
    Anorthosite,
    Water,
    Ice,
    Hydrogen,
    Helium,
    Nitrogen,
    Oxygen,
    CarbonDioxide,
    Argon,
    Methane,
    Ammonia,
}

impl Species {
    /// Chemical formula (rocks by name).
    #[must_use]
    pub fn formula(self) -> &'static str {
        match self {
            Self::Iron => "Fe",
            Self::Nickel => "Ni",
            Self::Sulfur => "S",
            Self::Silicate => "MgSiO₃",
            Self::Basalt => "basalt",
            Self::Granite => "granite",
            Self::Anorthosite => "anorthosite",
            Self::Water => "H₂O",
            Self::Ice => "H₂O ice",
            Self::Hydrogen => "H₂",
            Self::Helium => "He",
            Self::Nitrogen => "N₂",
            Self::Oxygen => "O₂",
            Self::CarbonDioxide => "CO₂",
            Self::Argon => "Ar",
            Self::Methane => "CH₄",
            Self::Ammonia => "NH₃",
        }
    }

    /// Molar mass, g/mol. Rocks use a representative mineral's formula unit; only ratios
    /// within one shell matter (mole fractions of an atmosphere).
    #[must_use]
    pub fn molar_mass_g_mol(self) -> f64 {
        match self {
            Self::Iron => 55.845,
            Self::Nickel => 58.693,
            Self::Sulfur => 32.06,
            Self::Silicate => 100.39,   // MgSiO₃
            Self::Basalt => 250.0,      // plagioclase + pyroxene
            Self::Granite => 200.0,     // feldspar + quartz
            Self::Anorthosite => 278.2, // CaAl₂Si₂O₈
            Self::Water | Self::Ice => 18.015,
            Self::Hydrogen => 2.016,
            Self::Helium => 4.003,
            Self::Nitrogen => 28.014,
            Self::Oxygen => 31.998,
            Self::CarbonDioxide => 44.01,
            Self::Argon => 39.948,
            Self::Methane => 16.043,
            Self::Ammonia => 17.031,
        }
    }
}

/// Species fractions by mass, sorted, summing to 1.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "Vec<(Species, f64)>", into = "Vec<(Species, f64)>")]
pub struct Composition(Vec<(Species, f64)>);

impl Composition {
    /// # Errors
    /// Fractions must be non-negative, species unique, and sum to 1 within 1e-6.
    pub fn new(mut parts: Vec<(Species, f64)>) -> Result<Self, crate::Error> {
        parts.sort_by_key(|a| a.0);
        if parts.windows(2).any(|w| w[0].0 == w[1].0) {
            return Err(crate::Error::Invalid("composition lists a species twice".into()));
        }
        if parts.iter().any(|p| p.1 < 0.0 || !p.1.is_finite()) {
            return Err(crate::Error::Invalid("composition fractions must be >= 0".into()));
        }
        let sum: f64 = parts.iter().map(|p| p.1).sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(crate::Error::Invalid(format!("composition sums to {sum}, not 1")));
        }
        Ok(Self(parts))
    }
    /// Normalize arbitrary positive mass weights.
    #[must_use]
    pub fn from_weights(parts: &[(Species, f64)]) -> Self {
        let sum: f64 = parts.iter().map(|p| p.1).sum();
        let mut v: Vec<_> = parts.iter().map(|&(s, w)| (s, w / sum)).collect();
        v.sort_by_key(|a| a.0);
        Self(v)
    }
    /// From amounts by volume (mole fractions), as atmospheres are usually quoted.
    #[must_use]
    pub fn from_moles(parts: &[(Species, f64)]) -> Self {
        let mass: Vec<_> = parts.iter().map(|&(s, x)| (s, x * s.molar_mass_g_mol())).collect();
        Self::from_weights(&mass)
    }
    #[must_use]
    pub fn single(s: Species) -> Self {
        Self(vec![(s, 1.0)])
    }
    #[must_use]
    pub fn parts(&self) -> &[(Species, f64)] {
        &self.0
    }
    /// Mass fraction of `s`.
    #[must_use]
    pub fn fraction(&self, s: Species) -> f64 {
        self.0.iter().find(|p| p.0 == s).map_or(0.0, |p| p.1)
    }
    /// Mole (volume) fraction of `s`.
    #[must_use]
    pub fn mole_fraction(&self, s: Species) -> f64 {
        let moles = |&(sp, w): &(Species, f64)| w / sp.molar_mass_g_mol();
        let total: f64 = self.0.iter().map(moles).sum();
        self.0.iter().find(|p| p.0 == s).map_or(0.0, |p| moles(p) / total)
    }
}

impl TryFrom<Vec<(Species, f64)>> for Composition {
    type Error = crate::Error;
    fn try_from(v: Vec<(Species, f64)>) -> Result<Self, Self::Error> {
        Self::new(v)
    }
}
impl From<Composition> for Vec<(Species, f64)> {
    fn from(c: Composition) -> Self {
        c.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Phase {
    Solid,
    Liquid,
    Gas,
    Supercritical,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn composition_checks() {
        assert!(Composition::new(vec![(Species::Water, 0.5), (Species::Ice, 0.5)]).is_ok());
        assert!(Composition::new(vec![(Species::Water, 0.5)]).is_err());
        assert!(Composition::new(vec![(Species::Water, 0.5), (Species::Water, 0.5)]).is_err());
    }

    #[test]
    fn moles_round_trip() {
        // Earth's air by volume: 78.1% N₂, 20.9% O₂, 0.93% Ar, 0.04% CO₂.
        let air = Composition::from_moles(&[(Species::Nitrogen, 0.781), (Species::Oxygen, 0.209), (Species::Argon, 0.0093), (Species::CarbonDioxide, 0.0004)]);
        assert!((air.fraction(Species::Nitrogen) - 0.755).abs() < 0.002, "mass N₂ {}", air.fraction(Species::Nitrogen));
        assert!((air.mole_fraction(Species::CarbonDioxide) - 0.0004).abs() < 1e-6);
        assert_eq!(air.mole_fraction(Species::Methane), 0.0);
    }
}
