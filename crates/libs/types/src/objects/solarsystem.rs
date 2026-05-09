use std::collections::HashMap;

use crate::coordinates::orbital::OrbitalCoordinates;
use crate::objects::orbital_group::OrbitalGroup;
use crate::objects::planetary_body::PlanetaryBody;
use crate::objects::star::Star;

#[derive(Debug, Clone, PartialEq)]
pub struct SolarSystem {
    pub stars: Vec<Star>,
    pub planetary_bodies: HashMap<OrbitalCoordinates, PlanetaryBody>,
    pub orbital_groups: HashMap<OrbitalCoordinates, OrbitalGroup>,
}
