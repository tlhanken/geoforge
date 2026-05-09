use std::collections::HashMap;

use crate::coordinates::OrbitalCoordinates;
use crate::objects::orbital_group::OrbitalGroup;
use crate::objects::planetary_body::PlanetaryBody;

#[derive(Debug, Clone, PartialEq)]
pub struct PlanetarySystem {
    pub central_planetary_body: PlanetaryBody,
    pub moons: HashMap<OrbitalCoordinates, PlanetaryBody>,
    pub rings: HashMap<OrbitalCoordinates, OrbitalGroup>,
}
