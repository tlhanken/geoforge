use std::collections::HashMap;
use super::supercluster::SuperclusterLayer;
use crate::coordinates::supercluster::SuperclusterCoordinates;
use crate::objects::Galaxy;

#[derive(Debug, Clone, PartialEq)]
pub struct GalaxyLayer {
    pub supercluster_layer: SuperclusterLayer,
    pub galaxies: HashMap<SuperclusterCoordinates, Galaxy>,
}
