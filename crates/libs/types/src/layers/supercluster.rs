use std::collections::HashMap;
use super::universe::UniverseLayer;
use crate::coordinates::universe::UniverseCoordinates;
use crate::objects::Supercluster;

#[derive(Debug, Clone, PartialEq)]
pub struct SuperclusterLayer {
    pub universe_layer: UniverseLayer,
    pub superclusters: HashMap<UniverseCoordinates, Supercluster>,
}
