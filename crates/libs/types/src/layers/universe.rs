use super::Seed;
use crate::objects::Universe;

#[derive(Debug, Clone, PartialEq)]
pub struct UniverseLayer {
    pub seed: Seed,
    pub universe: Universe,
}
