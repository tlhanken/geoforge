//! Products, sources, provenance and the generator contract.
//!
//! The unit that is declared, stored, injected, exported and hashed is a *product*: one
//! field, table, partition, network or timeline. Each product records where it comes from.

use crate::body::BodyContext;
use crate::grid::Grid;
use crate::seed::Seed;
use crate::shapes::{Field, Network, Partition, Timeline};
use crate::surface::{BiomeClass, CrustType, Lake, Plate, Province, RiverNode, RiverReach, Rock, Watershed};
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

/// Layers, in canonical order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LayerId {
    Tectonics,
    Geology,
    Elevation,
    Circulation,
    Climate,
    Erosion,
    Hydrology,
    Biomes,
}

impl LayerId {
    pub const ALL: [LayerId; 8] = [
        LayerId::Tectonics,
        LayerId::Geology,
        LayerId::Elevation,
        LayerId::Circulation,
        LayerId::Climate,
        LayerId::Erosion,
        LayerId::Hydrology,
        LayerId::Biomes,
    ];
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            LayerId::Tectonics => "tectonics",
            LayerId::Geology => "geology",
            LayerId::Elevation => "elevation",
            LayerId::Circulation => "circulation",
            LayerId::Climate => "climate",
            LayerId::Erosion => "erosion",
            LayerId::Hydrology => "hydrology",
            LayerId::Biomes => "biomes",
        }
    }
    #[must_use]
    pub fn from_name(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|l| l.name() == s)
    }
    #[must_use]
    pub fn products(self) -> Vec<ProductId> {
        ProductId::ALL.into_iter().filter(|p| p.layer() == self).collect()
    }
}

/// The data shape of a product.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Shape {
    Table,
    Field,
    Partition,
    Network,
    Timeline,
}

/// Element type of a field (or a partition's owner field).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ValueType {
    F32,
    U8,
    U32,
    /// Two f32 per cell: an (east, north) vector.
    V2F32,
}

/// Every primary product. Derived values (boundaries, slope, sea level, flow) are not here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ProductId {
    Plates,
    CrustType,
    CrustAge,
    TectonicEvents,
    Provinces,
    GeologyEvents,
    Lithology,
    RockAge,
    RoughElevation,
    WindJan,
    WindJul,
    OceanCurrents,
    SstAnomaly,
    Temperature,
    TemperatureRange,
    Precipitation,
    PrecipitationWarmShare,
    Elevation,
    Rivers,
    Lakes,
    Watersheds,
    Biome,
}

impl ProductId {
    pub const ALL: [ProductId; 22] = [
        ProductId::Plates,
        ProductId::CrustType,
        ProductId::CrustAge,
        ProductId::TectonicEvents,
        ProductId::Provinces,
        ProductId::GeologyEvents,
        ProductId::Lithology,
        ProductId::RockAge,
        ProductId::RoughElevation,
        ProductId::WindJan,
        ProductId::WindJul,
        ProductId::OceanCurrents,
        ProductId::SstAnomaly,
        ProductId::Temperature,
        ProductId::TemperatureRange,
        ProductId::Precipitation,
        ProductId::PrecipitationWarmShare,
        ProductId::Elevation,
        ProductId::Rivers,
        ProductId::Lakes,
        ProductId::Watersheds,
        ProductId::Biome,
    ];

    /// Stable name: used in file names and seed derivation. Frozen.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            ProductId::Plates => "plates",
            ProductId::CrustType => "crust_type",
            ProductId::CrustAge => "crust_age",
            ProductId::TectonicEvents => "tectonic_events",
            ProductId::Provinces => "provinces",
            ProductId::GeologyEvents => "geology_events",
            ProductId::Lithology => "lithology",
            ProductId::RockAge => "rock_age",
            ProductId::RoughElevation => "rough_elevation",
            ProductId::Elevation => "elevation",
            ProductId::WindJan => "wind_jan",
            ProductId::WindJul => "wind_jul",
            ProductId::OceanCurrents => "ocean_currents",
            ProductId::SstAnomaly => "sst_anomaly",
            ProductId::Temperature => "temperature",
            ProductId::TemperatureRange => "temperature_range",
            ProductId::Precipitation => "precipitation",
            ProductId::PrecipitationWarmShare => "precipitation_warm_share",
            ProductId::Rivers => "rivers",
            ProductId::Lakes => "lakes",
            ProductId::Watersheds => "watersheds",
            ProductId::Biome => "biome",
        }
    }

    #[must_use]
    pub fn from_name(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|p| p.name() == s)
    }

    #[must_use]
    pub fn layer(self) -> LayerId {
        match self {
            ProductId::Plates | ProductId::CrustType | ProductId::CrustAge | ProductId::TectonicEvents => LayerId::Tectonics,
            ProductId::Provinces | ProductId::GeologyEvents | ProductId::Lithology | ProductId::RockAge => LayerId::Geology,
            ProductId::RoughElevation => LayerId::Elevation,
            ProductId::Elevation => LayerId::Erosion,
            ProductId::WindJan | ProductId::WindJul | ProductId::OceanCurrents | ProductId::SstAnomaly => LayerId::Circulation,
            ProductId::Temperature | ProductId::TemperatureRange | ProductId::Precipitation | ProductId::PrecipitationWarmShare => LayerId::Climate,
            ProductId::Rivers | ProductId::Lakes | ProductId::Watersheds => LayerId::Hydrology,
            ProductId::Biome => LayerId::Biomes,
        }
    }

    #[must_use]
    pub fn shape(self) -> Shape {
        match self {
            ProductId::Plates | ProductId::Provinces | ProductId::Lithology | ProductId::Lakes | ProductId::Watersheds | ProductId::Biome => {
                Shape::Partition
            }
            ProductId::CrustType
            | ProductId::CrustAge
            | ProductId::RockAge
            | ProductId::RoughElevation
            | ProductId::Elevation
            | ProductId::WindJan
            | ProductId::WindJul
            | ProductId::OceanCurrents
            | ProductId::SstAnomaly
            | ProductId::Temperature
            | ProductId::TemperatureRange
            | ProductId::Precipitation
            | ProductId::PrecipitationWarmShare => Shape::Field,
            ProductId::TectonicEvents | ProductId::GeologyEvents => Shape::Timeline,
            ProductId::Rivers => Shape::Network,
        }
    }

    /// A partition whose table is a catalog of kinds (rocks, biomes) rather than of regions:
    /// every row is kept even when no cell uses it, and each row carries its own colour.
    #[must_use]
    pub fn is_catalog(self) -> bool {
        matches!(self, ProductId::Lithology | ProductId::Biome)
    }

    /// Element type for fields and partition owners.
    #[must_use]
    pub fn value_type(self) -> Option<ValueType> {
        match self.shape() {
            Shape::Partition => Some(ValueType::U32),
            Shape::Field => Some(match self {
                ProductId::CrustType => ValueType::U8,
                ProductId::WindJan | ProductId::WindJul | ProductId::OceanCurrents => ValueType::V2F32,
                _ => ValueType::F32,
            }),
            _ => None,
        }
    }

    /// Unit of a continuous field.
    #[must_use]
    pub fn unit(self) -> &'static str {
        match self {
            ProductId::CrustAge | ProductId::RockAge => "Ma",
            ProductId::RoughElevation | ProductId::Elevation => "m",
            ProductId::WindJan | ProductId::WindJul | ProductId::OceanCurrents => "m/s",
            ProductId::Temperature | ProductId::TemperatureRange | ProductId::SstAnomaly => "degC",
            ProductId::Precipitation => "mm/yr",
            ProductId::PrecipitationWarmShare => "fraction",
            _ => "",
        }
    }
}

impl fmt::Display for ProductId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// 32-byte content hash, hex-encoded in files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct Hash(pub [u8; 32]);

impl fmt::Display for Hash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for b in self.0 {
            write!(f, "{b:02x}")?;
        }
        Ok(())
    }
}

impl Serialize for Hash {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Hash {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let s = String::deserialize(d)?;
        Hash::from_hex(&s).ok_or_else(|| serde::de::Error::custom("hash must be 64 hex digits"))
    }
}

impl Hash {
    #[must_use]
    pub fn from_hex(s: &str) -> Option<Self> {
        if s.len() != 64 {
            return None;
        }
        let mut out = [0u8; 32];
        for (i, b) in out.iter_mut().enumerate() {
            *b = u8::from_str_radix(&s[2 * i..2 * i + 2], 16).ok()?;
        }
        Some(Self(out))
    }
}

/// Where a product comes from. Recorded per product in the world manifest.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Source {
    /// Seed only, or seed plus settings.
    Generate {
        #[serde(default, skip_serializing_if = "serde_json::Value::is_null")]
        params: serde_json::Value,
    },
    /// Generate, then apply your stored edits. `base` is the hash of the generated product
    /// the patch was made against.
    Patched {
        #[serde(default, skip_serializing_if = "serde_json::Value::is_null")]
        params: serde_json::Value,
        base: Hash,
        patch: String,
    },
    /// Your data replaces the product; the generator is skipped.
    Authored { file: String },
}

impl Source {
    #[must_use]
    pub fn generate() -> Self {
        Source::Generate { params: serde_json::Value::Null }
    }
    #[must_use]
    pub fn is_user_made(&self) -> bool {
        !matches!(self, Source::Generate { .. })
    }
    #[must_use]
    pub fn params(&self) -> &serde_json::Value {
        static NULL: serde_json::Value = serde_json::Value::Null;
        match self {
            Source::Generate { params } | Source::Patched { params, .. } => params,
            Source::Authored { .. } => &NULL,
        }
    }
    #[must_use]
    pub fn label(&self) -> &'static str {
        match self {
            Source::Generate { .. } => "generate",
            Source::Patched { .. } => "patched",
            Source::Authored { .. } => "authored",
        }
    }
}

/// Wrapper around every stored product document.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Envelope<T> {
    pub schema_version: u32,
    pub kind: ProductId,
    /// Hash of everything this product was built from.
    pub provenance: Hash,
    pub body: T,
}

/// Why a layer does not apply to a body.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Precondition(pub String);

impl fmt::Display for Precondition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// State of a product in a world.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "state", rename_all = "snake_case")]
pub enum Status {
    Fresh,
    Stale,
    /// A patch was re-applied to a regenerated base.
    Rebased,
    /// Some of a patch's edits refer to cells or keys that no longer exist.
    Orphaned { count: usize },
    NotApplicable { reason: String },
    Missing,
}

/// A product's data, typed.
#[derive(Debug, Clone, PartialEq)]
pub enum ProductData {
    Plates(Partition<Plate>),
    CrustType(Field<CrustType>),
    CrustAge(Field<f32>),
    TectonicEvents(Timeline),
    Provinces(Partition<Province>),
    GeologyEvents(Timeline),
    Lithology(Partition<Rock>),
    RockAge(Field<f32>),
    RoughElevation(Field<f32>),
    WindJan(Field<[f32; 2]>),
    WindJul(Field<[f32; 2]>),
    OceanCurrents(Field<[f32; 2]>),
    SstAnomaly(Field<f32>),
    Temperature(Field<f32>),
    TemperatureRange(Field<f32>),
    Precipitation(Field<f32>),
    PrecipitationWarmShare(Field<f32>),
    Elevation(Field<f32>),
    Rivers(Network<RiverNode, RiverReach>),
    Lakes(Partition<Lake>),
    Watersheds(Partition<Watershed>),
    Biome(Partition<BiomeClass>),
}

impl ProductData {
    #[must_use]
    pub fn id(&self) -> ProductId {
        match self {
            ProductData::Plates(_) => ProductId::Plates,
            ProductData::CrustType(_) => ProductId::CrustType,
            ProductData::CrustAge(_) => ProductId::CrustAge,
            ProductData::TectonicEvents(_) => ProductId::TectonicEvents,
            ProductData::Provinces(_) => ProductId::Provinces,
            ProductData::GeologyEvents(_) => ProductId::GeologyEvents,
            ProductData::Lithology(_) => ProductId::Lithology,
            ProductData::RockAge(_) => ProductId::RockAge,
            ProductData::RoughElevation(_) => ProductId::RoughElevation,
            ProductData::Elevation(_) => ProductId::Elevation,
            ProductData::WindJan(_) => ProductId::WindJan,
            ProductData::WindJul(_) => ProductId::WindJul,
            ProductData::OceanCurrents(_) => ProductId::OceanCurrents,
            ProductData::SstAnomaly(_) => ProductId::SstAnomaly,
            ProductData::Temperature(_) => ProductId::Temperature,
            ProductData::TemperatureRange(_) => ProductId::TemperatureRange,
            ProductData::Precipitation(_) => ProductId::Precipitation,
            ProductData::PrecipitationWarmShare(_) => ProductId::PrecipitationWarmShare,
            ProductData::Rivers(_) => ProductId::Rivers,
            ProductData::Lakes(_) => ProductId::Lakes,
            ProductData::Watersheds(_) => ProductId::Watersheds,
            ProductData::Biome(_) => ProductId::Biome,
        }
    }

    /// The values of a continuous (f32) field product.
    #[must_use]
    pub fn scalar(&self) -> Option<&Field<f32>> {
        match self {
            ProductData::CrustAge(f)
            | ProductData::RockAge(f)
            | ProductData::RoughElevation(f)
            | ProductData::Elevation(f)
            | ProductData::SstAnomaly(f)
            | ProductData::Temperature(f)
            | ProductData::TemperatureRange(f)
            | ProductData::Precipitation(f)
            | ProductData::PrecipitationWarmShare(f) => Some(f),
            _ => None,
        }
    }

    #[must_use]
    pub fn scalar_mut(&mut self) -> Option<&mut Field<f32>> {
        match self {
            ProductData::CrustAge(f)
            | ProductData::RockAge(f)
            | ProductData::RoughElevation(f)
            | ProductData::Elevation(f)
            | ProductData::SstAnomaly(f)
            | ProductData::Temperature(f)
            | ProductData::TemperatureRange(f)
            | ProductData::Precipitation(f)
            | ProductData::PrecipitationWarmShare(f) => Some(f),
            _ => None,
        }
    }

    /// Build a continuous field product from its id.
    #[must_use]
    pub fn from_scalar(id: ProductId, f: Field<f32>) -> Option<Self> {
        Some(match id {
            ProductId::CrustAge => ProductData::CrustAge(f),
            ProductId::RockAge => ProductData::RockAge(f),
            ProductId::RoughElevation => ProductData::RoughElevation(f),
            ProductId::Elevation => ProductData::Elevation(f),
            ProductId::SstAnomaly => ProductData::SstAnomaly(f),
            ProductId::Temperature => ProductData::Temperature(f),
            ProductId::TemperatureRange => ProductData::TemperatureRange(f),
            ProductId::Precipitation => ProductData::Precipitation(f),
            ProductId::PrecipitationWarmShare => ProductData::PrecipitationWarmShare(f),
            _ => return None,
        })
    }

    /// The values of a vector field product.
    #[must_use]
    pub fn vector(&self) -> Option<&Field<[f32; 2]>> {
        match self {
            ProductData::WindJan(f) | ProductData::WindJul(f) | ProductData::OceanCurrents(f) => Some(f),
            _ => None,
        }
    }

    #[must_use]
    pub fn vector_mut(&mut self) -> Option<&mut Field<[f32; 2]>> {
        match self {
            ProductData::WindJan(f) | ProductData::WindJul(f) | ProductData::OceanCurrents(f) => Some(f),
            _ => None,
        }
    }

    /// Build a vector field product from its id.
    #[must_use]
    pub fn from_vector(id: ProductId, f: Field<[f32; 2]>) -> Option<Self> {
        Some(match id {
            ProductId::WindJan => ProductData::WindJan(f),
            ProductId::WindJul => ProductData::WindJul(f),
            ProductId::OceanCurrents => ProductData::OceanCurrents(f),
            _ => return None,
        })
    }

    /// Grid of a field or partition product.
    #[must_use]
    pub fn grid(&self) -> Option<crate::grid::GridSpec> {
        if let Some(f) = self.scalar() {
            return Some(f.grid);
        }
        if let Some(f) = self.vector() {
            return Some(f.grid);
        }
        match self {
            ProductData::CrustType(f) => Some(f.grid),
            ProductData::Plates(p) => Some(p.owner.grid),
            ProductData::Lithology(p) => Some(p.owner.grid),
            ProductData::Biome(p) => Some(p.owner.grid),
            ProductData::Provinces(p) => Some(p.owner.grid),
            ProductData::Lakes(p) => Some(p.owner.grid),
            ProductData::Watersheds(p) => Some(p.owner.grid),
            _ => None,
        }
    }
}

/// A bag of products by id.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Products {
    map: BTreeMap<ProductId, ProductData>,
}

macro_rules! getters {
    ($($fn:ident => $var:ident : $ty:ty),+ $(,)?) => {
        impl Products {
            $(
                #[must_use]
                pub fn $fn(&self) -> Option<&$ty> {
                    match self.map.get(&ProductId::$var) { Some(ProductData::$var(x)) => Some(x), _ => None }
                }
            )+
        }
    };
}

getters! {
    plates => Plates: Partition<Plate>,
    crust_type => CrustType: Field<CrustType>,
    crust_age => CrustAge: Field<f32>,
    tectonic_events => TectonicEvents: Timeline,
    provinces => Provinces: Partition<Province>,
    geology_events => GeologyEvents: Timeline,
    lithology => Lithology: Partition<Rock>,
    rock_age => RockAge: Field<f32>,
    rough_elevation => RoughElevation: Field<f32>,
    elevation => Elevation: Field<f32>,
    wind_jan => WindJan: Field<[f32; 2]>,
    wind_jul => WindJul: Field<[f32; 2]>,
    ocean_currents => OceanCurrents: Field<[f32; 2]>,
    sst_anomaly => SstAnomaly: Field<f32>,
    temperature => Temperature: Field<f32>,
    temperature_range => TemperatureRange: Field<f32>,
    precipitation => Precipitation: Field<f32>,
    precipitation_warm_share => PrecipitationWarmShare: Field<f32>,
    rivers => Rivers: Network<RiverNode, RiverReach>,
    lakes => Lakes: Partition<Lake>,
    watersheds => Watersheds: Partition<Watershed>,
    biome => Biome: Partition<BiomeClass>,
}

impl Products {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
    pub fn insert(&mut self, data: ProductData) {
        self.map.insert(data.id(), data);
    }
    #[must_use]
    pub fn get(&self, id: ProductId) -> Option<&ProductData> {
        self.map.get(&id)
    }
    pub fn remove(&mut self, id: ProductId) -> Option<ProductData> {
        self.map.remove(&id)
    }
    #[must_use]
    pub fn contains(&self, id: ProductId) -> bool {
        self.map.contains_key(&id)
    }
    pub fn extend(&mut self, other: Products) {
        self.map.extend(other.map);
    }
    pub fn iter(&self) -> impl Iterator<Item = (&ProductId, &ProductData)> {
        self.map.iter()
    }
    pub fn into_values(self) -> impl Iterator<Item = ProductData> {
        self.map.into_values()
    }
}

/// What a generator gets when it runs.
pub struct RunContext<'a> {
    pub root_seed: Seed,
    pub body: &'a BodyContext,
    pub grid: &'a dyn Grid,
}

impl RunContext<'_> {
    /// The seed for one of this generator's products.
    #[must_use]
    pub fn seed(&self, product: ProductId) -> Seed {
        self.root_seed.for_product(product, self.body.key)
    }
    /// Kilometres per radian on this body.
    #[must_use]
    pub fn radius_km(&self) -> f64 {
        self.body.exterior.radius_km
    }
}

/// One implementation of a layer: declares inputs and applicability, and fills the layer's
/// products. A layer can have several implementations (plate tectonics, a rigid lid, …); the
/// pipeline runs the first whose [`applies`](Self::applies) holds for the body.
pub trait Generator {
    type Params: Default + Serialize + DeserializeOwned;
    /// Unique name; keys the implementation's parameters (`[params.<name>]`) and provenance.
    const NAME: &'static str;
    /// Bump when the algorithm changes output for the same inputs: stored products of this
    /// implementation, and everything downstream, then go stale.
    const VERSION: u32;
    const LAYER: LayerId;

    /// Required inputs.
    fn needs() -> &'static [ProductId];
    /// Optional inputs, used when present.
    fn wants() -> &'static [ProductId] {
        &[]
    }
    /// # Errors
    /// Explains why the layer does not apply to this body.
    fn applies(body: &BodyContext) -> Result<(), Precondition>;
    /// Produce every product of [`Self::LAYER`] that is not already in `inputs`
    /// (authored or patched products are passed in and must be respected).
    ///
    /// # Errors
    /// Missing required inputs or invalid parameters.
    fn run(ctx: &RunContext<'_>, inputs: &Products, params: &Self::Params) -> Result<Products, crate::Error>;
}
