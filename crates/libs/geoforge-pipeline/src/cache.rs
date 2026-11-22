//! Derived fields worth keeping between commands — hazard rates and the like, which take
//! seconds on a large grid — in `cache/derived/`, reused while the stored products and the
//! body they come from are unchanged. Best effort: an unreadable or stale entry is simply
//! worked out again.

use geoforge_core::product::ProductId;
use geoforge_store::{hash_parts, World};

/// The field `name` (bump `version` when its formula changes), from the cache when the
/// world is unchanged since it was kept, otherwise from `compute`, then kept.
pub fn field(world: &World, name: &str, version: u32, compute: impl FnOnce() -> Vec<f64>) -> Vec<f64> {
    let key = key(world, name, version);
    let path = world.root.join("cache").join("derived").join(format!("{name}.bin"));
    if let Ok(bytes) = std::fs::read(&path) {
        if bytes.len() >= 32 && bytes[..32] == key && (bytes.len() - 32) % 8 == 0 {
            return bytes[32..].chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().expect("8 bytes"))).collect();
        }
    }
    let values = compute();
    let mut bytes = key.to_vec();
    bytes.extend(values.iter().flat_map(|v| v.to_le_bytes()));
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir).and_then(|()| std::fs::write(&path, bytes));
    }
    values
}

/// What the field depends on: its name and version, the body, and every stored product's
/// provenance.
fn key(world: &World, name: &str, version: u32) -> [u8; 32] {
    let body = std::fs::read(world.body_dir().join("body.json")).unwrap_or_default();
    let provenance: Vec<u8> = ProductId::ALL.iter().flat_map(|&p| world.stored_provenance(p).map_or([0; 32], |h| h.0)).collect();
    hash_parts(&[name.as_bytes(), &version.to_le_bytes(), &body, &provenance]).0
}
