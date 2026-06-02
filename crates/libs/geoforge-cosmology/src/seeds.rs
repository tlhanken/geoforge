//! Seed paths for cosmology hierarchy.

use geoforge_types::cosmology::RegionId;
use geoforge_types::Seed;

/// Root → galaxy `g`.
#[must_use]
pub fn galaxy_seed(root: Seed, galaxy_index: u32) -> Seed {
    root.child(b"galaxy").child_bytes(&galaxy_index.to_le_bytes())
}

/// Galaxy → region.
#[must_use]
pub fn region_seed(galaxy: Seed, region: RegionId) -> Seed {
    let bytes = region_bytes(region);
    galaxy.child(b"region").child_bytes(&bytes)
}

/// Region → system `index`.
#[must_use]
pub fn system_seed(region: Seed, index_in_region: u32) -> Seed {
    region.child(b"system").child_bytes(&index_in_region.to_le_bytes())
}

fn region_bytes(region: RegionId) -> Vec<u8> {
    match region {
        RegionId::Cylindrical { ring, wedge, layer } => {
            let mut v = vec![0u8];
            v.extend_from_slice(&ring.to_le_bytes());
            v.extend_from_slice(&wedge.to_le_bytes());
            v.extend_from_slice(&layer.to_le_bytes());
            v
        }
        RegionId::Spherical {
            shell,
            theta_wedge,
            phi_wedge,
        } => {
            let mut v = vec![1u8];
            v.extend_from_slice(&shell.to_le_bytes());
            v.extend_from_slice(&theta_wedge.to_le_bytes());
            v.extend_from_slice(&phi_wedge.to_le_bytes());
            v
        }
    }
}
