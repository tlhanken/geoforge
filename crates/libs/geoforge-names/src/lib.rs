//! Pronounceable names from a seed, deterministically: the same seed always gives the same
//! name, on any machine. No dependencies, so it can name things in any project.
//!
//! ```
//! assert_eq!(geoforge_names::name(12345), geoforge_names::name(12345));
//! ```
//!
//! For now there is one sound set: 16 onsets and 8 nuclei, two or three open syllables
//! ("Thumi", "Kibedo"). Languages with their own sound inventories, syllable shapes and
//! spelling, per world, culture or region, are planned.
//!
//! The derivation (SplitMix64 mixing, the labels, the draw order) is **frozen**: changing it
//! renames everything already generated. The golden test pins it.

const ONSETS: [&str; 16] = ["k", "t", "v", "s", "r", "m", "n", "l", "d", "th", "z", "b", "h", "sh", "g", "p"];
const NUCLEI: [&str; 8] = ["a", "e", "i", "o", "u", "ae", "ai", "or"];

/// A capitalised name of two or three syllables for `seed`.
#[must_use]
pub fn name(seed: u64) -> String {
    let base = child(seed, b"name");
    let mut n = 0u64;
    let mut u = || {
        n += 1;
        unit(child_u64(base, b"draw", n))
    };
    let syllables = 2 + (u() * 2.0) as usize;
    let mut s = String::new();
    for _ in 0..syllables {
        s.push_str(ONSETS[(u() * ONSETS.len() as f64) as usize % ONSETS.len()]);
        s.push_str(NUCLEI[(u() * NUCLEI.len() as f64) as usize % NUCLEI.len()]);
    }
    let mut c = s.chars();
    c.next().map(|f| f.to_uppercase().collect::<String>() + c.as_str()).unwrap_or_default()
}

/// SplitMix64 finalizer.
fn mix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn child(seed: u64, label: &[u8]) -> u64 {
    mix(label.iter().fold(seed, |h, &b| mix(h ^ u64::from(b))))
}

fn child_u64(seed: u64, label: &[u8], n: u64) -> u64 {
    mix(child(seed, label) ^ mix(n))
}

fn unit(seed: u64) -> f64 {
    (mix(seed) >> 11) as f64 / (1u64 << 53) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Frozen: these are names of generated worlds.
    #[test]
    fn golden_names() {
        assert_eq!(name(12345), "Thumi");
        assert_eq!(name(42), "Kibedo");
    }

    #[test]
    fn names_are_pronounceable() {
        for s in 0..1000 {
            let n = name(s);
            assert!((4..=12).contains(&n.len()) && n.chars().next().unwrap().is_uppercase(), "{n}");
        }
    }
}
