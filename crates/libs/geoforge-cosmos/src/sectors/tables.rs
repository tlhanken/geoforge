//! What kinds of star a population holds, by brightness bin: how many per M☉ and drawn
//! with what weights.
//!
//! A bright star is rare and a faint one common, so kinds are sampled by importance: masses
//! log-uniform and ages within each star's shining life, each weighted by the IMF and the
//! star formation history; a second stratum samples the dead (remnants), which carry mass
//! and, as white dwarfs, faint light. Each kind may have an unresolved companion of the same
//! age, as often as the galaxy's systems do (mass ratios flat; Duchêne & Kraus 2013).

use super::bins::{bin_of, BINS};
use crate::draw::Draw;
use crate::galaxy::Imf;
use crate::population::multiple_fraction;
use crate::sky::absolute_v;
use crate::stars::{end_of_life_gyr, evolve};
use geoforge_core::seed::Seed;
use rayon::prelude::*;

/// The Sun's absolute visual magnitude.
const SUN_V: f64 = 4.83;

/// One kind of star (a system, with its companion).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kind {
    pub initial_mass: f64,
    pub age_gyr: f64,
    /// Standard deviations off its population's metallicity relation.
    pub scatter: f64,
    pub companion: Option<f64>,
    /// Absolute V of the pair, at the population's reference metallicity.
    pub absolute_v: f64,
}

/// The kinds of one bin, with cumulative weights.
#[derive(Debug, Clone, Default)]
pub struct BinKinds {
    kinds: Vec<Kind>,
    cumulative: Vec<f64>,
}

impl BinKinds {
    fn push(&mut self, k: Kind, w: f64) {
        let total = self.cumulative.last().copied().unwrap_or(0.0) + w;
        self.kinds.push(k);
        self.cumulative.push(total);
    }

    /// The kind at `u` (0 to 1) of the bin's weight.
    #[must_use]
    pub fn pick(&self, u: f64) -> Option<&Kind> {
        let total = *self.cumulative.last()?;
        let i = self.cumulative.partition_point(|&c| c < u * total).min(self.kinds.len() - 1);
        self.kinds.get(i)
    }
}

/// A population's stars by bin.
#[derive(Debug, Clone)]
pub struct Table {
    /// Systems per M☉ of the population (present mass for the field, mass formed for an age
    /// table), by bin.
    pub per_msun: [f64; BINS],
    /// V-band light per M☉ in each bin, L☉.
    pub light: [f64; BINS],
    pub bins: Vec<BinKinds>,
}

/// The heaviest star still shining `age_gyr` after its birth, M☉.
fn heaviest_alive(age_gyr: f64, (lo, hi): (f64, f64)) -> f64 {
    if end_of_life_gyr(hi) > age_gyr {
        return hi;
    }
    let (mut a, mut b) = (lo.ln(), hi.ln());
    for _ in 0..60 {
        let mid = 0.5 * (a + b);
        if end_of_life_gyr(mid.exp()) > age_gyr {
            a = mid;
        } else {
            b = mid;
        }
    }
    a.exp()
}

/// A sample: its weight per star born, its kind, present mass and bin (if it shines).
type Sample = (f64, Kind, f64, Option<usize>);

/// One sampled star of `m` born `t` ago, with weight `w` per star born.
fn sample(d: &mut Draw, m: f64, t: f64, w: f64, feh: &(dyn Fn(f64, f64) -> f64 + Sync), lo: f64) -> Sample {
    let scatter = d.gauss();
    let z = feh(t, scatter);
    let star = evolve(m, t, z);
    let mut present = star.mass_solar;
    let mut v = absolute_v(&star);
    let mut companion = None;
    if d.u() < multiple_fraction(m) {
        let q_min = (lo / m).clamp(0.1, 1.0);
        let mc = m * d.range(q_min, 1.0);
        let mate = evolve(mc, t, z);
        present += mate.mass_solar;
        if let Some(mv) = absolute_v(&mate) {
            v = Some(v.map_or(mv, |pv| -2.5 * (10f64.powf(-0.4 * pv) + 10f64.powf(-0.4 * mv)).log10()));
        }
        companion = Some(mc);
    }
    let bin = v.map(bin_of);
    (w, Kind { initial_mass: m, age_gyr: t, scatter, companion, absolute_v: v.unwrap_or(f64::INFINITY) }, present, bin)
}

impl Table {
    /// The stars of a population born between `ages` (Gyr ago) at lookback density
    /// `density` (per Gyr, normalised over those ages), [Fe/H] `feh(age, scatter)`: `alive`
    /// samples of shining stars and a tenth as many of remnants. Counts are per M☉ of the
    /// population as it is now.
    #[must_use]
    pub fn field(imf: &Imf, ages: (f64, f64), density: &(dyn Fn(f64) -> f64 + Sync), feh: &(dyn Fn(f64, f64) -> f64 + Sync), seed: Seed, alive: usize) -> Self {
        let (lo, hi) = imf.range();
        let (a0, a1) = ages;
        let dead = alive / 10;
        let samples: Vec<Sample> = (0..alive + dead)
            .into_par_iter()
            .filter_map(|i| {
                let mut d = Draw::new(seed.child_u64(b"kind", i as u64));
                if i < alive {
                    let span = (hi / lo).ln();
                    let m = lo * (d.u() * span).exp();
                    let top = end_of_life_gyr(m).min(a1);
                    if top <= a0 {
                        return None;
                    }
                    let t = a0 + d.u() * (top - a0);
                    let w = imf.density(m) * m * span * density(t) * (top - a0) / alive as f64;
                    Some(sample(&mut d, m, t, w, feh, lo))
                } else {
                    // The dead: born heavy enough to have died within the ages.
                    let m_lo = heaviest_alive(a1, (lo, hi));
                    if m_lo >= hi {
                        return None;
                    }
                    let span = (hi / m_lo).ln();
                    let m = m_lo * (d.u() * span).exp();
                    let from = end_of_life_gyr(m).max(a0);
                    if from >= a1 {
                        return None;
                    }
                    let t = from + d.u() * (a1 - from);
                    let w = imf.density(m) * m * span * density(t) * (a1 - from) / dead.max(1) as f64;
                    Some(sample(&mut d, m, t, w, feh, lo))
                }
            })
            .collect();
        let present: f64 = samples.iter().map(|s| s.0 * s.2).sum();
        Self::from_samples(&samples, present.max(1e-300))
    }

    /// The stars of a group born together `age_gyr` ago with [Fe/H] `feh(scatter)`: counts per
    /// M☉ formed.
    #[must_use]
    pub fn at_age(imf: &Imf, age_gyr: f64, feh: &(dyn Fn(f64, f64) -> f64 + Sync), seed: Seed, n: usize) -> Self {
        let (lo, hi) = imf.range();
        let top = heaviest_alive(age_gyr, (lo, hi));
        let span = (top / lo).ln();
        let samples: Vec<Sample> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut d = Draw::new(seed.child_u64(b"kind", i as u64));
                let m = lo * (d.u() * span).exp();
                let w = imf.density(m) * m * span / n as f64;
                sample(&mut d, m, age_gyr, w, feh, lo)
            })
            .collect();
        Self::from_samples(&samples, imf.mean_mass())
    }

    fn from_samples(samples: &[Sample], per: f64) -> Self {
        let mut t = Self { per_msun: [0.0; BINS], light: [0.0; BINS], bins: vec![BinKinds::default(); BINS] };
        for &(w, k, _, bin) in samples {
            if let Some(b) = bin {
                if w > 0.0 {
                    t.per_msun[b] += w / per;
                    t.light[b] += w / per * 10f64.powf(-0.4 * (k.absolute_v - SUN_V));
                    t.bins[b].push(k, w);
                }
            }
        }
        t
    }

    /// Systems per M☉ in all bins.
    #[must_use]
    pub fn total_per_msun(&self) -> f64 {
        self.per_msun.iter().sum()
    }

    /// A kind across all bins at `u` and `v` (0 to 1): a bin by count, then a kind in it.
    #[must_use]
    pub fn pick_any(&self, u: f64, v: f64) -> Option<(usize, &Kind)> {
        let total = self.total_per_msun();
        let mut left = u * total;
        let b = self.per_msun.iter().position(|&n| {
            left -= n;
            left <= 0.0
        })?;
        self.bins[b].pick(v).map(|k| (b, k))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::galaxy::StarFormation;

    /// A disk population like the Sun's neighbourhood: about two systems per M☉, most of them
    /// red dwarfs too faint to see, a few per thousand brighter than the Sun.
    #[test]
    fn a_neighbourhood_population() {
        let imf = Imf::kroupa();
        let h = StarFormation::between(11.0, 0.3);
        let t = Table::field(&imf, (0.3, 11.0), &|a| h.lookback_density(a) / h.share_younger(11.0).max(1e-9), &|a, s| -0.02 * a + 0.2 * s, Seed(1), 20_000);
        let total = t.total_per_msun();
        assert!((1.0..3.5).contains(&total), "{total} systems per M☉");
        let brighter_than_sun: f64 = t.per_msun[..16].iter().sum();
        assert!((0.003..0.05).contains(&(brighter_than_sun / total)), "{}", brighter_than_sun / total);
        let faint: f64 = t.per_msun[20..].iter().sum();
        assert!(faint > 0.6 * total, "{}", faint / total);
        // A young group holds blue stars an old population does not.
        let young = Table::at_age(&imf, 0.01, &|_, _| 0.0, Seed(2), 8000);
        assert!(young.per_msun[..6].iter().sum::<f64>() > 0.0 && t.per_msun[..6].iter().sum::<f64>() < 1e-6);
        assert!(young.pick_any(0.5, 0.5).is_some());
    }
}
