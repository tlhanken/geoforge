//! Galaxy models: where a galaxy's stars are, when they were born, how metal-rich.
//!
//! A galaxy is a list of components ([`Component`]): a nucleus, a bulge or bar, a thin disk
//! whose stars thicken with age, a thick disk, a halo, or one spheroid for an elliptical.
//! Each has a mass, a shape ([`shape`]), a star formation history and a metallicity relation;
//! disks may have spiral arms ([`arms`]) that gather their young stars, and dust. Stars of
//! every component are born by the same initial mass function. Morphologies are presets of
//! components ([`presets`]): spirals Sa to Sd with or without a bar, lenticular, elliptical,
//! dwarf spheroidal, dwarf irregular and starburst; [`Galaxy::milky_way`] is the default.
//! A galaxy is a plain document, so a world can carry its own.

pub mod arms;
pub mod presets;
pub mod shape;

pub use arms::Arms;
pub use presets::{GalaxyPreset, Morphology, SpiralStage};
pub use shape::{Cell, Height, PlaneLeaf, Prepared, Shape};

use crate::draw::Draw;
use serde::{Deserialize, Serialize};

/// The stars of a galaxy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Galaxy {
    pub name: String,
    pub morphology: Morphology,
    pub imf: Imf,
    pub components: Vec<Component>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub arms: Option<Arms>,
    pub dust: Dust,
    /// Where the metallicity relations take their reference values (the Sun's orbit, for the
    /// Milky Way), and near which Sun-like systems are placed, kpc.
    pub reference_kpc: f64,
    /// Disks end here, kpc.
    pub edge_kpc: f64,
}

/// One part of a galaxy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Component {
    pub kind: ComponentKind,
    /// Stellar mass now (stars and remnants), M☉.
    pub mass_msun: f64,
    pub shape: Shape,
    pub history: StarFormation,
    pub metals: Metallicity,
    /// How strongly the arms gather its youngest stars (0 to 1); older stars feel them less.
    #[serde(default)]
    pub arm_contrast: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ComponentKind {
    Nucleus,
    Bulge,
    ThinDisk,
    ThickDisk,
    Halo,
    Spheroid,
}

/// Dust: a smooth layer, and clouds where stars form.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Dust {
    /// Dimming in the plane at the reference radius from the smooth layer, mag per kpc.
    pub mag_per_kpc: f64,
    pub height_pc: f64,
    pub scale_length_kpc: f64,
    /// Gas in clouds per M☉ of stars younger than 300 Myr.
    pub cloud_gas: f64,
}

impl Dust {
    /// No dust.
    #[must_use]
    pub fn none() -> Self {
        Self { mag_per_kpc: 0.0, height_pc: 100.0, scale_length_kpc: 3.0, cloud_gas: 0.0 }
    }
}

/// How the arms' pull on a component's stars falls with age: a share that stays, and the
/// rest fading over this many Gyr (spiral structure is younger stars' and gas's).
const ARM_MEMORY_GYR: f64 = 0.5;
const ARM_LASTING: f64 = 0.15;

impl Galaxy {
    /// The arms' contrast for stars of `component` born `age_gyr` ago.
    #[must_use]
    pub fn arm_contrast(&self, component: usize, age_gyr: f64) -> f64 {
        let c = self.components[component].arm_contrast;
        if self.arms.is_none() || c <= 0.0 {
            return 0.0;
        }
        c * (ARM_LASTING + (1.0 - ARM_LASTING) * (-age_gyr / ARM_MEMORY_GYR).exp())
    }

    /// How the arms scale the density of `component`'s stars of age `age_gyr` at (`x`, `y`).
    #[must_use]
    pub fn arm_factor(&self, component: usize, age_gyr: f64, x: f64, y: f64) -> f64 {
        self.arms.as_ref().map_or(1.0, |a| a.factor(x, y, self.arm_contrast(component, age_gyr)))
    }

    /// [Fe/H] of a star of `component` born `age_gyr` ago at `p`, `scatter` standard
    /// deviations off the relation.
    #[must_use]
    pub fn metallicity(&self, component: usize, p: [f64; 3], age_gyr: f64, scatter: f64) -> f64 {
        let m = &self.components[component].metals;
        let r = p[0].hypot(p[1]) / 1000.0;
        (m.at_reference_dex + m.gradient_dex_per_kpc * (r - self.reference_kpc) + m.per_gyr_dex * age_gyr + m.scatter_dex * scatter).clamp(-3.0, 0.6)
    }

    /// The total stellar mass, M☉.
    #[must_use]
    pub fn mass_msun(&self) -> f64 {
        self.components.iter().map(|c| c.mass_msun).sum()
    }

    /// A star drawn from the whole galaxy: its component (by mass), age and place.
    pub(crate) fn draw_star(&self, d: &mut Draw) -> (usize, f64, [f64; 3]) {
        let total = self.mass_msun();
        let mut u = d.u() * total;
        let c = self.components.iter().position(|c| {
            u -= c.mass_msun;
            u <= 0.0
        });
        let c = c.unwrap_or(self.components.len() - 1);
        let comp = &self.components[c];
        let age = comp.history.draw_lookback(d.u(), d.u());
        let shape = Prepared::new(&comp.shape, age, self.edge_kpc);
        let peak = self.arms.as_ref().map_or(1.0, |_| 1.0 + 8.0 * self.arm_contrast(c, age));
        let p = shape.draw(d, &|p| self.arm_factor(c, age, p[0], p[1]) / peak);
        (c, age, p)
    }
}

impl Default for Galaxy {
    fn default() -> Self {
        Self::milky_way()
    }
}

/// Initial mass function: a broken power law dN/dm ∝ m^−α, continuous at the breaks.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Imf {
    /// Segment edges, M☉, rising: one more than `slopes`.
    pub breaks: Vec<f64>,
    /// α of each segment.
    pub slopes: Vec<f64>,
}

impl Imf {
    /// Kroupa (2001) over hydrogen-burning stars: α = 1.3 from 0.08 to 0.5 M☉, 2.3 above, to
    /// 150 M☉. Brown dwarfs (α = 0.3 below 0.08) are not drawn.
    #[must_use]
    pub fn kroupa() -> Self {
        Self { breaks: vec![0.08, 0.5, 150.0], slopes: vec![1.3, 2.3] }
    }

    /// A galaxy's field stars taken together: Kroupa (2001) to 1 M☉, then α = 2.6, near the
    /// field's present-day slope of 2.5–2.7 (Scalo 1986). Clusters are born by the canonical 2.3, but
    /// most stars form in small groups that cannot make the heaviest, so a galaxy integrated
    /// over its groups is steeper above a solar mass (Kroupa & Weidner 2003).
    #[must_use]
    pub fn galactic() -> Self {
        Self { breaks: vec![0.08, 0.5, 1.0, 150.0], slopes: vec![1.3, 2.3, 2.6] }
    }

    /// Lightest and heaviest masses born, M☉.
    #[must_use]
    pub fn range(&self) -> (f64, f64) {
        (self.breaks[0], self.breaks[self.breaks.len() - 1])
    }

    /// Probability density of initial mass `m` (per M☉), zero outside the breaks.
    #[must_use]
    pub fn density(&self, m: f64) -> f64 {
        let (mut k, mut total, mut at) = (1.0, 0.0, 0.0);
        for (i, &alpha) in self.slopes.iter().enumerate() {
            if i > 0 {
                k *= self.breaks[i].powf(alpha - self.slopes[i - 1]);
            }
            let (lo, hi) = (self.breaks[i], self.breaks[i + 1]);
            total += k * power_integral(lo, hi, alpha);
            if (lo..=hi).contains(&m) && at == 0.0 {
                at = k * m.powf(-alpha);
            }
        }
        at / total
    }

    /// Mean initial mass of a star, M☉.
    #[must_use]
    pub fn mean_mass(&self) -> f64 {
        let (lo, hi) = self.range();
        let n = 2000;
        let step = (hi / lo).ln() / f64::from(n);
        (0..n).map(|i| (lo.ln() + (f64::from(i) + 0.5) * step).exp()).map(|m| self.density(m) * m * m * step).sum()
    }
}

/// ∫ m^(−α) dm from `lo` to `hi`.
pub(crate) fn power_integral(lo: f64, hi: f64, alpha: f64) -> f64 {
    if (alpha - 1.0).abs() < 1e-9 {
        (hi / lo).ln()
    } else {
        (hi.powf(1.0 - alpha) - lo.powf(1.0 - alpha)) / (1.0 - alpha)
    }
}

/// When a component's stars were born.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StarFormation {
    /// Age of the oldest stars, Gyr.
    pub age_gyr: f64,
    /// Age of the youngest, Gyr: 0 while stars still form.
    #[serde(default)]
    pub youngest_gyr: f64,
    /// E-folding time of the star formation rate's decline, Gyr: small means most stars are
    /// old; negative means the rate has been rising. `None`: a constant rate.
    pub decline_gyr: Option<f64>,
    /// Share of all stars born in a burst ending with the youngest, and the burst's length,
    /// Gyr.
    pub burst_share: f64,
    pub burst_gyr: f64,
}

impl StarFormation {
    /// A steady rate from `age_gyr` ago to `youngest_gyr` ago.
    #[must_use]
    pub fn between(age_gyr: f64, youngest_gyr: f64) -> Self {
        Self { age_gyr, youngest_gyr, decline_gyr: None, burst_share: 0.0, burst_gyr: 0.1 }
    }

    /// Probability density of a star's birth `lookback_gyr` ago (per Gyr): the burst's share
    /// spread evenly over it, the rest following the rate's decline.
    #[must_use]
    pub fn lookback_density(&self, lookback_gyr: f64) -> f64 {
        let (y, a) = (self.youngest_gyr, self.age_gyr);
        if !(y..=a).contains(&lookback_gyr) {
            return 0.0;
        }
        let burst = if lookback_gyr < y + self.burst_gyr { self.burst_share / self.burst_gyr.min(a - y) } else { 0.0 };
        let steady = match self.decline_gyr {
            // Rate ∝ e^(−t/τ) since the first stars: lookback density ∝ e^(lookback/τ).
            Some(tau) if tau.abs() > 1e-6 => {
                let k = ((a - y) / tau).clamp(-700.0, 700.0);
                ((lookback_gyr - y) / tau).clamp(-700.0, 700.0).exp() / (tau * k.exp_m1())
            }
            _ => 1.0 / (a - y),
        };
        burst + (1.0 - self.burst_share) * steady
    }

    /// The share of stars born less than `lookback_gyr` ago.
    #[must_use]
    pub fn share_younger(&self, lookback_gyr: f64) -> f64 {
        let (y, a) = (self.youngest_gyr, self.age_gyr);
        let t = lookback_gyr.clamp(y, a);
        let burst = self.burst_share * ((t - y) / self.burst_gyr.min(a - y)).min(1.0);
        let steady = match self.decline_gyr {
            Some(tau) if tau.abs() > 1e-6 => (((t - y) / tau).clamp(-700.0, 700.0).exp_m1()) / (((a - y) / tau).clamp(-700.0, 700.0).exp_m1()),
            _ => (t - y) / (a - y),
        };
        burst + (1.0 - self.burst_share) * steady
    }

    /// A birth lookback drawn from two uniform numbers, Gyr.
    #[must_use]
    pub fn draw_lookback(&self, u: f64, v: f64) -> f64 {
        let (y, a) = (self.youngest_gyr, self.age_gyr);
        if v < self.burst_share {
            return y + u * self.burst_gyr.min(a - y);
        }
        match self.decline_gyr {
            Some(tau) if tau.abs() > 1e-6 => {
                let k = ((a - y) / tau).clamp(-700.0, 700.0);
                (y + tau * (1.0 + u * k.exp_m1()).ln()).clamp(y, a)
            }
            _ => y + u * (a - y),
        }
    }
}

/// How metal-rich a star is born: `[Fe/H] = at_reference + gradient·(R − R_ref) +
/// per_gyr·age`, scattered by `scatter_dex`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Metallicity {
    /// [Fe/H] of stars born now at the reference radius, dex.
    pub at_reference_dex: f64,
    /// Radial gradient, dex per kpc.
    pub gradient_dex_per_kpc: f64,
    /// Change per Gyr of age (negative: older stars are poorer), dex.
    pub per_gyr_dex: f64,
    /// Scatter at a given age and radius, dex.
    pub scatter_dex: f64,
}

impl Metallicity {
    /// A constant [Fe/H] with a scatter.
    #[must_use]
    pub fn flat(dex: f64, scatter_dex: f64) -> Self {
        Self { at_reference_dex: dex, gradient_dex_per_kpc: 0.0, per_gyr_dex: 0.0, scatter_dex }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The IMF and every component's history integrate to one, and draws follow them.
    #[test]
    fn densities_are_normalised() {
        let imf = Imf::kroupa();
        let n = 20_000;
        let (lo, hi) = (0.08f64.ln(), 150f64.ln());
        let step = (hi - lo) / f64::from(n);
        let mass: f64 = (0..n).map(|i| (lo + (f64::from(i) + 0.5) * step).exp()).map(|m| imf.density(m) * m * step).sum();
        assert!((mass - 1.0).abs() < 1e-3, "{mass}");
        for preset in GalaxyPreset::ALL {
            for c in preset.galaxy().components {
                let h = c.history;
                let span = h.age_gyr - h.youngest_gyr;
                let total: f64 = (0..n).map(|i| h.lookback_density(h.youngest_gyr + (f64::from(i) + 0.5) / f64::from(n) * span) * span / f64::from(n)).sum();
                assert!((total - 1.0).abs() < 1e-2, "{} {:?} {total}", preset.name(), c.kind);
                assert!((h.share_younger(h.age_gyr) - 1.0).abs() < 1e-9 && h.share_younger(h.youngest_gyr).abs() < 1e-9);
                let mut d = Draw::new(geoforge_core::seed::Seed(2));
                let mid = h.youngest_gyr + 0.5 * span;
                let drawn = (0..4000).filter(|_| h.draw_lookback(d.u(), d.u()) < mid).count() as f64 / 4000.0;
                assert!((drawn - h.share_younger(mid)).abs() < 0.03, "{} {:?}: {drawn} vs {}", preset.name(), c.kind, h.share_younger(mid));
            }
        }
        assert!((imf.mean_mass() - 0.6).abs() < 0.1, "{}", imf.mean_mass());
    }

    /// The Milky Way's arms gather young thin-disk stars, not old ones or the bulge's.
    #[test]
    fn arms_pull_young_stars() {
        let g = Galaxy::milky_way();
        let thin = g.components.iter().position(|c| c.kind == ComponentKind::ThinDisk).unwrap();
        let bulge = g.components.iter().position(|c| c.kind == ComponentKind::Bulge).unwrap();
        assert!(g.arm_contrast(thin, 0.01) > 0.6 && g.arm_contrast(thin, 8.0) < 0.2 && g.arm_contrast(bulge, 0.01) == 0.0);
    }
}
