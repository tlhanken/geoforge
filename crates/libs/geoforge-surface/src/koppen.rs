//! Köppen–Geiger climate classes, computed from the climate layer's products.
//!
//! The 30 classes of Beck et al. (2018), by Peel et al. (2007)'s rules. Arid (B) first, where
//! annual rain is under 20·T + 280 mm when 70% of it falls in the warm half-year, 20·T when
//! 70% falls in the cold half, 20·T + 140 otherwise (desert under half of that); then polar
//! (E) where the warmest month is under 10 °C; tropical (A) where the coldest is at least
//! 18 °C; temperate (C) or continental (D) either side of a coldest-month boundary, 0 °C by
//! default (Peel; Beck) or −3 °C (Köppen; Kottek et al. 2006).
//!
//! The climate layer stores annual means, the warmest-to-coldest month range and the rain's
//! warm-half share, not twelve months, so months are rebuilt as sinusoids in phase with each
//! other: temperature swings by the stored range, rain by whatever amplitude gives the stored
//! warm-half share. The rules that count months or need a driest or wettest month read these,
//! so a peaky continental summer or a two-season rain regime reads as a smoother one.

use geoforge_core::product::Products;
use rayon::prelude::*;
use std::f64::consts::PI;

/// A Köppen–Geiger class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[rustfmt::skip]
pub enum Koppen {
    Af, Am, Aw,
    BWh, BWk, BSh, BSk,
    Csa, Csb, Csc, Cwa, Cwb, Cwc, Cfa, Cfb, Cfc,
    Dsa, Dsb, Dsc, Dsd, Dwa, Dwb, Dwc, Dwd, Dfa, Dfb, Dfc, Dfd,
    ET, EF,
}

impl Koppen {
    #[rustfmt::skip]
    pub const ALL: [Self; 30] = [
        Self::Af, Self::Am, Self::Aw,
        Self::BWh, Self::BWk, Self::BSh, Self::BSk,
        Self::Csa, Self::Csb, Self::Csc, Self::Cwa, Self::Cwb, Self::Cwc, Self::Cfa, Self::Cfb, Self::Cfc,
        Self::Dsa, Self::Dsb, Self::Dsc, Self::Dsd, Self::Dwa, Self::Dwb, Self::Dwc, Self::Dwd,
        Self::Dfa, Self::Dfb, Self::Dfc, Self::Dfd,
        Self::ET, Self::EF,
    ];

    #[must_use]
    pub fn code(self) -> &'static str {
        self.info().0
    }

    #[must_use]
    pub fn name(self) -> &'static str {
        self.info().1
    }

    /// The colour of Beck et al. (2018)'s maps.
    #[must_use]
    pub fn color(self) -> [u8; 3] {
        self.info().2
    }

    /// The major group, with B split into BW and BS.
    #[must_use]
    pub fn group(self) -> &'static str {
        let c = self.code();
        &c[..if c.starts_with('B') { 2 } else { 1 }]
    }

    #[must_use]
    pub fn from_code(code: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|k| k.code() == code)
    }

    fn info(self) -> (&'static str, &'static str, [u8; 3]) {
        match self {
            Self::Af => ("Af", "Tropical rainforest", [0, 0, 255]),
            Self::Am => ("Am", "Tropical monsoon", [0, 120, 255]),
            Self::Aw => ("Aw", "Tropical savanna", [70, 170, 250]),
            Self::BWh => ("BWh", "Hot desert", [255, 0, 0]),
            Self::BWk => ("BWk", "Cold desert", [255, 150, 150]),
            Self::BSh => ("BSh", "Hot steppe", [245, 165, 0]),
            Self::BSk => ("BSk", "Cold steppe", [255, 220, 100]),
            Self::Csa => ("Csa", "Mediterranean, hot summer", [255, 255, 0]),
            Self::Csb => ("Csb", "Mediterranean, warm summer", [200, 200, 0]),
            Self::Csc => ("Csc", "Mediterranean, cool summer", [150, 150, 0]),
            Self::Cwa => ("Cwa", "Monsoon subtropical, hot summer", [150, 255, 150]),
            Self::Cwb => ("Cwb", "Subtropical highland", [100, 200, 100]),
            Self::Cwc => ("Cwc", "Subtropical highland, cool summer", [50, 150, 50]),
            Self::Cfa => ("Cfa", "Humid subtropical", [200, 255, 80]),
            Self::Cfb => ("Cfb", "Oceanic", [100, 255, 80]),
            Self::Cfc => ("Cfc", "Subpolar oceanic", [50, 200, 0]),
            Self::Dsa => ("Dsa", "Continental, dry hot summer", [255, 0, 255]),
            Self::Dsb => ("Dsb", "Continental, dry warm summer", [200, 0, 200]),
            Self::Dsc => ("Dsc", "Subarctic, dry summer", [150, 50, 150]),
            Self::Dsd => ("Dsd", "Subarctic, dry summer, severe winter", [150, 100, 150]),
            Self::Dwa => ("Dwa", "Continental monsoon, hot summer", [170, 175, 255]),
            Self::Dwb => ("Dwb", "Continental monsoon, warm summer", [90, 120, 220]),
            Self::Dwc => ("Dwc", "Subarctic monsoon", [75, 80, 180]),
            Self::Dwd => ("Dwd", "Subarctic monsoon, severe winter", [50, 0, 135]),
            Self::Dfa => ("Dfa", "Humid continental, hot summer", [0, 255, 255]),
            Self::Dfb => ("Dfb", "Humid continental, warm summer", [55, 200, 255]),
            Self::Dfc => ("Dfc", "Subarctic", [0, 125, 125]),
            Self::Dfd => ("Dfd", "Subarctic, severe winter", [0, 70, 95]),
            Self::ET => ("ET", "Tundra", [178, 178, 178]),
            Self::EF => ("EF", "Ice cap", [102, 102, 102]),
        }
    }
}

/// Classification settings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KoppenRules {
    /// Coldest month separating temperate (C) from continental (D), °C.
    pub cold_boundary_c: f64,
}

impl Default for KoppenRules {
    fn default() -> Self {
        Self { cold_boundary_c: 0.0 }
    }
}

/// One month's mean temperature (°C) and rain (mm).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Month {
    pub t: f64,
    pub rain: f64,
}

/// Month k's phase from the warmest month (0). As in a real calendar, the warm half-year runs
/// from three months before the warmest to two after (April to September in the north).
fn phase(k: usize) -> f64 {
    2.0 * PI * k as f64 / 12.0
}

fn is_warm(k: usize) -> bool {
    !(3..9).contains(&k)
}

fn rain_weights(amplitude: f64) -> [f64; 12] {
    std::array::from_fn(|k| (1.0 + amplitude * phase(k).cos()).max(0.0))
}

fn warm_share_of(w: &[f64; 12]) -> f64 {
    let warm: f64 = (0..12).filter(|&k| is_warm(k)).map(|k| w[k]).sum();
    warm / w.iter().sum::<f64>()
}

/// Twelve months rebuilt from a climate's annual mean `t` (°C), warmest-to-coldest `range`
/// (°C), annual `precip` (mm) and its `warm_share` in the warm half-year (months 9–11 and
/// 0–2). Month 0 is the warmest, month 6 the coldest.
#[must_use]
pub fn months(t: f64, range: f64, precip: f64, warm_share: f64) -> [Month; 12] {
    // Share rises monotonically with the rain amplitude, towards 0 and 1 as the amplitude grows
    // (all the rain in the cold or warm months). Bisect for it.
    let (mut lo, mut hi) = (-1e3, 1e3);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if warm_share_of(&rain_weights(mid)) < warm_share {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let w = rain_weights(0.5 * (lo + hi));
    let total: f64 = w.iter().sum();
    std::array::from_fn(|k| Month { t: t + range / 2.0 * phase(k).cos(), rain: precip * w[k] / total })
}

/// The class of a climate: annual mean `t` (°C), warmest-to-coldest month `range` (°C), annual
/// rain `precip` (mm) and the share of it that falls in the warm half-year.
#[must_use]
pub fn classify(t: f64, range: f64, precip: f64, warm_share: f64, rules: &KoppenRules) -> Koppen {
    use Koppen::*;
    let (hot, cold) = (t + range / 2.0, t - range / 2.0);
    let threshold = 20.0 * t
        + if warm_share >= 0.7 {
            280.0
        } else if warm_share <= 0.3 {
            0.0
        } else {
            140.0
        };
    if precip < threshold {
        return match (precip < threshold / 2.0, t >= 18.0) {
            (true, true) => BWh,
            (true, false) => BWk,
            (false, true) => BSh,
            (false, false) => BSk,
        };
    }
    if hot < 10.0 {
        return if hot > 0.0 { ET } else { EF };
    }
    let m = months(t, range, precip, warm_share);
    let rain = |warm: bool| (0..12).filter(move |&k| is_warm(k) == warm).map(|k| m[k].rain);
    let fold = |warm: bool, f: fn(f64, f64) -> f64, from: f64| rain(warm).fold(from, f);
    let (summer_dry, summer_wet) = (fold(true, f64::min, f64::MAX), fold(true, f64::max, 0.0));
    let (winter_dry, winter_wet) = (fold(false, f64::min, f64::MAX), fold(false, f64::max, 0.0));
    if cold >= 18.0 {
        let driest = summer_dry.min(winter_dry);
        return if driest >= 60.0 {
            Af
        } else if driest >= 100.0 - precip / 25.0 {
            Am
        } else {
            Aw
        };
    }
    let second = if summer_dry < 40.0 && summer_dry < winter_wet / 3.0 {
        's'
    } else if winter_dry < summer_wet / 10.0 {
        'w'
    } else {
        'f'
    };
    let first = if cold > rules.cold_boundary_c { 'C' } else { 'D' };
    let third = if hot >= 22.0 {
        'a'
    } else if m.iter().filter(|x| x.t >= 10.0).count() >= 4 {
        'b'
    } else if first == 'D' && cold < -38.0 {
        'd'
    } else {
        'c'
    };
    Koppen::from_code(&format!("{first}{second}{third}")).expect("every C and D combination is a class")
}

/// The class of every cell, `None` at sea (below 0 m, on a world with `sea`); `None` overall
/// until the climate layer has run.
#[must_use]
pub fn map(p: &Products, sea: bool, rules: &KoppenRules) -> Option<Vec<Option<Koppen>>> {
    let (t, range) = (&p.temperature()?.values, &p.temperature_range()?.values);
    let (precip, share) = (&p.precipitation()?.values, &p.precipitation_warm_share()?.values);
    let elev = p.elevation().or(p.rough_elevation()).map(|f| f.values.as_slice());
    let wet = |i: usize| sea && elev.is_some_and(|e| e[i] < 0.0);
    Some(
        (0..t.len())
            .into_par_iter()
            .map(|i| (!wet(i)).then(|| classify(f64::from(t[i]), f64::from(range[i]), f64::from(precip[i]), f64::from(share[i]), rules)))
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn k(t: f64, range: f64, p: f64, share: f64) -> &'static str {
        classify(t, range, p, share, &KoppenRules::default()).code()
    }

    /// Climates like these stations' (annual mean, range, rain, warm-half share).
    #[test]
    fn stations() {
        assert_eq!(k(26.7, 1.5, 2300.0, 0.55), "Af", "Manaus");
        assert_eq!(k(27.0, 3.0, 3000.0, 0.85), "Am", "Conakry");
        assert_eq!(k(27.0, 3.0, 1100.0, 0.85), "Aw", "a savanna");
        assert_eq!(k(25.0, 13.0, 30.0, 0.5), "BWh", "Riyadh");
        assert_eq!(k(10.0, 30.0, 120.0, 0.6), "BWk", "a basin desert");
        assert_eq!(k(13.0, 26.0, 450.0, 0.75), "BSk", "Dodge City");
        assert_eq!(k(17.5, 16.0, 450.0, 0.1), "Csa", "Athens");
        assert_eq!(k(14.0, 7.0, 500.0, 0.1), "Csb", "San Francisco");
        assert_eq!(k(11.3, 13.0, 600.0, 0.5), "Cfb", "London");
        assert_eq!(k(16.0, 22.0, 1200.0, 0.62), "Cfa", "Shanghai");
        assert_eq!(k(12.9, 30.0, 570.0, 0.85), "Dwa", "Beijing");
        assert_eq!(k(5.8, 26.0, 700.0, 0.62), "Dfb", "Moscow");
        assert_eq!(k(-4.0, 36.0, 400.0, 0.6), "Dfc", "a taiga town");
        assert_eq!(k(-14.5, 64.0, 180.0, 0.7), "Dfd", "Verkhoyansk");
        assert_eq!(k(-11.0, 31.0, 140.0, 0.65), "ET", "Utqiagvik");
        assert_eq!(k(-30.0, 30.0, 250.0, 0.5), "EF", "Summit");
    }

    #[test]
    fn months_keep_the_figures() {
        for share in [0.05, 0.3, 0.5, 0.8, 0.97] {
            let m = months(10.0, 20.0, 900.0, share);
            let total: f64 = m.iter().map(|x| x.rain).sum();
            let warm: f64 = (0..12).filter(|&k| is_warm(k)).map(|k| m[k].rain).sum();
            assert!((total - 900.0).abs() < 1e-6 && (warm / total - share).abs() < 1e-6, "{share}: {}", warm / total);
        }
        let m = months(10.0, 20.0, 900.0, 0.5);
        assert!((m[0].t - 20.0).abs() < 1e-9 && (m[6].t - 0.0).abs() < 1e-9);
    }

    #[test]
    fn cold_boundary_is_a_rule() {
        let continental = classify(8.0, 20.0, 700.0, 0.5, &KoppenRules::default());
        let temperate = classify(8.0, 20.0, 700.0, 0.5, &KoppenRules { cold_boundary_c: -3.0 });
        assert_eq!((continental.group(), temperate.group()), ("D", "C"));
    }

    #[test]
    fn codes_round_trip() {
        for k in Koppen::ALL {
            assert_eq!(Koppen::from_code(k.code()), Some(k));
        }
    }
}
