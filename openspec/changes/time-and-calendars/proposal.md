# time-and-calendars

## Why

Time is a bare `f64` of Earth days in every API (`from_day`, `t_days`), durations carry
their unit only in their names (`rotation_h`, `year_days`, ages in Ma), and the one
calendar (`ephemeris::Calendar`) knows local days and years but only one kind of month.
The almanac, droughts, weather and sky all build on time, and worldbuilders need their own
calendars: lunar months, twelfths of the year, or a fantasy calendar with named months.

## What Changes

- An `Instant` and a `Span` in SI seconds (`geoforge-core::time`), with Earth-day and hour
  constructors, as the one representation of time under the hood.
- A world's `Clock`: its local day (mean solar day) and local year (one orbit) as constant
  spans. The second, the local day and the local year never change with the calendar.
- A `Calendar` trait that only divides a clock's years into named months, with three
  instantiations: `Twelfths` (twelve equal months), `Lunar` (a moon's lunations, new moon
  to new moon) and `Custom` (named months of whole local days, read from a world's
  `calendar.json`, for fantasy calendars).
- Dates in whole days from the epoch, each in the year it begins in (leap years fall out).
- Climate statistics (drought months, seasons) use twelfths of the local year; dates shown
  to people use the chosen calendar.
- The almanac and its command take a calendar; dates print in it.
- Deep time (rock and crust ages, plate motion) stays in Ma.

## Capabilities

- `time` — instants, spans, a world's clock, calendars

## Impact

`geoforge-core` (new `time` module), `geoforge-almanac` (`ephemeris::Calendar` built on
the clock; the almanac and droughts on `Instant` and `Calendar`), `geoforge-store` (reads
`calendar.json`), `geoforge-cli` (`almanac --calendar`). The weather, storms and sky
internals move to `Instant` in later tasks.
