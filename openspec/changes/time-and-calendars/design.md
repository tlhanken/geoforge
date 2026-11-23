# Design — time-and-calendars

## One time, many calendars

    Instant ── SI seconds from the epoch (the focus's northern spring equinox, "the present")
    Span    ── SI seconds
    Clock   ── the world's constants: local day (mean solar day; none if locked), local year
    Calendar ─ months on top of a clock: Twelfths | Lunar(moon) | Custom(calendar.json)

`f64` seconds: relative precision 1e-16, so ±10,000 years from the epoch resolve better
than a millisecond; integer kiloseconds would add nothing but rounding. Deep time stays in
Ma (`f64`), its own unit.

## Clock

The local day is the mean solar day, `1 / (1/spin − 1/year)` (retrograde spins add), as
`ephemeris::Calendar` computes it today; a world locked to its star has none and counts
its time in Earth days. The local year is the orbital period. Year `n` starts at
`n × year` from the epoch; dates before the epoch have negative years.

## Calendars

A calendar answers one question: how is year `n` divided into months?

- `Twelfths`: twelve months of a twelfth of the year each (30.44 local days on Earth).
  What the climate statistics use, since seasons follow the orbit.
- `Lunar`: one moon's lunations, new moon to new moon, the new moons found from the
  ephemeris (where the moon's lit share is least). Month `n` begins at the year's `n`-th
  new moon; the days before its first, if any begin there, are month 0. A moon whose month
  is longer than the year, or no longer than the day (Phobos), falls back to twelfths.
- `Custom`: named months of whole local days, from `calendar.json` beside the world's
  `body.json`; the last month takes the days left over, so the year stays the clock's.
  It may name a week (counted from the epoch), the present year's number and an era.

## Dates are whole days

A date is (year, month, day of month, share of the day gone). Days are counted from the
epoch, the same instant everywhere (like UTC at the prime meridian), so the time of day
keeps step with the sun. As the year is rarely a whole number of days, a day belongs to the
year and month in which it begins: calendar years differ by a day (Earth's are 365 or 366,
its leap years fall out of the orbit), months keep their whole days, and the clock's year
stays one orbit.

## What does not change

The second, the local day and the local year: a calendar never changes how long anything
is, only how it is named.
