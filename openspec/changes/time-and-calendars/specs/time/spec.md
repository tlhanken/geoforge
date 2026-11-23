## ADDED Requirements

### Requirement: One time under every calendar
Time MUST be represented as SI seconds: an `Instant` from the epoch (the focus world's
northern spring equinox) and a `Span` for lengths. Earth days and hours MUST be conversions
of seconds, never a separate representation.

#### Scenario: Earth days convert exactly
- **WHEN** an instant is built from 1.5 Earth days
- **THEN** it holds 129,600 seconds, and reads back as 1.5 Earth days

### Requirement: The clock's lengths are constants
A world's local day (its mean solar day) and local year (its orbital period) MUST be fixed
by its body and orbit, and MUST NOT depend on the calendar used to name dates.

#### Scenario: calendars agree on the year
- **WHEN** the same instant is dated in two calendars of one world
- **THEN** both give the same year and the same local day of the year

#### Scenario: a locked world has no local day
- **WHEN** the body is tidally locked to its star
- **THEN** its clock has no local day, and dates count Earth days

### Requirement: Calendars divide years into months
A calendar MUST divide each local year into months whose lengths sum to the year, and MUST
map an instant to a date (year, month, day of month, fraction of the day) and a date back to
the instant it begins.

#### Scenario: twelfths
- **WHEN** a year is divided by `Twelfths`
- **THEN** it has twelve months of a twelfth of the year each

#### Scenario: lunar months
- **WHEN** a year is divided by `Lunar` for a moon whose synodic month is longer than the
  day and shorter than the year
- **THEN** each month after month 0 begins at a new moon, and month 0 holds the days before
  the year's first new moon, if any begin there

#### Scenario: a moon that does not fit
- **WHEN** the moon's month is longer than the year or no longer than the day
- **THEN** the year is divided into twelfths

#### Scenario: a fantasy calendar
- **WHEN** a world's `calendar.json` names its months and their lengths in whole local days
- **THEN** the year is divided into those months, the last taking the days left over

#### Scenario: a month of part days
- **WHEN** a month before the last is not a whole number of days, or they fill the year
- **THEN** the calendar is refused

#### Scenario: round trip
- **WHEN** an instant is dated and the date's start is taken back
- **THEN** the start is at or before the instant and less than one local day earlier

### Requirement: Dates are whole days from the epoch
Days MUST be counted from the epoch at the same instant everywhere, and a day MUST belong
to the year and month in which it begins, so the time of day keeps step with the sun.

#### Scenario: leap years
- **WHEN** Earth's years are dated in whole days
- **THEN** each has 365 or 366 days, and eight of them hold eight orbits' days to within one

### Requirement: Climate statistics in twelfths
Statistics over months that follow the seasons (drought indices, monthly normals) MUST use
twelfths of the local year, whatever calendar dates are shown in.

#### Scenario: drought months
- **WHEN** the almanac is printed in a lunar calendar
- **THEN** its drought index is still read over twelfths of the year
