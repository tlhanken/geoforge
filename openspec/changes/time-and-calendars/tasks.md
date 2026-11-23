# Tasks — time-and-calendars

## 1. Types
- [x] 1.1 `Instant`, `Span` in SI seconds with Earth-day and hour constructors
- [x] 1.2 `Clock` from a body: local day, local year, year and day of an instant
- [x] 1.3 `Calendar` trait; `Twelfths`, `Lunar`, `Custom`; dates both ways

## 2. Use
- [x] 2.1 `ephemeris::Calendar` built on the clock
- [x] 2.2 The almanac on `Instant` and a `Calendar`; drought months as twelfths
- [x] 2.3 `calendar.json` read by the store; `almanac --calendar twelfths|lunar|world`

## 3. Migrate the internals
- [ ] 3.1 Weather and storms take `Instant` and `Span` instead of Earth days
- [ ] 3.2 Sky, ephemeris and eclipses take `Instant`
- [ ] 3.3 Commands take dates in a calendar (`--date 3/5/12`) as well as Earth days
