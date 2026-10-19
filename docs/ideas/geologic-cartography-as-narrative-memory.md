# Geologic Cartography as Narrative Memory

## Status
Idea — not yet designed or scheduled.

## Core Concept
A semantic interpretation layer that reads Geoforge's simulated geologic history and translates it into the mythology the people of that world would tell. Not randomly generated lore bolted onto a heightmap — a derivation from the same causal chain that produced the terrain.

## The Insight
Real mythologies often encode geologic reality as narrative:
- Norse creation: Ymir's body becomes the world — geology rendered as myth.
- Aboriginal songlines: terrain IS the song; land and story are inseparable.

This layer would discover the generative rule underneath those patterns rather than inventing fiction. The myths would be load-bearing — they'd encode real information about the world's geologic history in a form a person in that world could carry and transmit.

## Mapping Examples (Geologic Event -> Archetypal Narrative)
- **Caldera** -> creation myth: the mountain that gave birth to the world, or the site where a god died.
- **Fault line** -> betrayal or unhealing wound: the land itself remembers the fracture.
- **Glacial erratic** -> giant or exile tale: a stone that traveled, carrying the memory of a distant homeland.
- **River meander migration** -> story of indecision: a spirit that couldn't settle, carving its name over a thousand years.
- **Unconformity (missing geologic time)** -> gap myth: the age the gods destroyed, the time before time — known but unnamable.

## Architectural Shape
Geoforge already simulates the geologic processes. This would add:
1. A semantic interpretation layer mapping geologic events to archetypal narrative structures.
2. An output alongside terrain: the mythopoetic record derived from the same simulation.

The causality runs one direction: terrain -> geologic history -> the stories the land would tell -> the culture that forms around those stories.

## Philosophical Claim
Mythology is a lossy compression of environmental knowledge. A sufficiently detailed geologic simulation contains the information a culture would naturally compress into myth. The story is already in the rock — this layer just makes it readable.

## Open Questions
- Where in the pipeline does this layer sit? Post-terrain generation, or integrated into the province system?
- How granular should the event-to-myth mapping be?
- Should myths reference each other (a fault-line myth that mentions the caldera creation myth), forming an interconnected mythopoetic web?
- Output format: structured data, prose, or both?

## Origin
Emerged from a conversation about cross-domain pattern recognition between procedural world generation and narrative structure. June 2026.
