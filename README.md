# Julia_ELFIN_Tools
A package for processing and interacting with data from the UCLA Electron Losses and Fields Investigation (ELFIN) CubeSats.

## Installation
download zip
unpack
place wherever you like

(NEED PYTHON 3.9 TO RUN DATA DOWNLOAD SCRIPT. ALSO IT DEFAULTS TO USING 8 THREADS NEED TO CHANGE that at top of script if not wanted)

## Example Usage
```
include("<path to Julia_ELFIN_Tools>/Events.jl")
include("<path to Julia_ELFIN_Tools>/Visualization.jl")
```

Creating an event ... creates a variable of type `Event`, which is the container type for ELFIN data recorded during the specified time period.

Integrate flux

visualize

## Description of `Event` Type
description of type with all the fields

## Provided Methods
every user-facing function and its purpose and usage

### Event Creation Methods
`create_event(DateTime start, DateTime stop, String sat; Bool warn = false)` \
Creates an Event object from a start datetime, stop datetime, and satellite ID (either "a", "A", "b", or "B")
* `start`: DateTime of the desired event start time. If the satellite was not active at this time, the nearest active time is selected as the start time.
* `stop`: DateTime of the desired event start time. If the satellite was not active at this time, the nearest active time is selected as the start time.
* `sat`: String determining whether to use ELFIN-A or ELFIN-B data. Allowed inputs are `"A"`, `"a"`, `"B"`, or `"b"`
* `warn`: Optional keyword. If `true`, print non-fatal issues to terminal for debugging purposes.

TODO

### Data Processing Methods
`integrate_flux` TODOOO

### Event Visualization Methods
---
Julia Claxton (she/her) \
Contact: julia.claxton@colorado.edu \
Lightning, Atmosphere, Ionosphere, and Radiation Belts (LAIR) Group, University of Colorado Boulder

Developed for Python 3.9.18 and Julia 1.9.2