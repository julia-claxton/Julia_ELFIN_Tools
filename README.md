# Julia_ELFIN_Tools
A package for processing and interacting with data from the [UCLA Electron Losses and Fields Investigation (ELFIN)](https://elfin.igpp.ucla.edu/) CubeSats.

## Installation
To use this library, download the .zip from GitHub and unpack it. This should give you a folder named `Julia_ELFIN_Tools`. Move this folder to a location of your convenience. Add these lines to any Julia script to access the functions:

```
include("<path to Julia_ELFIN_Tools folder>/Events.jl")
include("<path to Julia_ELFIN_Tools folder>/Visualization.jl")
```

This library requires ELFIN data to work. ELFIN data is not included in the download of this library due to its size (approximately 15 GB). To download this data yourself, follow these instructions. If these instructions do not work for you, please email me (julia.claxton@colorado.edu) for alternatives.
1. In a terminal, navigate to the `Julia_ELFIN_Tools` folder.
2. Create a python virtual environment with `python3.9 -m venv ./python_environment`. **Important: This requires Python 3.9 in order for `spacepy` (a required library) to work!**
3. Activate the virtual environment with `source ./python_environment/bin/activate`
4. Install required Python packages with `pip install -r requirements.txt`
5. Run the ELFIN data download script with `python3.9 elfin_data_download.py`. This script uses 8 threads by default. This can be changed by editing the `N_THREADS` variable at the top of `elfin_data_download.py`.

This is the only part of the library that uses Python. After this step, Julia is used exclusively.

## Example Usage
This library's main function is providing the `Event` type, a container for ELFIN data. The `Event` type is explained in detail in later sections. For now, we will just look at the basics of creating an event and viewing it. First, add the library's functions to your Julia script with the following commands:

```
include("<path to Julia_ELFIN_Tools>/Events.jl")
include("<path to Julia_ELFIN_Tools>/Visualization.jl")
```

Next, we will create an event. All you need is the date and times of the period you wish to view ELFIN data for.

```
start_time = DateTime("2020-09-02T14:21:00")
stop_time = DateTime("2020-09-02T14:25:00")
satellite = "B"

my_event = create_event(start_time, stop_time, satellite)
```

We can now view our event:

```
quicklook(my_event)
```
![image](./readme_files/quicklook.png)

Hooray!

## Fields of the `Event` Type
The Event type contains several fields with metadata, location/time, and science data. It serves as a container for the data recorded by ELFIN during a given period of time. Below is a description of the fields provided by the `Event` type.

### Metadata Fields
| Field Name | Description | Units |
|---:|:---|:---:|
| `satellite` | Which satellite the data was recorded on, ELFIN-A or -B. Can be either "A" or "B". |  |
| `name` | Single string that can be used to re-create an event. Format = `"yyyy-mm-dd_HH:MM:SS_<duration>_EL<satellite>"` |  |
| `duration` | Event duration in seconds | s |
| `n_datapoints` | Number of datapoints recorded during the event |  |
| `n_observations` | Number of distinct periods where instruments were recording in the event |  |
| `observation_start_idxs` | Indices where observations periods begin |  |
| `observation_stop_idxs` | Indices where observations periods end |  |
| `kp` | 3-hour Kp index for each datapoint |  |
| `dst` | Hourly Dst index for each datapoint | nT |

### Time & Location Fields
| Field Name | Description | Units |
|---:|:---|:---:|
| `time` | Time since event start for each datapoint | s |
| `time_datetime` | Date and time of each datapoint as DateTime objects |  |
| `position` | ELFIN's position at each time point in GEI coordinates. Access coordinates with `position["x"][t]`, `position["y"][t]`, `position["z"][t]` | km |
| `L` | L value of ELFIN's location at each data point. Calculated with T87 model |  |
| `MLT` | MLT of ELFIN's location at each data point | hr |

### Science Data Fields
| Field Name | Description | Units |
|---:|---|:---:|
| `pitch_angles` | Pitch angle of each of ELFIN's look directions at each timestep. Dimensions = [time, look direction] | degrees |
| `energy_bins_min` | Minimum value of each of ELFIN's energy channels | keV |
| `energy_bins_mean` | Mean value of each of ELFIN's energy channels | keV |
| `energy_bins_max` | Maximum value of each of ELFIN's energy channels | keV |
| `lc_idxs` | Look direction indices inside the loss cone at each time step |  |
| `alc_idxs` | Look direction indices inside the anti loss cone at each time step |  |
| `loss_cone_angles` | Angle of the loss cone at each time step | degrees |
| `anti_loss_cone_angles` | Angle of the anti-loss cone at each time step | degrees |
| `avg_pitch_angles` | Pitch angle of each of ELFIN's look directions averaged over all timesteps. Mainly useful for short events. | degrees |
| `avg_loss_cone_angle` | Average loss cone angle over full event. | degrees |
| `avg_anti_loss_cone_angle` | Average anti loss cone angle over full event. | degrees |
| `e_flux` | Electron energy flux at each timestep, energy channel, and look direction. Dimensions = [time, energy channel, pitch angle] | kev/(cm^2 str s MeV) |
| `n_flux` | Electron number flux at each timestep, energy channel, and look direction. Dimensions = [time, energy channel, pitch angle] | #/(cm^2 str s MeV) |
| `Jprec_over_Jperp` | Ratio of electron flux in the loss cone to flux in the trapped area at each time step and energy channel. Dimensions = [time, energy channel] |  |


## Provided Methods
This library also provides a number of methods for interacting with the `Event` type. Methods are provided to create events, integrate the flux data, and visualize the data.

### Event Creation Methods
`create_event(DateTime start, DateTime stop, String sat; Bool warn = false)` \
Creates an Event object from a start datetime, stop datetime, and satellite ID (either "a", "A", "b", or "B")

Arguments:
* `start`: DateTime of the desired event start time. If the satellite was not active at this time, the nearest active time is selected as the start time.

* `stop`: DateTime of the desired event stop time. If the satellite was not active at this time, the nearest active time is selected as the stop time.

* `sat`: String determining whether to use ELFIN-A or ELFIN-B data. Allowed inputs are `"A"`, `"a"`, `"B"`, or `"b"`.

* `warn`: Optional keyword. If `true`, print non-fatal issues to terminal for debugging purposes.

Returns:
* Variable of type `Event` containing ELFIN data between specified times.

Example Usage:
```
event = create_event(DateTime("2020-09-02T14:21:00"), DateTime("2020-09-02T14:25:00"), "b")
```
---
<br>

`create_event(Date date, String sat; Bool warn = false)` \
Creates an Event object containing all data recorded by an ELFIN satellite on a given date.

Arguments:
* `date`: Date to create event for.

* `sat`: String determining whether to use ELFIN-A or ELFIN-B data. Allowed inputs are `"A"`, `"a"`, `"B"`, or `"b"`

* `warn`: Optional keyword. If `true`, print non-fatal issues to terminal for debugging purposes.

Returns:
* Variable of type `Event` containing ELFIN data on the specified date.

Example Usage:
```
event = create_event(Date("2020-09-02"), "b")
```
---
<br>

`create_event(String name; Bool warn = false)` \
Creates an Event object based on the namestring of another event. You can access the name of an event using the `name` field.

Arguments:
* `name`: String generated in the `name` field of an event. Can also be user-generated using the format `"yyyy-mm-dd_HH:MM:SS_<duration>_EL<satellite>"`

* `warn`: Optional keyword. If `true`, print non-fatal issues to terminal for debugging purposes.

Returns:
* Variable of type `Event` containing ELFIN data specified by namestring.

Example Usage:
```
event = create_event("2020-09-02_14:21:00_229.138_ELB")
```
---

### Data Processing Methods
```
integrate_flux(Event event; Bool time = false, Bool pitch_angle = false, Bool energy = false, 
               UnitRange time_range = 1:event.n_datapoints,
               Tuple{Float64} energy_range_keV = (-Inf, Inf),
               String pitch_angle_range = "full",
               )
```
Integrates the flux recorded by ELFIN with respect to time, energy, pitch angle, or any combination thereof. Specify which dimensions to integrate along setting the keyword arguments `time`, `energy`, `pitch angle` to `true` according to the dimensions to integrate along.

Arguments:
* `event`: Event to integrate the flux for.

* `time`: Optional keyword argument, set to `true` to integrate flux in time. Removes the 1/s dimension from the flux data.

* `time_range`: Range of time indices to integrate flux over, if only a subset of the event is to be integrated. Defaults to integrate over the full event.

* `energy`: Optional keyword argument, set to `true` to integrate flux over energy. Removes the 1/MeV dimension from the flux data.

* `energy_range_keV`: Range of energies to integrate over, if only a subset of energies is needed. Defaults to integrating all energy channels.

* `pitch_angle`: Optional keyword argument, set to `true` to integrate flux over pitch angle. Removes the 1/str dimension from the flux data.

* `pitch_angle_range`: Region to integrate pitch angle. Set to `"loss cone"` to integrate over the loss cone, `"anti loss cone"` for the anti loss cone, `"trapped"` for the trapped region, and `"full"` (default) to integrate over all available pitch angles.

Returns:
* `e_flux`: Energy flux integrated as specified by input arguments.

* `n_flux`: Number flux integrated as specified by input arguments.

Example Usage:
```
total_energy_recorded_in_event, total_electrons_recorded_in_event = integrate_flux(event, time = true, energy = true, pitch_angle = true)

omnidirectional_energy_flux, omnidirectional_number_flux = integrate_flux(event, pitch_angle = true)
```
---

### Event Visualization Methods
`quicklook(Event event; String by = "index")`
Generates a quick-look chart showing an overview of the data recorded during the event. Chart generation can sometimes be slow due to the number of subplots.

Arguments:
* `event`: Event to generate quicklook chart for.

* `by`: For graphs in the quicklook that have a time axis, this optional keyword determines what the time axis will be. `"index"` (default) plots data by the timestep's index. `"time"` plots data by number of seconds elapsed since event start time, and `"date"` plots time series data by absolute date and time.

Returns:
* None

Example Usage:
```
event = create_event(DateTime("2020-09-02T14:21:00"), DateTime("2020-09-02T14:25:00"), "b")
quicklook(event)
```

---
Julia Claxton (she/her) \
Contact: julia.claxton@colorado.edu \
Lightning, Atmosphere, Ionosphere, and Radiation Belts (LAIR) Group, University of Colorado Boulder

Developed for Python 3.9.18 and Julia 1.9.2