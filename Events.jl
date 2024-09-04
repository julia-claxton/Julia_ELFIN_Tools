# Housekeeping
const Julia_ELFIN_Tools_TOP_LEVEL = @__DIR__
@assert endswith(Julia_ELFIN_Tools_TOP_LEVEL, "Julia_ELFIN_Tools")

# Library includes
using Statistics
using LinearAlgebra
using Dates
using PyCall
using DelimitedFiles
using NumericalIntegration
using BasicInterpolators

"""
Events.jl

This file provides the type Event used for storing and representing ELFIN data.

To create an Event, use the function create_event(). This function has several methods for creating events in different ways. Refer to
the comments on each method of create_event() for the details of creating an event. 

This file also provides the method integrate_flux(), which operates on Events. It integrates the flux data stored in the event object
in time, energy, pitch angle, or any combination of the three.

Written by Julia Claxton. Contact: julia.claxton@colorado.edu
Released under MIT License (see LICENSE.txt for full license)
"""

struct Event
#============== FIELD NAME ===================|=========================== FIELD DESCRIPTION ==================================|============= FIELD UNITS =============#
    # Metadata fields         
    satellite::String                         # Which satellite the data was recorded on, ELFIN -A or -B ..................... |               
    name::String                              # Single string that can be used to re-create an event.                          |
    duration::Float64                         # Event duration in seconds .................................................... | seconds
    n_datapoints::Int                         # Number of datapoints recorded during the event
    n_observations::Int                       # Number of distinct periods where instruments were recording in the event ..... |
    observation_edge_idxs::Vector{Int}        # Indices of where observations periods begin. Includes last index of data       |
    kp::Vector{Float64}                       # 3-hour Kp index for each datapoint ........................................... |
    dst::Vector{Float64}                      # Hourly Dst index for each datapoint                                            | nT

    # Spatial variable fields 
    time::Vector{Float64}                     # Time since event start for each datapoint                                      | seconds
    time_datetime::Vector{DateTime}           # Date and time of each datapoint as DateTime objects .......................... | 
    position::Dict                            # ELFIN's position at each time point in GEI coordinates. Access coordinates with
                                              # position["x"][t], position["y"][t], position["z"][t] ......................... | km
    L::Vector{Float64}                        # L value of ELFIN's location at each data point. Calculated with T87 model .... | 
    MLT::Vector{Float64}                      # MLT of ELFIN's location at each data point                                     | hour

    # Science fields
    pitch_angles::Matrix{Float64}             # Pitch angle of each of ELFIN's look directions at each timestep .............. | degrees
    energy_bins_min::Vector{Float64}          # Minimum value of each of ELFIN's energy channels ............................. | keV
    energy_bins_mean::Vector{Float64}         # Mean value of each of ELFIN's energy channels ................................ | keV
    energy_bins_max::Vector{Float64}          # Maximum value of each of ELFIN's energy channels ............................. | keV
    lc_idxs::Vector{Any}                      # Pitch angle indices inside the loss cone at each time step                     |
    alc_idxs::Vector{Any}                     # Pitch angle indices inside the anti loss cone at each time step                |
    loss_cone_angles::Vector{Float64}         # Angle of the loss cone at each time step ..................................... | degrees
    anti_loss_cone_angles::Vector{Float64}    # Angle of the anti-loss cone at each time step ................................ | degrees
    avg_pitch_angles::Vector{Float64}         # Pitch angle of each of ELFIN's look directions averaged over all timesteps ... | degrees
    avg_loss_cone_angle::Float64              # Average loss cone angle over full event ...................................... | degrees
    avg_anti_loss_cone_angle::Float64         # Average anti loss cone angle over full event ................................. | degrees

    e_flux::Array{Float64}                    # Electron energy flux at each timestep, energy channel, and look direction .... | kev/(cm^2 str s MeV)
                                                  # Dimensions: [time, energy channel, pitch angle]                            |
    n_flux::Array{Float64}                    # Electron number flux at each timestep, energy channel, and look direction .... | electrons/(cm^2 str s MeV)
                                                  # Dimensions: [time, energy channel, pitch angle]                            |
    
    Jprec_over_Jtrap::Array{Float64}          # Ratio of electron flux in the loss cone to flux in the trapped area at each time step and energy channel.  | 
                                                  # Dimensions: [time, energy channel]
end



########################################################################
#                       ~~~~~~~~~~~~~~~~~~~~~~~~                       #
#                       EVENT CREATION FUNCTIONS                       #
#                       ~~~~~~~~~~~~~~~~~~~~~~~~                       #
########################################################################
function create_event(name::String; warn = false)
# Creates an Event based on the name string of another event
    separator_idxs = findall('_', name)

    date = Date(name[1:separator_idxs[1]-1])
    time = Time(name[separator_idxs[1]+1:separator_idxs[2]-1])
    duration = Millisecond(parse(Float64, name[separator_idxs[2]+1:separator_idxs[3]-1]) * 1000)
    sat = name[end]

    start_datetime = DateTime(date, time)
    stop_datetime   = DateTime(date, time + duration)

    return create_event(start_datetime, stop_datetime, sat, warn = warn)
end

function create_event(date::Date, sat; warn = false)
# Creates an Event object from a date and satellite name. Event covers all data recorded on the given date.
    start_datetime = DateTime(date, Time("00:00:00"))
    stop_datetime   = DateTime(date, Time("23:59:59"))

    return create_event(start_datetime, stop_datetime, sat, warn = warn)
end

function create_event(start_datetime::DateTime, stop_datetime::DateTime, sat; warn = false)
# Creates an Event object from a start datetime, stop datetime, and satellite ID (either "a", "A", "b", or "B")

    #################### GUARD BLOCK ####################
    global warnings = warn # Sets whether we print non-fatal warnings to the terminal.
    year = string(Year(start_datetime).value, pad = 4)
    month = string(Month(start_datetime).value, pad = 2)
    day = string(Day(start_datetime).value, pad = 2)

    if lowercase(sat) ∉ ["a", "b", 'a', 'b']
        error("Satellite \"$(sat)\" not recognized.")
    end

    science_data_path = "$Julia_ELFIN_Tools_TOP_LEVEL/data/processed_scientific_data/$(year)$(month)$(day)_el$(sat).npz"
    position_data_path = "$Julia_ELFIN_Tools_TOP_LEVEL/data/processed_position_data/$(year)$(month)$(day)_el$(sat).npz"
    if !isfile(science_data_path)
        # Always warn as this is more critical
        @warn "\033[93mCould not find $(science_data_path).\033[0m"
        return nothing
    end
    if !isfile(position_data_path)
        # Always warn as this is more critical
        @warn "\033[93mCould not find $(position_data_path).\033[0m"
        return nothing
    end
    data = _load_data(science_data_path, position_data_path)

    # Check bounds
    time_from_start = data["fs_time"] .- start_datetime
    time_from_end   = data["fs_time"] .- stop_datetime
    # Warn if either the start or stop time is outside a recording period and set to nearest available time
    if min(abs.(time_from_start)...) >= Second(2) # If the nearest time is more than 2 seconds away. This means ~10 missed observations at ELFIN's sampling about every .2 seconds.
        _, idx = findmin(abs.(time_from_start))
        old_start = start_datetime
        start_datetime = data["fs_time"][idx] # Set start to nearest time we have data for.
        difference = time_from_start[idx].value / 1000
        if difference > 0
            difference = "+" * string(difference)
        end
        difference = string(difference)

        if warnings == true; @warn "\033[93mcreate_event(): No data exists for start time $(old_start). Starting event at nearest observation ($(difference) s).\033[0m"; end
    end
    if min(abs.(time_from_end)...) >= Second(2) # If the nearest time is more than 2 seconds away. This means ~10 missed observations at ELFIN's sampling about every .2 seconds.
        _, idx = findmin(abs.(time_from_end))
        old_end = stop_datetime
        stop_datetime = data["fs_time"][idx] # Set start to nearest time we have data for.
        difference = time_from_end[idx].value / 1000
        if difference > 0
            difference = "+" * string(difference)
        end
        difference = string(difference)

        if warnings == true; @warn "\033[93mcreate_event(): No data exists for end time $(old_end). Ending event at nearest observation ($(difference) s).\033[0m"; end
    end
    #################### END GUARD BLOCK ####################



    #################### BEGIN EVENT CREATION ####################
    indices_of_interest = (data["fs_time"] .>= start_datetime) .& (data["fs_time"] .<= stop_datetime)
    data["indices_of_interest"] = indices_of_interest # Save to data dict in case needed for debugging. Data dict will keep it hidden from end user.
    n_datapoints = sum(indices_of_interest)
    if n_datapoints == 0
        @warn "\033[93mNo data exists between start time and end time (Δt = $(stop_datetime - start_datetime)), no event created.\033[0m"
        return nothing
    end
        

    # Get science time grid
    time_datetime = data["fs_time"][indices_of_interest]
    duration = 3 + ((time_datetime[end] - time_datetime[1]).value) / 1000 # Units: seconds. Adding three because each datapoint represents the start of a ~three second recording interval
    time = time_datetime .- time_datetime[1] # Float time since start
    time  = [time[i].value for i = eachindex(time)] ./ 1000 # Convert to float. Units: seconds
    data["time"] = time

    # Get geomagnetic indices
    dst, kp = _get_geomagnetic_indices(time_datetime)


    # Get observation periods
    observation_edge_idxs = [1]
    append!(observation_edge_idxs, findall(diff(time_datetime) .>= Second(60)) .+ 1) # Plus one since diff() drops an index
    append!(observation_edge_idxs, [n_datapoints + 1]) # Plus one because obs periods only last until edge - 1 index to not overlap periods
    observation_edge_idxs = unique(observation_edge_idxs) # In case first index is the last before a long gap
    n_observations = length(observation_edge_idxs) - 1

    # Create event name
    name = Dates.format(start_datetime, dateformat"yyyy-mm-dd_HH:MM:SS") * "_" * string(duration) * "_EL$(data["satellite"])"


    # WARNINGS BLOCK ------
    # Ensure there is data available during the time specified
    if n_datapoints == 0
        date       = Dates.format(start_datetime, dateformat"yyyy-mm-dd")
        start_time = Dates.format(start_datetime, dateformat"HH:MM:SS")
        end_time   = Dates.format(stop_datetime, dateformat"HH:MM:SS")
        @warn "\033[93mcreate_event(): No data on $(date) between $(start_time) and $(end_time). Event empty.\033[0m"
        return nothing
    end

    ############################################################################
    #                    INTERPOLATION OF SPATIAL VARIABLES                    #
    ############################################################################
    # Get state time grid for interpolation
    state_datetime = data["state_time"]
    state_time_grid = cumsum([state_datetime[1] - time_datetime[1]; diff(state_datetime)])
    state_time_grid = [state_time_grid[i].value for i = eachindex(state_time_grid)] ./ 1000 # Convert to float. Units: seconds

    # Interpolate position from state grid to time grid
    position = Dict()
    position_state_grid = data["position"]
    axes = ["x", "y", "z"]
    for i = 1:3
        position_interpolation_object = LinearInterpolator(state_time_grid, position_state_grid[:,i], NoBoundaries())
        position[axes[i]] = position_interpolation_object.(time)
    end

    # Interpolate L-shell to science time grid
    L_state_grid = data["L"]
    L_interpolation_object = LinearInterpolator(state_time_grid, L_state_grid, NoBoundaries())
    L = L_interpolation_object.(time)

    # Interpolate MLT
    MLT_state_grid = data["MLT"]
    # MLT has the unique issue that it wraps at 24, causing erroneous interpolation results in the period between wraps. So we check for those and treat them seperately.
    midnight_crossing_idxs = findall(abs.(diff(MLT_state_grid)) .> 23.5) # This index is on the state grid
    continuous_MLT_edge_idxs = cat(1, midnight_crossing_idxs .+ 1, length(MLT_state_grid) + 1, dims = 1)

    # Now, for each period of continuous MLTs, we interpolate the MLT and then add it to the final MLT result on the science grid
    MLT = Float64[]
    n_continuous_MLT_periods = length(continuous_MLT_edge_idxs)-1
    for i = 1:n_continuous_MLT_periods
        slice = continuous_MLT_edge_idxs[i]:continuous_MLT_edge_idxs[i+1]-1
        
        # Get the datapoitns that we are going to interpolate to (ie science grid)
        interpolation_start_time = state_time_grid[slice[1]]
        if i == n_continuous_MLT_periods # Last period only
            interpolation_stop_time  = 0 # Doesn't matter for last period
        else
            # Normal behaviour
            interpolation_stop_time  = state_time_grid[slice[end]+1]
        end

        # Create the interpolant from the raw data
        interpolation_time = state_time_grid[slice]
        interpolation_MLT = MLT_state_grid[slice]

        # If the slice is only 1 index long (ie wrap happens on the second or the second to last index), we add a ficticious point to the interpolant so it can interpolate.
        # This should be okay to do, as the most it will be off by is ~.1 hr MLT.
        if length(slice) == 1
            append!(interpolation_time, interpolation_time[end] + mean(diff(state_time_grid))) # Because the state grid is uniform, the mean(diff()) call should get the state grid timestep
            append!(interpolation_MLT, round(interpolation_MLT[1])) # round() will get either 0 or 24 depending on which direction we are crossing midnight in
        end

        MLT_interpolation_object = LinearInterpolator(interpolation_time, interpolation_MLT, NoBoundaries())


        # Normal behaviour
        science_grid_continuous_MLT_slice = interpolation_start_time .<= time .< interpolation_stop_time
        # First index only in case there's science data before the first state data
        if i == 1
            science_grid_continuous_MLT_slice = time .< interpolation_stop_time
        end
        # Last index only in case there's science data after the last state data
        if i == n_continuous_MLT_periods
            science_grid_continuous_MLT_slice = interpolation_start_time .<= time
        end
        # If there's only one period for the whole event
        if n_continuous_MLT_periods == 1
            science_grid_continuous_MLT_slice = 1:length(time)
        end

        # Do the interpolation and add it to vector
        append!(MLT, MLT_interpolation_object.(time[science_grid_continuous_MLT_slice]))
    end

    # Assert they're the same length for sanity
    @assert(length(MLT) == length(time),
        "MLT interpolation error. length(MLT) = $(length(MLT)), length(time) = $(length(time)).")
    @assert(length(L) == length(time),
        "L interpolation error. length(L) = $(length(MLT)), length(time) = $(length(time)).")

    ##########################################################################
    # Bin data into ELFIN's discrete look directions
    binned = _bin_data(data, n_datapoints, indices_of_interest)
    if binned == nothing
        # Warning contained in _bin_data()
        return nothing
    end

    # Calculate loss/anti-loss cones over the observation
    _calculate_loss_cones(data, binned, n_datapoints, indices_of_interest) # Mutates data dict
    _calculate_Jprec_over_Jtrap(data, binned)

    # Get average angles
    avg_pitch_angles = dropdims(mean(binned["pitch_angles"], dims = 1), dims = 1) # Average pitch angle each bin is facing over time
    avg_loss_cone_angle = mean(data["loss_cone_angles"]) # Average loss cone angle over the observation period
    avg_anti_loss_cone_angle = mean(data["anti_loss_cone_angles"]) # Average anti loss cone angle over the observation period


    #############################################################################
    #                           FINISH EVENT CREATION                           #
    #############################################################################
    # Create the Event object and return it
    return Event(
        data["satellite"], name, duration, n_datapoints, n_observations, observation_edge_idxs, kp, dst,

        time, time_datetime, position, L, MLT, 
        
        binned["pitch_angles"], data["energy_bins_min"], data["energy_bins_mean"], data["energy_bins_max"], data["loss_cone_idxs"], data["anti_loss_cone_idxs"], data["loss_cone_angles"], data["anti_loss_cone_angles"], avg_pitch_angles, avg_loss_cone_angle, avg_anti_loss_cone_angle,
                 
        binned["e_flux"], binned["n_flux"], data["Jprec_over_Jtrap"]
    )
end

function _load_data(science_data_path, position_data_path)
# Load ELFIN data from cleaned .npz into a dictionary
    np = pyimport("numpy")
    science_data = np.load(science_data_path,   allow_pickle=true)
    position_data = np.load(position_data_path, allow_pickle=true)

    # Determine if data is from ELFIN-A or ELFIN-B
    # Doing this because Julia can't read keys that are a single string >:(
    ela = get(science_data, :ela)
    if ela[1] == true
        satellite = "A"
    else
        satellite = "B"
    end

    return Dict("satellite"        => satellite,
                "et_time"          => _timestamp_to_DateTime(get(science_data, :et_time)),
                "hs_time"          => _timestamp_to_DateTime(get(science_data, :hs_time)),
                "fs_time"          => _timestamp_to_DateTime(get(science_data, :fs_time)),
                "Et_nflux"         => get(science_data, :Et_nflux),
                "Et_eflux"         => get(science_data, :Et_eflux),
                "Et_dfovf"         => get(science_data, :Et_dfovf),
                "energy_bins_min"  => get(science_data, :energy_bins_min),
                "energy_bins_mean" => get(science_data, :energy_bins_mean),
                "energy_bins_max"  => get(science_data, :energy_bins_max),
                "pa"               => get(science_data, :pa),
                "spinphase"        => get(science_data, :spinphase),
                "sectnum"          => get(science_data, :sectnum),
                "Tspin"            => get(science_data, :Tspin),
                "hs_Epat_nflux"    => get(science_data, :hs_Epat_nflux),
                "hs_Epat_eflux"    => get(science_data, :hs_Epat_eflux),
                "hs_Epat_dfovf"    => get(science_data, :hs_Epat_dfovf),
                "hs_LCdeg"         => get(science_data, :hs_LCdeg),
                "hs_antiLCdeg"     => get(science_data, :hs_antiLCdeg),
                "hs_epa_spec"      => get(science_data, :hs_epa_spec),
                "fs_Epat_nflux"    => get(science_data, :fs_Epat_nflux),
                "fs_Epat_eflux"    => get(science_data, :fs_Epat_eflux),
                "fs_Epat_dfovf"    => get(science_data, :fs_Epat_dfovf),
                "fs_LCdeg"         => get(science_data, :fs_LCdeg),
                "fs_antiLCdeg"     => get(science_data, :fs_antiLCdeg),
                "fs_epa_spec"      => get(science_data, :fs_epa_spec),
                "nspinsinsum"      => get(science_data, :nspinsinsum),
                "nsectors"         => get(science_data, :nsectors),
                "sect2add"         => get(science_data, :sect2add),
                "spinph2add"       => get(science_data, :spinph2add),
                "state_time"       => _timestamp_to_DateTime(get(position_data, :state_time)),
                "L"                => abs.(vec(get(position_data, :L))),
                "MLT"              => get(position_data, :MLT),
                "altitude"         => get(position_data, :altitude),
                "position"         => get(position_data, :position) .* 6378 # Convert Re back to km
                )
end

function _timestamp_to_DateTime(times)
# Convert timestamp string to DateTime object
    date_time = String.(collect(times))
    date_time = DateTime.(date_time, dateformat"yyyy.mm.dd.HH.MM.SS.sss")
end

function _get_geomagnetic_indices(time_datetime)
    results = Dict()
    for index_name in ["dst", "kp"]
        path = "$Julia_ELFIN_Tools_TOP_LEVEL/data/$(index_name).csv"
        index_data = readdlm(path, ',', skipstart = 1)
        index_times = DateTime.(index_data[:,1], dateformat"yyyy-mm-dd/HH:MM:SS")
        geomag_index = float.(index_data[:,3])

        results[index_name] = fill(0.0, length(time_datetime))
        for t = 1:length(time_datetime)
            _, nearest_time_idx = findmin(abs.(index_times .- time_datetime[t]))
            results[index_name][t] = geomag_index[nearest_time_idx]
        end
    end
    return results["dst"], results["kp"]
end

function _bin_data(data, n_datapoints, indices_of_interest)
# Bin pitch angle and flux data into ELFIN's discrete look directions
    total_observations = length(data["fs_Epat_nflux"][:,1,1]) # Dimensions of data["fs_Epat_nflux"]: [time, pitch angle, energy channel]

    nflux        = zeros(total_observations, 16, 16) # Dimensions: [time, energy channel, pitch angle] Units: electrons / (cm^2 str s)
    eflux        = zeros(total_observations, 16, 16) # Dimensions: [time, energy channel, pitch angle] Units: keV / (cm^2 str s)
    pitch_angles = zeros(total_observations, 16)     # Dimensions: [time, pitch angle] Units: 

    for i = 1:total_observations
        nflux[i,:,:] = (data["fs_Epat_nflux"][i, 2:end-1, :])' # Trim NaNs that pad the pitch angle dimension on the first & last columns. Transpose as well. Dimensions: [energy channel, pitch angle] This makes plotting easier: cols = pitch angle & rows = energy channel
        eflux[i,:,:] = (data["fs_Epat_eflux"][i, 2:end-1, :])'    
        pitch_angles[i,:] = data["fs_epa_spec"][i,2:end-1] # Trim corresponding pitch angles
    end


    if sum(isnan.(nflux)) ≠ 0
        if warnings == true; @warn "\033[93m_bin_data(): n_flux data contains NaN.\033[0m"; end
        nflux[isnan.(nflux)] .= 0
    end
    if sum(isnan.(eflux)) ≠ 0
        if warnings == true; @warn "\033[93m_bin_data(): e_flux data contains NaN.\033[0m"; end
        eflux[isnan.(eflux)] .= 0
    end
    if sum(isnan.(pitch_angles)) ≠ 0
        if warnings == true; @warn "\033[93m_bin_data(): Pitch angle data contains NaN.\033[0m"; end
        return nothing
    end


    # Trim cleaned data to only the time of interest
    nflux        = nflux[indices_of_interest, :, :]
    eflux        = eflux[indices_of_interest, :, :]
    pitch_angles = pitch_angles[indices_of_interest, :]

    # Collect data into one structure for easy access elsewhere
    binned = Dict("pitch_angles" => pitch_angles,
                  "n_flux" => nflux,
                  "e_flux" => eflux)
    return binned
end

function _calculate_loss_cones(data, binned, n_datapoints, indices_of_interest)
# Calculate loss cones for an observation period
    # Get loss cone angles - loss cone angle is only given at each half rotation, so we have to propagate it out to all data points
    loss_cone_angles = data["fs_LCdeg"][indices_of_interest]
    anti_loss_cone_angles = data["fs_antiLCdeg"][indices_of_interest]

    # Get indices of loss and anti-loss cones. The direction of LC/ALC depends on hemisphere.
    loss_cone_idxs      = []
    anti_loss_cone_idxs = []

    for i = 1:n_datapoints
        lc_smaller_than_alc = loss_cone_angles[i] < anti_loss_cone_angles[i]
        if lc_smaller_than_alc # Northern hemisphere, loss cone line is on left side of flux chart
            # We want smaller angles for loss cone, larger angles for anti-loss cone
            append!(loss_cone_idxs,      [findall(binned["pitch_angles"][i,:] .<= loss_cone_angles[i])])
            append!(anti_loss_cone_idxs, [findall(binned["pitch_angles"][i,:] .>= anti_loss_cone_angles[i])])        
        else # Southern hemisphere, loss cone line is on right side of flux chart
            # We want larger angles for the loss cone and smaller angles for anti-loss
            append!(loss_cone_idxs,      [findall(binned["pitch_angles"][i,:] .>= loss_cone_angles[i])])
            append!(anti_loss_cone_idxs, [findall(binned["pitch_angles"][i,:] .<= anti_loss_cone_angles[i])])        
        end
    end

    data["loss_cone_angles"] = loss_cone_angles
    data["loss_cone_idxs"] = loss_cone_idxs
    data["anti_loss_cone_angles"] = anti_loss_cone_angles
    data["anti_loss_cone_idxs"] = anti_loss_cone_idxs
end

function _calculate_Jprec_over_Jtrap(data, binned)
    time = data["time"]
    pitch_angles = binned["pitch_angles"]

    # Get fluxes. Doesn't matter if we use number flux or energy flux since they're related by a scalar, which disappears when we take the ratio. Either one gives the exact same result.
    n_flux = binned["n_flux"]
    loss_cone_n_flux = zeros(length(time), 16)
    trapped_n_flux   = zeros(length(time), 16)

    for t = 1:length(time)
        # Get loss cone and trapped indexes, angles
        lc_idxs = data["loss_cone_idxs"][t]
        lc_pitch_angles = pitch_angles[t,:][lc_idxs]

        trapped_idxs = setdiff(1:16, data["loss_cone_idxs"][t], data["anti_loss_cone_idxs"][t])
        trapped_pitch_angles = pitch_angles[t,:][trapped_idxs]

        # Calculate coverage of the loss cone
        lc_angle = data["loss_cone_angles"][t]
        if lc_angle < 90 # Northern hemisphere
            lc_size_deg = lc_angle
            max_distance_into_lc = lc_angle - min(pitch_angles[t,:]...)
        else # Southern hemisphere, lc_angle > 90
            lc_size_deg = 180 - lc_angle
            max_distance_into_lc = max(pitch_angles[t,:]...) - lc_angle
        end

        # NaN out if there are no loss cone/trapped observations, OR we have less than 15º of the loss cone (since less than that is not representative)
        if (max_distance_into_lc <= 15) || (lc_idxs == []) || (trapped_idxs == [])
            loss_cone_n_flux[t,:] .= NaN
            continue
        end

        # Average the flux in the loss cone and trapped region for each energy bin at this time. Do this for both energy and number flux.
        loss_cone_n_flux[t,:] = [mean(n_flux[t, E, lc_idxs])      for E = 1:16]
        trapped_n_flux[t,:]   = [mean(n_flux[t, E, trapped_idxs]) for E = 1:16]
        
        
    end

    # Noise removal - removing places that have low trapped counts that would artificially spike the ratio
    trapped_n_flux[trapped_n_flux .<= 10^2.5] .= 0

    # Get ratio. NaN => prec and trapped both zero, Inf => trapped is zero
    Jprec_over_Jtrap = loss_cone_n_flux ./ trapped_n_flux

    # Mutate data dict to return results
    data["Jprec_over_Jtrap"] = Jprec_over_Jtrap
end


########################################################################
#                       ~~~~~~~~~~~~~~~~~~~~~~~~                       #
#                            EVENT METHODS                             #
#                       ~~~~~~~~~~~~~~~~~~~~~~~~                       #
########################################################################
function integrate_flux(event::Event; time = false, energy = false, pitch_angle = false,
                        time_range = 1:event.n_datapoints,
                        pitch_angle_range = "full",
                        energy_range_keV = (-Inf, Inf)
                        )
    if pitch_angle_range == "full"
        pitch_angle_idxs = [1:16 for i = 1:event.n_datapoints]
    elseif pitch_angle_range == "loss cone"
        pitch_angle_idxs = event.lc_idxs
    elseif pitch_angle_range == "anti loss cone"
        pitch_angle_idxs = event.alc_idxs
    elseif pitch_angle_range == "trapped"
        pitch_angle_idxs = [setdiff(1:16, event.lc_idxs[t], event.alc_idxs[t]) for t = 1:event.n_datapoints]
    else
        error("Keyword argument pitch_angle_range=$(pitch_angle_range) not recognized. Please use \"full\", \"loss cone\", or \"anti loss cone\".")
    end

    energy_range_idx = energy_range_keV[1] .<= event.energy_bins_mean .<= energy_range_keV[2]

    e_flux = copy(event.e_flux)
    n_flux = copy(event.n_flux)

    # Matrix slice indices. We use the full dimension if it hasn't been integrated. Slice is changed to [1] if it has in order to flatten result.
    pa = 1:16
    e = 1:16
    t = 1:event.n_datapoints

    # Do integrations in reverse order of their dimension (ie pitch angle, energy, then time). This prevents one integration from changing another.
    if pitch_angle == true
       e_flux, n_flux = _integrate_over_pitch_angle(event, e_flux, n_flux, pitch_angle_idxs)
       pa = 1
    end

    if energy == true
        e_flux, n_flux = _integrate_over_energy(event, e_flux, n_flux, energy_range_idx)
        e = 1
    end

    if time == true
        e_flux, n_flux = _integrate_over_time(event, e_flux, n_flux, time_range)
        t = 1
    end

    return e_flux[t, e, pa], n_flux[t, e, pa]
end

function _integrate_over_pitch_angle(event::Event, e_flux, n_flux, pitch_angle_idxs)
    for t = 1:event.n_datapoints
        # Don't continue if we don't have loss cone coverage
        if length(pitch_angle_idxs[t]) < 2
            e_flux[t, :, :] .= NaN
            continue
        end

        α = event.pitch_angles[t, pitch_angle_idxs[t]]
        for E = 1:16
            e_flux[t, E, :] .= 2π * integrate(α, e_flux[t, E, pitch_angle_idxs[t]] .* sind.(α), Trapezoidal())
            n_flux[t, E, :] .= 2π * integrate(α, n_flux[t, E, pitch_angle_idxs[t]] .* sind.(α), Trapezoidal())
        end
    end
    return e_flux, n_flux
end

function _integrate_over_energy(event::Event, e_flux, n_flux, energy_range)
    for t = 1:event.n_datapoints
        for α = 1:16
            e_flux[t, :, α] .= integrate(event.energy_bins_mean[energy_range] ./ 1000, e_flux[t, energy_range, α], Trapezoidal())
            n_flux[t, :, α] .= integrate(event.energy_bins_mean[energy_range] ./ 1000, n_flux[t, energy_range, α], Trapezoidal())
        end
    end
    return e_flux, n_flux
end

function _integrate_over_time(event::Event, e_flux, n_flux, time_range)
    # Different behavior for 1-datapoint events or input where time range is only one index
    if length(time_range) == 1
        for t = 1:length(e_flux[:, 1, 1])
            e_flux[t, :, :] = event.e_flux[time_range, :, :] .* 3 # T_spin ~= 3 seconds
            n_flux[t, :, :] = event.n_flux[time_range, :, :] .* 3
        end
        return e_flux, n_flux
    end

    for E = 1:16
        for α = 1:16
            e_flux[:, E, α] .= integrate(event.time[time_range], e_flux[time_range, E, α], Trapezoidal())
            n_flux[:, E, α] .= integrate(event.time[time_range], n_flux[time_range, E, α], Trapezoidal())
        end
    end
    return e_flux, n_flux
end