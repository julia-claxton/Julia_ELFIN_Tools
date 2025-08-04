#=
I was originally saving time as a string and using PyCall in Julia to access them, as NPZ doesn't
have support for strings. However, PyCall is not thread-safe and my analysis takes way too long
without multithreading. This converts the original files produced by the python ELFIN downloader
to files safe for the Julia NPZ library to read. It's hacky, but it works.
=#

using Statistics, LinearAlgebra          # Core math
using BenchmarkTools, Profile, TickTock  # Debugging
using NPZ, DelimitedFiles                # File interactions
using Dates
using PyCall
    np = pyimport("numpy")

function split_and_save(data_dict, prefix, datetimes)
    h = [Hour(t).value for t in datetimes]
    m = [Minute(t).value for t in datetimes]
    s = [Second(t).value for t in datetimes]
    ms = [Millisecond(t).value for t in datetimes]

    data_dict["$(prefix)_hour"] = h
    data_dict["$(prefix)_minute"] = m
    data_dict["$(prefix)_second"] = s
    data_dict["$(prefix)_millisecond"] = ms
end


function main_science()
    dir = "/Users/luna/Research/Julia_ELFIN_Tools/data/_backup_processed_scientific_data"
    files = readdir(dir)
    for filename in files
        if filename[end-3:end] ≠ ".npz"; continue; end # Skip non-data files
        original_path = "$(dir)/$(filename)"

        updated_data_dict = Dict{String, Any}()

        year = parse(Int, original_path[end-15:end-12])
        month = parse(Int, original_path[end-11:end-10])
        day = parse(Int, original_path[end-9:end-8])

        updated_data_dict["year"] = year
        updated_data_dict["month"] = month
        updated_data_dict["day"] = day

        data = np.load(original_path, allow_pickle = true)
        data_keys = data[:files]

        for key in data_keys
            if key in ["et_time", "hs_time", "fs_time"]; continue; end

            to_write = collect(data[:__getitem__](key))
            if length(to_write) == 0; to_write = NaN; end # NPZ fails reading empty arrays
            updated_data_dict[key] = to_write
        end

        # Save datetime in new format
        et_datetime = DateTime.(get(data, :et_time), dateformat"yyyy.mm.dd.HH.MM.SS.sss")
        split_and_save(updated_data_dict, "et", et_datetime)

        hs_datetime = DateTime.(get(data, :hs_time), dateformat"yyyy.mm.dd.HH.MM.SS.sss")
        split_and_save(updated_data_dict, "hs", hs_datetime)
        
        fs_datetime = DateTime.(get(data, :fs_time), dateformat"yyyy.mm.dd.HH.MM.SS.sss")
        split_and_save(updated_data_dict, "fs", fs_datetime)

        # Save new file
        destination_path = replace(original_path, "_backup_" => "")
        npzwrite(destination_path, updated_data_dict)
    end
end

function main_position()
    dir = "/Users/luna/Research/Julia_ELFIN_Tools/data/_backup_processed_position_data"
    files = readdir(dir)
    for filename in files
        if filename[end-3:end] ≠ ".npz"; continue; end # Skip non-data files
        original_path = "$(dir)/$(filename)"

        updated_data_dict = Dict{String, Any}()

        year = parse(Int, original_path[end-15:end-12])
        month = parse(Int, original_path[end-11:end-10])
        day = parse(Int, original_path[end-9:end-8])

        updated_data_dict["year"] = year
        updated_data_dict["month"] = month
        updated_data_dict["day"] = day

        data = np.load(original_path, allow_pickle = true)
        data_keys = data[:files]

        for key in data_keys
            if key in ["state_time"]; continue; end

            to_write = collect(data[:__getitem__](key))
            if length(to_write) == 0; to_write = NaN; end # NPZ fails reading empty arrays
            updated_data_dict[key] = to_write
        end

        # Save datetime in new format
        state_datetime = DateTime.(get(data, :state_time), dateformat"yyyy.mm.dd.HH.MM.SS.sss")
        split_and_save(updated_data_dict, "state_time", state_datetime)

        # Save new file
        destination_path = replace(original_path, "_backup_" => "")
        npzwrite(destination_path, updated_data_dict)
    end
end

main_science()
main_position()