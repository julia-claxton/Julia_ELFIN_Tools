# Housekeeping
const Julia_ELFIN_Tools_TOP_LEVEL = @__DIR__
@assert endswith(Julia_ELFIN_Tools_TOP_LEVEL, "Julia_ELFIN_Tools")

# Library includes
using Statistics
using LinearAlgebra
using Dates
using Plots
using Plots.PlotMeasures
using ColorSchemes

include("$(Julia_ELFIN_Tools_TOP_LEVEL)/Events.jl")

USE_VSCODE_PLOTS = false # This variable should be set to true if you are running Julia from VSCode. This tailors the plots to be a bit better when using VSCode.

"""
Visualization.jl

This file provides several methods for quickly visualizing the data stored in an Event object. The
most important method is quicklook(), which takes in an Event and displays a plot containing an
overview of the event.

Written by Julia Claxton (julia.claxton@colorado.edu)
Released under MIT License (see LICENSE.txt for full license)
"""

function quicklook(event::Event; by = nothing)
# `by` kwarg determines what the x-axis of the plots should be and can take values of "index", "time", or "date".
    if event.n_datapoints == 1; @warn "n_datapoints = 1, no plot"; return nothing; end
    
    # Auto-set x-axis if not user specified
    if (by == nothing) && (event.n_observations == 1)
        by = "date"
    elseif (by == nothing) && (event.n_observations > 1)
        by = "index"
    end

    # Create plots
    energy_heatmap = energy_time_series(          event, by = by, show_plot = false)
    lc_ratio       = Jprec_over_Jtrap_time_series(event, by = by, show_plot = false)
    pad_heatmap    = pad_time_series(             event, by = by, show_plot = false)
    L_shell        = L_MLT_time_series(           event, by = by, show_plot = false)
    abs_sunset     = absolute_sunset(             event, show_plot = false)
    rel_sunset     = relative_sunset(             event, show_plot = false)
    mlt_l_path     = MLT_L_track(                 event, show_plot = false)

    blank = plot(yticks = nothing, framestyle = :none)

    l = @layout [ plot1{.14h}
                  plot2{.14h}
                  plot3{.14h}
                  blank plot4{.93w} blank{.1125w}
                  plot5 plot6
    ]

    plot(lc_ratio, energy_heatmap, pad_heatmap, blank, L_shell, blank, abs_sunset, rel_sunset,
        layout = l,
        plot_title = "$(Date(event.time_datetime[1])) $(Time(event.time_datetime[1])) ELFIN-$(event.satellite)\nDuration = $(event.duration) s / $(event.n_datapoints) observations",
        background_color = :white,
        size = (1,1.4) .* 900,
        format = :png,
        dpi = 300
    )
    _display_plot()
end

function quicklook(event::Nothing; by = nothing)
    @warn "No event"
    return nothing
end

function plot_event(event::Event; colormap = :ice)
    xtick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 5)))
    heatmap(1:event.n_datapoints, event.energy_bins_mean, log10.(event.Jprec_over_Jtrap'),
        xlabel = "$(Date(event.time_datetime[1])) ELFIN-$(event.satellite)",
        xticks = (xtick_idxs, Dates.format.(Time.(event.time_datetime[xtick_idxs]), "HH:MM:SS")),
        xminorticks = true,

        ylabel = "Energy (keV)",
        ylims = (event.energy_bins_min[begin], event.energy_bins_max[end]),
        yscale = :log10,
        yminorticks = true,

        colorbar_title = "\nLog10 Jprec/Jtrap",
        colormap = colormap,
        clims = (-1.25, .25),

        leftmargin = 5mm,
        rightmargin = 5mm,
        
        background_color_inside = RGB(.8, .8, .8),
        background_color_outside = :transparent,
        framestyle = :box,
        tickdirection = :out,
        aspect_ratio = (event.n_datapoints/(event.energy_bins_mean[end]-event.energy_bins_mean[begin])) * .2,
        size = (2, .5) .* 350,
        dpi = 300
    )
    return plot!()
end

function _strip_heatmap(z; colormap = cgrad(:inferno))
    dims = size(z)
    x = 1:dims[2]
    y = 1:dims[1]

    return _strip_heatmap(x, y, z, colormap = colormap)
end

function _strip_heatmap(x, y, z; colormap = cgrad(:inferno), bg = :black)
    n_x = length(x)
    x_min = min(x...)
    x_max = max(x...)
    x_range = x_max - x_min

    n_y = length(y)
    y_min = min(y...)
    y_max = max(y...)
    y_range = y_max - y_min

    heatmap(x, y, z,
        background_color_inside = bg,
        colormap = colormap,
        framestyle = :grid,
        grid = true,
        tickdirection = :out,
        xlims = (x_min, x_max),
        ylims = (y_min, y_max),
        aspect_ratio = (x_range / y_range) * .15,
        dpi = 300
    )
    return plot!()
end

function energy_time_series(event::Event; by = "index", show_plot = true)
    if by == "index"
        x = eachindex(event.time)
        x_label = "Data Index"
        obs_edges = event.observation_start_idxs
    elseif by == "time"
        x = event.time
        x_label = "Time (s)"
        obs_edges = event.time[event.observation_start_idxs]
    elseif by == "date"
        x = event.time
        x_label = "Time (UTC), $(Date(event.time_datetime[1]))"
        obs_edges = event.time[event.observation_start_idxs]
    else
        error("Keyword argument by=\"$(by)\" not recognized.")
    end

    energy = event.energy_bins_mean ./ 1000 # MeV
    e_flux, _ = integrate_flux(event, pitch_angle = true)
    _strip_heatmap(x, log10.(energy), log10.(e_flux'), colormap = :ice)

    plot!(
        title = "Omnidirectional Flux",
        xlabel = x_label,
        ylabel = "Energy (MeV)",
        colorbar_title = "Log10 Energy Flux\n(keV/(cm² MeV s))",
        clims = (4, 8.5),
        yticks = (log10.(energy[1:3:16]), round.(energy[1:3:16], sigdigits = 2)),
        margin = 5mm,
        size = (1.3, .33) .* 500
    )
    # Add date ticks if needed
    if by == "date"
        tick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 6)))
        ticks = Dates.format.(event.time_datetime[tick_idxs], "HH:MM:SS")
        plot!(xticks = (event.time[tick_idxs], ticks))
    end
    vline!(obs_edges,
        linewidth = 1,
        linecolor = :white,
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display(plot!())
    end
    return plot!()
end

function Jprec_over_Jtrap_time_series(event::Event; by = "index", show_plot = true)
    by = lowercase(by)
    if by == "index"
        x = eachindex(event.time)
        x_label = "Data Index"
        obs_edges = event.observation_start_idxs
    elseif by == "time"
        x = event.time
        x_label = "Time (s)"
        obs_edges = event.time[event.observation_start_idxs]
    elseif by == "date"
        x = event.time
        x_label = "Time (UTC), $(Date(event.time_datetime[1]))"
        obs_edges = event.time[event.observation_start_idxs]
    else
        error("Keyword argument by=\"$(by)\" not recognized.")
    end    

    energy = event.energy_bins_mean ./ 1000 # MeV
    _strip_heatmap(x, log10.(energy), log10.(event.Jprec_over_Jtrap'), bg = RGB(.8,.8,.8), colormap = :ice)

    plot!(
        title = "Jprec/Jtrap",
        xlabel = x_label,
        ylabel = "Energy (MeV)",
        colorbar_title = "Log10 Jprec/Jtrap",
        background = :transparent,
        grid = false,
        clims = (-1.25, .25),
        yticks = (log10.(energy[1:3:16]), round.(energy[1:3:16], sigdigits = 2)),
        margin = 5mm,
        size = (1.3, .33) .* 500,
        dpi = 300
    )
    # Add date ticks if needed
    if by == "date"
        tick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 6)))
        ticks = Dates.format.(event.time_datetime[tick_idxs], "HH:MM:SS")
        plot!(xticks = (event.time[tick_idxs], ticks))
    end
    vline!(obs_edges,
        linewidth = 1,
        linecolor = :white,
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display("image/png", plot!())
    end
    return plot!()
end

function pad_time_series(event::Event; by = "index", show_plot = true)
    if by == "index"
        time = collect(eachindex(event.time))
        x_label = "Data Index"
        obs_edges = event.observation_start_idxs
    elseif by == "time"
        time = event.time
        x_label = "Time (s)"
        obs_edges = event.time[event.observation_start_idxs]
    elseif by == "date"
        time = event.time
        x_label = "Time (UTC), $(Date(event.time_datetime[1]))"
        obs_edges = event.time[event.observation_start_idxs]
    else
        error("Keyword argument by=\"$(by)\" not recognized.")
    end
    
    # Get pitch angle bin edges
    distance_between_bins = diff(event.pitch_angles, dims = 2)
    pitch_angle_bin_edges = event.pitch_angles[:, 1:end-1] .+ (distance_between_bins ./ 2)
    # Add first and last edge
    pitch_angle_bin_edges = cat(event.pitch_angles[:,1] .- (distance_between_bins[:,1] ./ 2), pitch_angle_bin_edges, dims = 2)
    pitch_angle_bin_edges = cat(pitch_angle_bin_edges, event.pitch_angles[:,end] .+ (distance_between_bins[:,end] ./ 2), dims = 2)

    # Get fluxes
    e_flux, n_flux = integrate_flux(event, energy = true)
    flux = log10.(e_flux) # Flux to plot

    # Create grid for the image
    deltas = diff(time)
    append!(time, time[end] + deltas[end]) # Add another index so that the time grid represents edges
    image_grid_time = LinRange(0, time[end], length(time) * 2)
    image_grid_pitch_angle = LinRange(0, 180, 180 * 2)
    image = fill(NaN, length(image_grid_time), length(image_grid_pitch_angle))

    # Fill the image with the flux data
    for t = 1:length(time)-1
        for α = 1:16
            image_t_slice = time[t] .<= image_grid_time .< time[t+1]
            image_pa_slice = pitch_angle_bin_edges[t, α] .<= image_grid_pitch_angle .< pitch_angle_bin_edges[t, α+1]

            image[image_t_slice, image_pa_slice] .= flux[t, α]
        end
    end

    # Plot image
    replace!(image, -Inf => -100) # Makes -Inf plot as black instead of not plotting at all. Lets us differentiate between no flux vs. no measurements
    plot(
        title = "Pitch Angle Distribution",
        xlabel = x_label,
        ylabel = "Pitch Angle (deg)",
        colorbar_title = "Log10 Energy Flux\nkev/(s cm² str)",

        xlims = (time[1], time[end-1]),
        ylims = (0, 180),

        framestyle = :box,
        background_color_inside = RGB(.8,.8,.8),
        background_color_outside = :transparent,
        grid = false,
        markershape = :square,
        markersize = (distance_between_bins ./ 5) .+ 1,
        markerstrokewidth = 0,
        tickdirection = :out,

        aspect_ratio = ((time[end-1] - time[1]) / 180) * .15,
        leftmargin = 8mm,
        size = (1, .3) .* 600,
        dpi = 300
    )
    heatmap!(image_grid_time, image_grid_pitch_angle, image',
        colormap = :ice,
        clims = (4.0, 8.5)
    )
    # Add date ticks if needed
    if by == "date"
        tick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 6)))
        ticks = Dates.format.(event.time_datetime[tick_idxs], "HH:MM:SS")
        plot!(xticks = (event.time[tick_idxs], ticks))
    end
    deleteat!(time, length(time)) # Remove extra time index we added before for the heatmap plot to work

    # Plot loss cones
    plot!(time, event.loss_cone_angles,
        linewidth = 2,
        linecolor = :white,
        linestyle = :solid,
        label = ""
    )
    plot!(time, event.anti_loss_cone_angles,
        linewidth = 2,
        linecolor = :white,
        linestyle = :dash,
        label = ""
    )
    vline!(obs_edges,
        linewidth = 1,
        linecolor = :white,
        linestyle = :dash,
        label = ""
    )

    if show_plot == true
        display("image/png", plot!())
    end
    return plot!()
end

function L_MLT_time_series(event::Event; by = "index", show_plot = true)
    if by == "index"
        x = eachindex(event.time)
        x_label = "Data Index"
        obs_edges = event.observation_start_idxs
    elseif by == "time"
        x = event.time
        x_label = "Time (s)"
        obs_edges = event.time[event.observation_start_idxs]
    elseif by == "date"
        x = event.time
        x_label = "Time (UTC), $(Date(event.time_datetime[1]))"
        obs_edges = event.time[event.observation_start_idxs]
    else
        error("Keyword argument by=\"$(by)\" not recognized.")
    end

    n_x = length(x)
    x_min = min(x...)
    x_max = max(x...)
    x_range = x_max - x_min

    y = event.L
    y_range = 15 - 0

    plot(
        xlabel = x_label,
        xlims = (x_min, x_max),

        ylabel = "L (T87)",
        ylims = (0, 15),
        yticks = 0:5:15,
        minorgrid = true,
        yminorticks = 5, 

        aspect_ratio = (x_range / y_range) * .15,
        framestyle = :box,
        leftmargin = -5.2mm,
        tickdirection = :out
    )
    # Add date ticks if needed
    if by == "date"
        tick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 6)))
        ticks = Dates.format.(event.time_datetime[tick_idxs], "HH:MM:SS")
        plot!(xticks = (event.time[tick_idxs], ticks))
    end

    # Plot each observation period
    for i = 1:event.n_observations
        slice = event.observation_start_idxs[i]:(event.observation_stop_idxs[i])
        plot!(x[slice], y[slice],
            label = "",
            linewidth = 1.5,
            linecolor = :black
        )
    end
    for i = 1:event.n_observations
        slice = event.observation_start_idxs[i]:(event.observation_stop_idxs[i])
        plot!(twinx(), x[slice], event.MLT[slice],
            xlims = (x_min, x_max),
            ylims = (0, 24),
            ylabel = "MLT",
            yguidefontcolor = :grey,
            y_foreground_color_axis = :grey,
            y_foreground_color_text = :grey,
            aspect_ratio = (x_range / 24) * .15,
            label = "",
            linewidth = 1.5,
            linecolor = RGBA(0,0,0,.3),
            tickdirection = :out
        )
    end

    # Plot belt L ranges
    plot!(Shape([x_min, x_min, x_max, x_max], [1, 2.5, 2.5, 1]),
        fillcolor = RGBA(0,0,0,.1),
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    plot!(Shape([x_min, x_min, x_max, x_max], [3, 6, 6, 3]),
        fillcolor = RGBA(0,0,0,.1),
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    vline!(obs_edges,
        linewidth = 1,
        linecolor = :black,
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display(plot!())
    end
    return plot!()
end

function absolute_sunset(event::Event; show_plot = true)
    x = event.avg_pitch_angles
    y = log10.(event.energy_bins_mean)
    e_fluence, n_fluence = integrate_flux(event, time = true)
    e_flux = e_fluence ./ event.duration

    n_x = length(x)
    x_min = min(x...)
    x_max = max(x...)
    x_range = 180

    n_y = length(y)
    y_min = min(y...)
    y_max = max(y...)
    y_range = y_max - y_min

    plot(
        title = "Absolute PAD",
        xlabel = "Pitch Angle (deg)",
        ylabel = "Energy (keV)",
        colorbar_title = "Avg. Flux, \nlog10(kev/(cm² str Mev))",
        xlims = (0, 180),
        ylims = (y_min, y_max),
        clims = (3.5, 8),
        yticks = (y[1:2:end], Int.(round.(10 .^ y[1:2:end], sigdigits = 2))),
        aspect_ratio = (x_range / y_range) * 1,
        background_color_inside = :black,
        background_color_outside = :transparent,
        framestyle = :grid,
        grid = false,
        dpi = 300
    )
    heatmap!(x, y, log10.(e_flux),
        colormap = :ice #:cherry
    )
    vline!([event.avg_loss_cone_angle],
        linewidth = 5,
        linecolor = RGB(.5,.5,.5),
        linestyle = :solid,
        label = ""
    )
    vline!([event.avg_anti_loss_cone_angle],
        linewidth = 5,
        linecolor = RGB(.5,.5,.5),
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display(plot!())
    end
    return plot!()
end

function relative_sunset(event::Event; show_plot = true)
    # Normalize fluence to 90º fluence (J/J90º)
    e_fluence, _ = integrate_flux(event, time = true)
    e_flux = e_fluence ./ event.duration
    _, normalizing_idx = findmin(abs.(90 .- event.avg_pitch_angles))
    for e = 1:16 # Normalize each energy bin (row)
        J90 = e_flux[e, normalizing_idx]
        if J90 < 10^5.25
            e_flux[e,:] .= 0
        else
            e_flux[e,:] ./= J90
        end
    end
    e_flux = log10.(e_flux)

    #=
    # Create diverging colorbar centered on 0
    colorbar_range = (-2.5, .5)
    color_anchors = [colorbar_range[1], 0, colorbar_range[2]]
    color_anchors = (color_anchors .- min(color_anchors...)) ./ (max(color_anchors...) - min(color_anchors...)) # Map to [0,1]
    colorbar = cgrad(:vik, color_anchors)
    =#
    colorbar = :ice
    colorbar_range = (-2, 0)

    x = event.avg_pitch_angles
    n_x = length(x)
    x_min = min(x...)
    x_max = max(x...)
    x_range = 180

    y = log10.(event.energy_bins_mean)
    finite_idxs = isfinite.(y)
    y_min = min(y[finite_idxs]...)
    y_max = max(y[finite_idxs]...)
    y_range = y_max - y_min

    plot(
        title = "J/J₉₀ PAD",
        xlabel = "Pitch Angle (deg)",
        ylabel = "Energy (keV)",
        colorbar_title = "log10(J/J₉₀)",
        xlims = (0, 180),
        ylims = (y_min, y_max),
        clims = colorbar_range,
        yticks = (y[1:2:end], Int.(round.(10 .^ y[1:2:end], sigdigits = 2))),
        aspect_ratio = (x_range / y_range) * 1,
        background_color_inside = :black,
        background_color_outside = :transparent,
        framestyle = :grid,
        grid = false,
        dpi = 300
    )
    heatmap!(x, y, e_flux,
        colormap = colorbar
    )
    vline!([event.avg_loss_cone_angle],
        linewidth = 5,
        linecolor = RGB(.5,.5,.5),
        linestyle = :solid,
        label = ""
    )
    vline!([event.avg_anti_loss_cone_angle],
        linewidth = 5,
        linecolor = RGB(.5,.5,.5),
        linestyle = :dash,
        label = ""
    )
    if show_plot == true
        display(plot!())
    end
    return plot!()
end

function MLT_L_track(event::Event; index = false, show_plot = true)
    net_e_flux, _ = integrate_flux(event, energy = true, pitch_angle = true)
    net_e_flux = log10.(net_e_flux)
    net_e_flux[isinf.(net_e_flux)] .= 0

    plot(
        proj = :polar,
        xaxis = false,
        yaxis = false,
        ylims = (0, 8),
        background_color = :transparent,
        grid = false,
        dpi = 300
    )
    # Black background
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, 0), reverse(Plots.partialcircle(0, 2π, 100, 1)))), fillcolor = :black, label = "")
    # L-shells
    for L = 2:2:8
        plot!(LinRange(0,2π,100), zeros(100) .+ L,
            linecolor = RGBA(.5,.5,.5),
            label = ""
        )
    end
    # Radiation belts
    inner_belt = [1, 2.5] ./ 8    # Because plotting circles only works with relative units, we need to get the radiation belt ranges in relative units too
    outer_belt = [3, 6]   ./ 8
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, inner_belt[1]), reverse(Plots.partialcircle(0, 2π, 100, inner_belt[2])))),
        fillcolor = RGBA(1,1,1,.25),
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, outer_belt[1]), reverse(Plots.partialcircle(0, 2π, 100, outer_belt[2])))),
        fillcolor = RGBA(1,1,1,.25),
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    # Earth
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, 0), reverse(Plots.partialcircle(0, 2π, 100, 1/15)))),
        fillcolor = :black,
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    plot!(Shape(vcat(Plots.partialcircle(3π/2, π/2, 100, 0), reverse(Plots.partialcircle(3π/2, π/2, 100, 1/15)))),
        fillcolor = :white,
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
    # Spacecraft path & flux
    MLT_radians = (event.MLT/24) * 2π
    scatter!(MLT_radians, event.L,
        linewidth = 4,
        zcolor = net_e_flux,
        color = :inferno,
        markerstrokewidth = 0,
        label = "",
        colorbar = false,
        #clims = (4,8)
    )
    # MLT Labels
    annotate!(1.4, -.07, text("0000", 12))
    annotate!(-1.35, -.07, text("1200", 12))
    annotate!(0, 1.1, text("0600", 12))
    annotate!(0, -1.1, text("1800", 12))

    if show_plot == true
        display(plot!())
    end
    return plot!()
end

function _display_plot()
    if USE_VSCODE_PLOTS == true
        display("image/png", plot!())
    else
        display(plot!())
    end
end