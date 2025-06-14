using StatsBase
using CairoMakie
using CairoMakie.Colors
using StyledStrings
using Distributions
using Printf
using SideKicks

#export get_results
#export createPlottingProps
#export MasterCornerPlot
#export PlotOverlappingHistograms

##

"""
    struct CornerPlot(fig, ranges, distributions_1d, distributions_2d)

Stores the elements of a CornerPlot to allow easy access to figure elements.
"fig" contains the Makie figure, "ranges" corresponds to the ranges used for each
variable in each plot, as well as to determine credible intervals. The Axes with
1D marginalized distributions are stored in the dictionary "distributions_1d",
using as keys the name of each variable. Similarly, "distributions_2d" provides
the Axes with 2D marginalized distributions as a dictionary of dictionaries. 

"""
struct CornerPlot
    fig
    ranges
    distributions_1d
    distributions_2d
end

function get_results(fname::String)
    results, observations, priors, metadata = SideKicks.ExtractResults(fname, transpose_results=false)
    if (size(results[:P_f])[1] < size(results[:P_f])[2])
        results, observations, priors, metadata = SideKicks.ExtractResults(fname, transpose_results=true)
    end
    return [results, observations, priors, metadata]
end

function createPlottingProps(props_matrix::Vector{Vector{Any}})
    props_matrix = stack(props_matrix) # make into actual matrix
    names   = convert(Vector{Symbol}, props_matrix[1,:])
    scaling = Dict(zip(names, props_matrix[2,:]))
    labels  = Dict(zip(names, props_matrix[4,:]))
    ranges  = Dict()
    for (name, range) in zip(names, props_matrix[3,:])
        if !ismissing(range)
            ranges[name] = range
        end
    end
    return [names, scaling, ranges, labels]
end


function default_range_for_params(param)

    master_default_ranges = Dict(
        :m1_f      => [0,100],  # m_sun
        :m2_f      => [0,100],  # m_sun
        :m1_i      => [0,100],  # m_sun
        :m2_i      => [0,100],  # m_sun
        :dm2       => [0,100],  # m_sun
        :P_f       => [0,1000], # day
        :P_i       => [0,1000], # day
        :P_circ    => [0,1000], # day
        :a_f       => [0,1000], # r_sun
        :a_i       => [0,1000], # r_sun
        :e_f       => [0,1],    # 1
        :e_i       => [0,1],    # 1
        :frac      => [0,1],    # 1
        :Ω_f       => [0,360],  # degree
        :ω_f       => [0,360],  # degree
        :i_f       => [0,90],   # degree
        :Ω_i       => [0,360],  # degree
        :ω_i       => [0,360],  # degree
        :i_i       => [0,90],   # degree
        :ϕ         => [0,360],  # degree
        :θ         => [0,180],  # degree
        :ν_i       => [0,180],  # degree
        :sum_ωi_νi => [0,360],  # degree
        :K1        => [0,300],  # km_per_s
        :K2        => [0,300],  # km_per_s
        :vkick     => [0,300],  # km_per_s
        :vi_α       => [0,300],  # km_per_s
        :vi_δ       => [0,300],  # km_per_s
        :vi_r       => [0,300],  # km_per_s
        :vf_α       => [0,300],  # km_per_s
        :vf_δ       => [0,300],  # km_per_s
        :vf_r       => [0,300],  # km_per_s
        :Δv_α       => [0,300],  # km_per_s
        :Δv_δ       => [0,300],  # km_per_s
        :Δv_r       => [0,300],  # km_per_s
        :Δv       => [0,300],  # km_per_s
        :Venv_α    => [0,300],  # km_per_s
        :Venv_δ    => [0,300],  # km_per_s
        :Venv_r    => [0,300],  # km_per_s
    )
    return master_default_ranges[param]
end





"""
    TODO: clean up the docs here
    MasterCornerPlot(results, names; ...)

Constructs a corner plot from the provided "results", using the variables given
in the "names" vector. 

# Arguments:
- results: Array containing the samples to be plotted. This can be either of type MCMCChains
or a Dictionary containing vectors with values for an individual chain, or arrays for multiple
chains (in this case each column represents a chain). 
- names: Vector of symbols containing the key needed to access each result from `results`.
The corner plot will only include the values specified in `names`.
- labels: Dictionary of strings containing the labels that should be used for each variable.
- ranges: Dictionary of two element vectors, containing the ranges that will be used for each
variable of the plot. If ranges are not provided for a variable these are determined based on the
`quantile_for_range` option.
- scaling: Dictionary containing scaling factors for variables. For any `name` in `names` that
is also a key of `scaling`, all values are divided by `scaling[name]`.
- fig: The Makie figure used for the plot. If not provided it is created.
- quantile_for_range: If ranges are not specified for an axis, then they are set to be between
the quantiles `quantile_for_range` and `1- quantile_for_range`. This is done using weighted
quantiles if `use_weights=true`
- use_weights: If true, then `results[:weights]` is expected to be defined to provide weights
for each sample.
- fraction_1d: Fraction of samples contained in the shown 1D credible intervals. Credible intervals
are determined using highest density intervals. By default 90% credible intervals are shown.
- fractions_2d: Similar to `fraction_1d`, but used to determine the contours in the 2D marginalized
distributions. Values are provided as a Vector of fractions.
- show_CIs: If true, credible intervals are shown in the corner plot.
- nbins: Number of bins in each axis used to plot the heatmaps and the 1D marginalized distributions
- nbins_contour: Number of bins used to plot the contours in the 2D marginalized distributions.
using `nbins_contour<nbins` allows for smoother contous.
- axis_size: The Makie axis will be set to have width and height equal to this value.

# Output:
Returns an instance of CornerPlot
"""

function MasterCornerPlot(
        flist_results, parameters::Vector{Symbol}, labels, scaling, ranges=Dict();
        figsize=10, fig=Figure(size=(figsize,figsize)), axis_size=100, 
        quantile_for_range=0.01, fraction_1D=0.9, fractions_2D=[0.9], 
        use_weights=true, extra_weights=missing, print_CIs=true, show_CIs=true, include_heatmap=true, 
        show_grid=false, show_shaded_bands=true, show_individual_chains=false,
        nbins_1d=100, nbins_contour=20, nbins_heatmap=20,
        colors=missing, linewidths=missing, linestyles=missing,
        invert_histogram_for_params=missing, true_values=missing,
        only_get_xvals_yvals_dict=false, fname_xval_yval_dict="xvals_yvals_dict.jld2",
        latex_strs=missing,
    )

    # =================================================================================
    # The following all must have the same shape, if used:
    #   {parameters, labels, ranges, scaling, invert_histogram_for_params}            
    # The following all must have the same shape, if used:
    #   {flist_results, colors, linewidths, linestyles}  
    # =================================================================================
   
    num_cols = length(parameters)
    num_mods = length(flist_results)
    if ismissing(fig)
        fig = Figure(size=(figsize,figsize))
    end
    if ismissing(colors)
        colors = fill(:black, num_mods)
    end
    if ismissing(linewidths)
        linewidths = fill(1, num_mods)
    end
    if ismissing(linestyles)
        linestyles = fill(:solid, num_mods)
    end
    if ismissing(invert_histogram_for_params)
        invert_histogram_for_params = fill(false, num_cols)
    end
    # Verify the results contain all the parameters
    for param in parameters
        if param ∉ keys(ranges)
            # If range not supplied, get the default range from above
            ranges[param] = default_range_for_params(param)
        end
    end

    # ====================
    # Create all the axes 
    # ====================
    distributions_1d = Dict() # this will store the 1D distributions in a dictionary 
    distributions_2d = Dict() # this will store the 2D distributions in a dictionary of dictionaries
    axes = Dict()
    for ii in 1:num_cols  # ii is the x-coord param
        param_x = parameters[ii]
        axes[param_x] = Dict()
        distributions_2d[param_x] = Dict()

        # Diagonals
        if show_CIs # CIs take up an extra cell above
            axis = Axis(fig[ii+1,ii], xlabel=labels[param_x], height=axis_size, width=axis_size, xgridvisible=show_grid, ygridvisible=show_grid) 
        else
            axis = Axis(fig[ii  ,ii], xlabel=labels[param_x], height=axis_size, width=axis_size, xgridvisible=show_grid, ygridvisible=show_grid) 
        end
        axes[param_x][param_x] = axis
        distributions_1d[param_x] = axis
        
        # Set limits and decorations
        if !invert_histogram_for_params[ii]                                             
            xlims!(axis, (ranges[param_x][1], ranges[param_x][2]))
            hideydecorations!(axis, grid=false)
            if ii != num_cols
                hidexdecorations!(axis, ticks=false, minorticks=false, grid=false)
            end
        else
            ylims!(axis, (ranges[param_x][1], ranges[param_x][2]))
            hidexdecorations!(axis, grid=false)
            hideydecorations!(axis, ticks=false, minorticks=false, grid=false)
        end

        # Off-diagonals
        for jj in ii+1:num_cols # jj is the y-coord param
            if ii >= num_cols
                continue
            end
            param_y = parameters[jj]
            if show_CIs # CIs take up an extra cell above
                axis = Axis(fig[jj+1,ii], xlabel=labels[param_x], ylabel=labels[param_y], height=axis_size, width=axis_size, xgridvisible=show_grid, ygridvisible=show_grid) 
            else
                axis = Axis(fig[jj  ,ii], xlabel=labels[param_x], ylabel=labels[param_y], height=axis_size, width=axis_size, xgridvisible=show_grid, ygridvisible=show_grid) 
            end
            axes[param_x][param_y] = axis
            distributions_2d[param_x][param_y] = axis
            if ii>1
                hideydecorations!(axis, grid=false, ticks=false, minorticks=false)
            end
            if jj!=num_cols
                hidexdecorations!(axis, grid=false, ticks=false, minorticks=false)
            end         
            xlims!(axis, (ranges[param_x][1], ranges[param_x][2]))
            ylims!(axis, (ranges[param_y][1], ranges[param_y][2]))
        end
        
    end

    # ====================
    # Fill all the axes, iterate over the models
    # ====================
    for (kk, fname) in enumerate(flist_results)

        println()
        println(fname)
        results, _, _, _ = get_results(fname)

        
        # Verify the results contain all the parameters
        for param in parameters
            if param ∉ keys(results)
                throw(ArgumentError("$param is not a valid key"))
            end
        end
        # Configure sample_weights
        if :weights ∈ keys(results) && use_weights
            sample_weights = results[:weights]
        else
            # create an array of ones with the needed size
            sample_weights = ones(size(results[parameters[1]]))
        end
        if !ismissing(extra_weights)
            sample_weights .*= extra_weights
        end
        if only_get_xvals_yvals_dict 
            # put the xval_yval_dict in the same dir as the results file
            path_dict = join(split(fname, '/')[1:end-1], '/')*'/'*fname_xval_yval_dict
            if isfile(path_dict)
                println("xval_yval_dict already exists for: "* path_dict)
                continue
            else
                println(path_dict)
                xvals_yvals_dict = Dict()
            end
        end

 
        # =================================
        # Create 1D PDFs along the diagonal
        # =================================

        for ii in 1:num_cols
            param_x = parameters[ii]
            axis = axes[param_x][param_x]
            if param_x ∉ keys(scaling)
                values_x = results[param_x]
            else
                values_x = results[param_x]./scaling[param_x]
            end
            (xmin, xmode, xmax), x, h, y, dx, frac_lost =
                plot_compound_1D_density(axis, param_x, values_x, ranges[param_x], sample_weights, 
                                            fraction_1D, nbins_1d, show_individual_chains=show_individual_chains,
                                            show_shaded_bands=show_shaded_bands, invert_histogram=invert_histogram_for_params[ii],
                                            color=colors[kk], linewidth=linewidths[kk], linestyle=linestyles[kk])

            # Configure confidence intervals
            xupp = xmax-xmode
            if xupp < 1e-3
                xupp = 0.0
            end
            xlow = xmode-xmin
            if xlow < 1e-3
                xlow = 0.0
            end
            # TODO: make the .2 a user-set variable...
            if ismissing(latex_strs)
                str_xmode = @sprintf("%.2f", xmode)
                str_xupp = @sprintf("%.2f", xupp)
                str_xlow = @sprintf("%.2f", xlow)
                latex_bounds = L"%$(str_xmode)^{+%$(str_xupp)}_{-%$(str_xlow)}"
            else
                latex_bounds = latex_strs[ii]
            end

            if print_CIs
                println(labels[param_x]*"="*latex_bounds)
            end
            if show_CIs
                Label(fig[ii,ii], latex_bounds, valign=:bottom)
                axis = Axis(fig[ii,ii])
                hidedecorations!(axis)
                hidespines!(axis)
            end
            # If only interested in xvals_yvals_dict, store values and skip the rest
            if only_get_xvals_yvals_dict
                xvals_yvals_dict[string(param_x)] = (x, y, latex_bounds)
                continue
            end
            GC.gc()
        end     
        resize_to_layout!(fig)
        GC.gc()
        if only_get_xvals_yvals_dict   
            save(path_dict, xvals_yvals_dict)
            continue # skip 2D plots if you only want histogram output for xvals_yvals_dict
        end

        # =======================
        # Create 2D density plots
        # =======================

        for ii in 1:num_cols-1         # ii is the x-coord param
            param_x = parameters[ii]
            for  jj in ii+1:num_cols   # jj is the y-coord param
                param_y = parameters[jj]
                axis = axes[param_x][param_y]

                if param_x ∉ keys(scaling)
                    values_x = results[param_x]
                else
                    values_x = results[param_x]./scaling[param_x]
                end
                if param_y ∉ keys(scaling)
                    values_y = results[param_y]
                else
                    values_y = results[param_y]./scaling[param_y]
                end

                plot_2D_density(axis, param_x, param_y, vec(values_x), ranges[param_x], vec(values_y), ranges[param_y], 
                                   vec(sample_weights), fractions_2D, include_heatmap=include_heatmap, 
                                   nbins_heatmap=nbins_heatmap, nbins_contour=nbins_contour, 
                                   color=colors[kk], linewidth=linewidths[kk], linestyle=linestyles[kk])
                GC.gc()
            end  
            GC.gc()
        end
    end
        
    # ===============
    # Add true values
    # ===============
    if !ismissing(true_values)
        # true value line configs
        color=:black
        linestyle=:dash
        linewidth=1

        for ii in 1:length(parameters)
            param_x = parameters[ii]
            if param_x ∉ keys(true_values)
                println("$param_x is not in true_values dict")
                continue
            end
            axis = axes[param_x][param_x]
            if !invert_histogram_for_params[ii]
                vlines!(axis, true_values[param_x], color=color, linestyle=linestyle, linewidth=linewidth)
            else
                hlines!(axis, true_values[param_x], color=color, linestyle=linestyle, linewidth=linewidth)
            end
            for jj in 1:ii-1
                param_y = parameters[jj]
                if param_y ∉ keys(true_values)
                    continue
                end
                axis = axes[param_y][param_x]
                hlines!(axis, true_values[param_x], color=color, linestyle=linestyle, linewidth=linewidth)
                vlines!(axis, true_values[param_y], color=color, linestyle=linestyle, linewidth=linewidth)
                #if add_point
                #    scatter!(ax, true_vals[param2], true_vals[param1], color=:red, markersize=10)
                #end
            end
        end
    end

    if !only_get_xvals_yvals_dict   
        return CornerPlot(fig, ranges, distributions_1d, distributions_2d)
    end
end


### HELPER FUNCTIONS BELOW ###

"""
    get_bounds_for_fractions(weights, fractions)

Calculate the probability values that corresponds to a given set of highest density intervals.

# Arguments:
- h: `h` is expected to be a vector or array containing the values of a PDF
within a regular grid. This can be either a 1D vector for 1D marginalized distributions,
or a 2D vector for the 2D marginalized distributions. Can also be used for higher dimensional
cases though.
- h is now weights! (what used to be h.weights)
- fractions: Vector containing the individual fractions that are contained within the HDIs.
for example, `fractions=[0.5,0.9]` means that the value of the PDF corresponding to the 50%
and 90% HDIs will be computed.

# Output:
- bounds: Vector containing the values of the PDF corresponding to the HDIs that contain
the fraction of the samples given by `fractions`.
"""
function get_bounds_for_fractions(weights, fractions)
    integral = sum(weights)
    bounds = zeros(length(fractions))
    for jj in eachindex(fractions)
        fraction = fractions[jj]
        minbound = 0
        maxbound = maximum(weights)
        newbound = 0
        # find HDI 
        for ii in 1:15
            newbound = 0.5*(minbound+maxbound)
            #print("  nb = ")
            #println(newbound)
            integral2 = sum(weights[weights.>newbound])
            newfraction = integral2/integral
            if newfraction>fraction
                minbound = newbound
            else
                maxbound = newbound
            end
        end
        bounds[jj] = newbound
    end 
    return sort(bounds)
end


"""
    plot_2D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

Creates a 1D marginalized distribution plot in `axis` from the sample values given in `values`.
The plotted line is normalized, such that it corresponds to the PDF followed by the samples.

# Arguments:
- axis: The Axis on which the distribution will be plotted. If set to `nothing`, no plot will be made,
but the values used for the plot will still be returned.
- name_x: Symbol containing the name for the variable on the x-axis. Only used to report if too
many samples are outside the range to determine a good HDI.
- name_y: Same as `name_x` but for the y-axis.
- values_x: Values for each sample corresponding to the variable `name_x`.
Either a vector or a matrix of numbers, corresponding to a single chain or multiple chains
respectively (in the latter case, columns represent individual samples). The heatmap will include
values from all chains.
- values_y: Same as `values_x` but for the y-axis.
- range_x: Vector with two values containing the range of the plot for the x-axis variable. It also corresponds to the range
that will be used to determine the HDI. The HDI does correct for values outside the range, assuming
that these do not fall within the highest density region.
- range_y: Same as `range_x` but for the y-axis.
- sample_weights: weights used for each sample, should have the same size as `values`.
- fractions: Contours will be plotted corresponding the HDIs containing each fraction contained in fractions.
- nbins: Number of bins used for the heatmap.
- nbins_contours: Number of bins used to determine the contours for the HDIs.

# Output:
- x_hm: Bin centers used in the heatmap in the x-axis
- y_hm: Bin centers used in the heatmap in the y-axis
- z_hm: Value of the PDF estimated at each bin of the heatmap
- dx_hm: Size of bins used for heatmap in the x-axis
- dy_hm: Size of bins used for heatmap in the y-axis
- x_ct: Same as `x_hm` but for the contours.
- y_ct: Same as `y_hm` but for the contours.
- z_ct: Same as `z_hm` but for the contours.
- dx_ct: Same as `dx_hm` but for the contours.
- dy_ct: Same as `dy_hm` but for the contours.
- bounds: probablity values corresponding to the HDIs
- frac_lost: Fraction of samples (including weights) that is outside of the ranges
"""
function plot_2D_density(axis, name_x, name_y, values_x, range_x, values_y, range_y, sample_weights, fractions; nbins_heatmap, nbins_contour=-1, color=(:black, 1), linewidth=1, linestyle=:solid, include_heatmap=true)

    # bounds_99percent are for the underlying plot. range just shows whats within the window...
    filter = (values_x .> range_x[1]) .&& (values_x .< range_x[2]) .&& (values_y .> range_y[1]) .&& (values_y .> range_y[1])
    total_weight = sum(sample_weights)
    missing_weight = sum(sample_weights[.!filter])
    frac_lost = missing_weight/total_weight
    sample_weights = weights(sample_weights)
    correction = 1
    
    # Add in heatmap
    if include_heatmap
        edges = (LinRange(range_x[1], range_x[2], nbins_heatmap+1),
                 LinRange(range_y[1], range_y[2], nbins_heatmap+1))
        h_hm = fit(Histogram, (values_x[filter], values_y[filter]), weights(sample_weights[filter]), edges)#, nbins=nbins_heatmap) 
        x_hm = (h_hm.edges[1][2:end] .+ h_hm.edges[1][1:end-1])./2
        y_hm = (h_hm.edges[2][2:end] .+ h_hm.edges[2][1:end-1])./2
        dx_hm = x_hm[2]-x_hm[1]
        dy_hm = x_hm[2]-x_hm[1]
        z_hm = h_hm.weights/(total_weight*dx_hm*dy_hm)
        if !isnothing(axis)
            heatmap!(axis, x_hm, y_hm, z_hm, colormap=:dense, colorrange=(0, maximum(z_hm)))
        end
    end

    # Add in contours
    binwidth_x = (range_x[2] - range_x[1])/nbins_contour
    binwidth_y = (range_y[2] - range_y[1])/nbins_contour

    # nbins_contour should be bins desired according to the block shape
    candidate_bounds_x = [quantile(values_x, sample_weights, 0.005), quantile(values_x, sample_weights, 0.995)]
    candidate_bounds_y = [quantile(values_y, sample_weights, 0.005), quantile(values_y, sample_weights, 0.995)] 
    # If "close to 0", just set to 0 - check this
    if candidate_bounds_x[1] < 0.01*candidate_bounds_x[2]
        candidate_bounds_x[1] = 0.0
    end
    # If "close to 0", just set to 0 - check this
    if candidate_bounds_y[1] < 0.01*candidate_bounds_y[2]
        candidate_bounds_y[1] = 0.0
    end
    bounds_x = (min(candidate_bounds_x[1], range_x[1]),
                max(candidate_bounds_x[2], range_x[2]))
    bounds_y = (min(candidate_bounds_y[1], range_y[1]),
                max(candidate_bounds_y[2], range_y[2]))

    nbins_x = Int(ceil((bounds_x[2] - bounds_x[1])/binwidth_x))
    nbins_y = Int(ceil((bounds_y[2] - bounds_y[1])/binwidth_y))
    edges = (LinRange(bounds_x[1], bounds_x[2], nbins_x+1),
             LinRange(bounds_y[1], bounds_y[2], nbins_y+1))

    h_ct = fit(Histogram, (values_x, values_y), sample_weights, edges)
    x_ct = (h_ct.edges[1][2:end] .+ h_ct.edges[1][1:end-1])./2
    y_ct = (h_ct.edges[2][2:end] .+ h_ct.edges[2][1:end-1])./2
    dx_ct = x_ct[2]-x_ct[1]
    dy_ct = y_ct[2]-y_ct[1]
    z_ct = h_ct.weights/(total_weight*dx_ct*dy_ct)

    bounds = get_bounds_for_fractions(z_ct, fractions./correction)
    if !isnothing(axis)
        contour!(axis, x_ct, y_ct, z_ct, levels=bounds, color=color, linewidth=linewidth, linestyle=linestyle)
    end

    return x_ct, y_ct, dx_ct, dy_ct, bounds, frac_lost

end  


"""
    plot_1D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

Creates a 1D marginalized distribution plot in `axis` from the sample values given in `values`.
The plotted line is normalized, such that it corresponds to the PDF followed by the samples.

# Arguments:
- axis: The Axis on which the distribution will be plotted. If set to `nothing`, no plot will be made,
but the values used for the plot will still be returned.
- values: Either a vector or a matrix of numbers, corresponding to a single chain or multiple chains
respectively (in the latter case, columns represent individual samples).
- range: Vector with two values containing the range of the plot. It also corresponds to the range
that will be used to determine the HDI. The HDI does correct for values outside the range, assuming
that these do not fall within the highest density region.
- chain_weights: weights used for each sample, should have the same size as `values`.
- fraction_1d: Fraction of samples (including weights) condained within the HDI for which the
credible interval is reported.
- nbins: Number of bins used to build the histogram and determine the HDI
- color: color used for the line plot.
- linewidth: linewidth used for the line plot

# Output:
- x: Bin centers used in the histogram
- h: Output of `Histogram` fit. This is not normalized, so it contains the sum of all weights
(or count of all samples) within each bin.
- y: Value of the PDF estimated at each bin
- frac_lost: Fraction of samples (including weights) that is outside of `range`
"""
function plot_1D_density(axis, values, range, chain_weights, nbins_1d; color, linewidth, linestyle, invert_histogram=false)

    # bounds_99percent are for the underlying plot. range just shows whats within the window...
    filter = values .> range[1] .&& values .< range[2]
    total_weight = sum(chain_weights)
    missing_weight = sum(chain_weights[.!filter])
    frac_lost = missing_weight/total_weight
    #chain_weights = weights(chain_weights) # weights is a StatsBase function

    # new things
    values = values[filter]
    chain_weights = weights(chain_weights[filter]) # weights is a StatsBase function

    binwidth_x = (range[2] - range[1])/nbins_1d
    candidate_bounds_x = [quantile(values, chain_weights, 0.005), quantile(values, chain_weights, 0.995)]
    # If "close to 0", just set to 0 - check this
    if candidate_bounds_x[1] < 0.01*candidate_bounds_x[2]
        candidate_bounds_x[1] = 0.0
    end
    bounds_x = (min(candidate_bounds_x[1], range[1]),
                max(candidate_bounds_x[2], range[2]))
    nbins_x = Int(ceil((bounds_x[2] - bounds_x[1])/binwidth_x))
    edges = LinRange(bounds_x[1], bounds_x[2], nbins_x+1)
    h = fit(Histogram, values, chain_weights, edges)

    x =(h.edges[1][2:end] .+ h.edges[1][1:end-1])./2
    dx = x[2]-x[1]
    y = h.weights /(total_weight*dx) # normalize
    if !isnothing(axis)
        if !invert_histogram
            lines!(axis, x, y, color=color, linewidth=linewidth, linestyle=linestyle)
        else
            lines!(axis, y, x, color=color, linewidth=linewidth, linestyle=linestyle)
        end
    end

    return x, h, y, dx, frac_lost
end


"""
    plot_compund_1D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

Creates a 1D marginalized distribution plot in `axis` from the sample values given in `values`.
If values contains more than one chain, it will plot the individual chains separetely, as well
as one line plot containing all of the samples from all chains.

# Arguments:
- axis: The Axis on which the distribution will be plotted. If not provided, no plot will be made,
but the values used for the plot will still be returned.
- name: Symbol containing the name of the variable.
- values: Either a vector or a matrix of numbers, corresponding to a single chain or multiple chains
respectively (in the latter case, columns represent individual samples).
- range: Vector with two values containing the range of the plot. It also corresponds to the range
that will be used to determine the HDI. The HDI does correct for values outside the range, assuming
that these do not fall within the highest density region.
- sample_weights: weights used for each sample, should have the same size as `values`.
- fraction_1d: Fraction of samples (including weights) condained within the HDI for which the
credible interval is reported.
- nbins: Number of bins used to build the histogram and determine the HDI

# Output:
- (xmin, xmode, xmax): Credible interval edges and mode.
- x: Bin centers used in the histogram. This will include all samples from all chains.
- h: Output of `Histogram` fit. This is not normalized, so it contains the sum of all weights
(or count of all samples) within each bin.
- y: Value of the PDF estimated at each bin
- frac_lost: Fraction of samples (including weights) that is outside of `range`
"""
function plot_compound_1D_density(axis, name, values_x, range_x, sample_weights, fraction_1D, nbins_1d; color=(:black, 1), linewidth=1, linestyle=:solid, show_individual_chains=false, show_shaded_bands=false, invert_histogram=false)

    # Plot once for all the values
    x, h, y, dx, frac_lost = plot_1D_density(axis, vec(values_x), range_x, vec(sample_weights), nbins_1d, color=color, linewidth=linewidth, linestyle=linestyle, invert_histogram=invert_histogram)

    # If we have multiple chains, iterate over them
    if show_individual_chains
        if ndims(values_x)==2 && !isnothing(axis)
            for ii in 1:size(values_x)[2] # nchains
                chain_values = @view values_x[:,ii]
                chain_weights = @view sample_weights[:,ii]
                if sum(chain_weights) == 0
                    continue
                end
                plot_1D_density(axis, chain_values, range_x, chain_weights, nbins_1d, color=(:gray, 0.85), linewidth=1, linestyle=:solid, invert_histogram=invert_histogram)
            end
        end
    end

    effective_fraction = [fraction_1D * (1/(1-frac_lost))] # need to rescale to the remaining part of the distribution
    bound = get_bounds_for_fractions(h.weights, effective_fraction)[1]
    xmin = minimum(x[h.weights .>= bound]) - dx/2 # get left most value of bin
    xmax = maximum(x[h.weights .>= bound]) + dx/2
    xmode = x[argmax(h.weights)]
    if xmode == x[1]
        xmode = range_x[1]
    elseif xmode == x[end]
        xmode = range_x[2]
    end

    if show_shaded_bands
        filter = x .>= xmin .&& x .<= xmax
        if !invert_histogram
            vlines!(axis, xmode, color=(:black, 1.0), linewidth=1)
            xlims!(axis, range_x[1], range_x[2])
            ylims!(axis, 0, 1.1*maximum(y))
            band!(axis, x[filter], zeros(length(x[filter])), y[filter], color=(:gray, 0.35)) #HSV(0, 0, 0.85)) 
        else
            hlines!(axis, xmode, color=(:black, 1.0), linewidth=1)
            ylims!(axis, range_x[1], range_x[2])
            xlims!(axis, 0, 1.1*maximum(y))

            # Need some fuckery for the horizontal bands, this doesn't exist in Cairomakie
            # Need to identify all the points that make up the two curves, interlace them,
            # And plot a bunch of filled triangles between them. Hacky as shit, but good enough for now.
            
            npoints = size(y[filter])[1]
            points = Array{Float64}(undef, 2, 2*npoints)
            buffer = 0.99
            for ii in 1:npoints
                point1 = [0 ,                   x[filter][ii]]
                point2 = [y[filter][ii]*buffer, x[filter][ii]] 
                points[:,2*ii-1] = point1
                points[:,2*ii]   = point2
            end
            for ii in 1:2*npoints-3
                triplet = Point2f[
                    (points[1,ii], points[2,ii]),
                    (points[1,ii+1], points[2,ii+1]),
                    (points[1,ii+2], points[2,ii+2]),
                ]
                poly!(axis, triplet, strokewidth=1, strokecolor=HSV(0, 0, 0.80)) 
            end

        end
    end
    # Plot once for all the values
    x, h, y, dx, frac_lost = plot_1D_density(axis, vec(values_x), range_x, vec(sample_weights), nbins_1d, color=color, linewidth=linewidth, linestyle=linestyle, invert_histogram=invert_histogram)

    return (xmin, xmode, xmax), x, h, y, dx, frac_lost
   
end





function PlotOverlappingHistograms( 
        flist_datadict, parameters::Vector{Symbol}, labels, ranges;
        figsize=10, fig=Figure(size=(figsize,figsize)), axis_size=100, ncols=3,
        linelabels=missing, colors=missing, linewidths=missing, linestyles=missing)

    function get_row_col(ii)
        row = ((ii-1) % ncols) + 1 
        col = floor(Int, (ii-1) / ncols) + 1 + 1  # extra +1 is to make room for the legend at the top
        return (row, col)
    end
    for (ii, param) in enumerate(parameters)
        (row, col) = get_row_col(ii)
        println(row, " ", col)
        # TODO: swap row and col?
        ax = Axis(fig[col, row], xlabel=labels[param], height=axis_size, width=axis_size, yticklabelsvisible=false, yticksvisible=true)
        if param ∈ keys(ranges)
            xlims!(ax, ranges[param][1], ranges[param][2])
        end
        for (jj, fdatadict) in enumerate(flist_datadict)
            data_dict_1D_hists = load(fdatadict) 
            (xvals, yvals, latex_str) = data_dict_1D_hists[String(param)]
            lin = lines!(ax, xvals, yvals, color=colors[jj], linestyle=linestyles[jj], linewidth=linewidths[jj])
            # Add text of the 90% CIs
            xpos = 0.5
            ypos = jj/6
            if jj == 5
                ypos = 0
            end
            text!(xpos, ypos + .16 , text=latex_str, space=:relative, color=colors[jj], fontsize=14)
        end
        GC.gc()
    end 


    return fig
end

