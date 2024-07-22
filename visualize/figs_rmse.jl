using CairoMakie, ProgressBars, HDF5, GeoMakie
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl") 


#get the latvec to be used later on
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
file3 = file_head*"ssp585/tas/r1i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas")
lonvec, latvec = ts3.lonvec[:], ts3.latvec[:]
lonvec2 = lonvec .-180.

using_two = true 
second_var = "huss" # "pr" or "huss"
non_dim = false  
use_metrics = false
if using_two
    if second_var == "pr"
        parent_folder = "temp_precip"
    else
        parent_folder = "temp_huss"
    end
else
    parent_folder = "temp"
end
if non_dim
    parent_folder = "nondim"
end
if use_metrics && using_two
    parent_folder = "metrics"
elseif use_metrics && !using_two
    parent_folder = "temp_metrics"
end

## 
scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
scenario_colors = Dict("historical" => :red4, "ssp585" => :red, "ssp245" => :magenta3, "ssp119" => :indigo)
time_history = [x for x in 1850:2014]
time_future = [x for x in 2015:2100]
L1, L2 = 1980, 1032 #for CMIP6
l1, l2 = 165, 86


## plots etc 

function plot_rmse(ax, variable, measure, numbers; testing_k=false, rel_error=false)
    variable = variable == "temp" ? "tas" : "pr"

    linestyles = testing_k ? [  :dash, :solid] : [ :dot, :solid]
    for (j, scenario) in enumerate(scenarios[2:end])
        for (i, number) in enumerate(numbers)
            hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$(scenario).hdf5", "r") #toggle "updated" here
            
            if testing_k
                ser = read(hfile, "rmse_$(measure)s_time_$(variable)_k$(number)")
            elseif rel_error
                ser = read(hfile, "rmse_$(measure)s_time_$(variable)_$(number)_rel")
            else    
                ser = read(hfile, "rmse_$(measure)s_time_$(variable)_$(number)")
            end

            lines!(ax, time_future, month_to_year_avg(ser), color=scenario_colors[scenario], alpha=0.6, linestyle=linestyles[i])
            close(hfile)
        end
    end

    elems = testing_k ? [LineElement(color=:black, linestyle=:dash),LineElement(color=:black, linestyle=:solid)] : [LineElement(color=:black, linestyle=:solid), LineElement(color=:black, linestyle=:dot)]
    elems_2 = [ LineElement(color=scenario_colors["ssp119"]), LineElement(color=scenario_colors["ssp245"]), LineElement(color=scenario_colors["ssp585"])]
    labels = testing_k ? ["k=1",  "k=2"] : [ "100 modes",  "10 modes"]
    labels_2 = [ "SSP119", "SSP245", "SSP585"]
    axislegend(ax, [elems..., elems_2...], [labels..., labels_2...], position=(measure=="var" && variable=="temp" ? :lb : :lt))
end


numbers = [10, 100]
ks = [x for x in 1:2]

variable = "pr" #temp/pr
begin 
    fig = Figure(resolution=(1500,1000)) #
    # lims = Dict("temp" => (0.15, 0.6), "pr" => (3e-6, 9e-6))

    rel_error = false # if true, remove ylims settings
    measure = "mean"
    ax = Axis(fig[1,1:4], title="a) Average RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") for varied # of modes", xlabel="Year", ylabel="RMSE")
    # ylims!(ax, lims[variable])
    plot_rmse(ax, variable, measure, numbers; rel_error=rel_error)
    ax = Axis(fig[1,5:8], title="b) Average RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") for varied degree of fit", xlabel="Year")
    # ylims!(ax, lims[variable])
    plot_rmse(ax, variable, measure, ks; rel_error=rel_error, testing_k=true)
    

    ax = GeoAxis(fig[2,1:5], title="d) RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") on SSP119")
    hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$("ssp119").hdf5", "r") #toggle here
    begin
        data =  read(hfile, "rmse_$(measure)s_$(variable == "temp" ? "tas" : variable)_100") #CHANGE BACK
        close(hfile)
        ext = (0., maximum(data))
        heatmap!(ax, lonvec2, latvec, data[:,:,1], colormap=:thermal, colorrange=ext)
        Colorbar(fig[2,6], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    end 

    measure = "std"
    ax = Axis(fig[1,9:12], title="c) Average RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") for varied # of modes", xlabel="Year")
    plot_rmse(ax, variable, measure, numbers; rel_error=rel_error)
    
    ax = GeoAxis(fig[2,7:11], title="e) RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") on SSP119")
    hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$("ssp119").hdf5", "r") #toggle here
    begin
        data = read(hfile, "rmse_$(measure)s_$(variable == "temp" ? "tas" : variable)_100") #CHANGE BACK
        close(hfile)
        ext = (0., maximum(data))
        heatmap!(ax, lonvec2, latvec, data[:,:,1], colormap=:thermal, colorrange=ext)
        Colorbar(fig[2,12], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    end
    colsize!(fig.layout, 12, Relative(1/11))
    # save("figs/$(parent_folder)/rmse_joint_$(variable).png", fig)
    fig
end 
