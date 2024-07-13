using CairoMakie, ProgressBars, HDF5, GeoMakie
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl") 

# fig 1: different scenarios

scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
scenario_colors = Dict("historical" => :red4, "ssp585" => :red, "ssp245" => :magenta3, "ssp119" => :indigo)
M, N = 192, 96
L1, L2 = 1980, 1032 #for CMIP6

using_precip = true 
non_dim = false  
use_metrics = false
if using_precip
    parent_folder = "temp_precip"
else
    parent_folder = "temp"
end
if non_dim
    parent_folder = "nondim"
end
if use_metrics && using_precip
    parent_folder = "metrics"
elseif use_metrics && !using_precip
    parent_folder = "temp_metrics"
end

#get sample timevecs
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
file3 = file_head*"ssp585/tas/r1i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas")
lonvec, latvec = ts3.lonvec[:], ts3.latvec[:]
lonvec2 = lonvec .-180.
time_history = [x for x in 1850:2014]
time_future = [x for x in 2015:2100]

begin
    fig = Figure(resolution=(800, 500))
    ax = Axis(fig[1,1], xlabel="Year", ylabel="GMT (K)", title="Global Mean Temperature in MPI-ESM1.2-LR", xticks=1850:50:2100)
    for scenario in scenarios
        hfile = scenario=="historical" ? h5open("data/$(scenario)_gmts_49ens.hdf5", "r") : h5open("data/$(scenario)_gmts_50ens.hdf5", "r")
        num_ens_members = read(hfile, "num_ens_members")
        ens_gmt = read(hfile, "ens_gmt")
        close(hfile)
        for i in 1:num_ens_members
            lines!(ax, (scenario=="historical" ? time_history : time_future), ens_gmt[i,:], color=scenario_colors[scenario], alpha=0.3, linestyle=:dot)
        end
        mean_gmt = mean(ens_gmt, dims=1)
        lines!(ax, (scenario=="historical" ? time_history : time_future), mean_gmt[1,:], color=scenario_colors[scenario], alpha=1, linestyle=:solid, label=scenario)
    end
    axislegend(ax, position=:lt)
    display(fig)
    # save("figs/gmt_scenarios.png", fig)
end

# fig 2: show the basis
hfile = h5open("data/$(parent_folder)/basis_2000d.hdf5", "r")
basis = read(hfile, "basis")
close(hfile)

#get the lines for the first basis modes
basis_trends_history = zeros((L1, 49, 20, 1))
basis_trends_future = zeros((L2, 50, 20, 3))

hfile = h5open("data/$(parent_folder)/projts_withpr_historical_100d_49ens.hdf5", "r")
projts = read(hfile, "projts")
close(hfile)
# mean_projts = mean(projts, dims=3)[:,:,1]
for i in 1:20
    basis_trends_history[:,:,i,1] = projts[i,:,:]
end
for (j,scenario) in enumerate(scenarios[2:end]) #this except it doesn't really exist
    hfile = h5open("data/temp_precip/projts_withpr_$(scenario)_100d_50ens.hdf5", "r")
    projts = read(hfile, "projts")
    close(hfile)
    for i in 1:20
        basis_trends_future[:,:,i,j] = projts[i,:,:]
    end
end



begin
    fig = Figure(resolution=(1500, 1000)) 
#question here : should it be using uniform coloring? or each their own 
# + should we show yearly averages for the modes? or for one month?
    for (n, i) in enumerate([x for x in 1:3]) #i in 1:3
        ax = GeoAxis(fig[1,n], title="Mode $i", ylabel=(i==1 ? "Temperature" : ""))
        tmp = reshape(basis[1:M*N,i], (M, N))
        heatmap!(ax, lonvec2, latvec, tmp, colormap=:thermometer)

        ex = Axis(fig[2,n], xlabel="Year", ylabel=(i==1 ? "Coefficient of mode (annual mean)" : ""), title="Mode $i")
        for x in 1:49
            lines!(ex, time_history, month_to_year_avg(basis_trends_history[:,x,i,1]), color=scenario_colors["historical"], alpha=0.1, linestyle=:dot)
        end
        mean_trend = month_to_year_avg(mean(basis_trends_history[:,:,i,1], dims=2))
        lines!(ex, time_history, mean_trend, color=scenario_colors["historical"], alpha=1, linestyle=:solid, linewidth=2, label="historical")
        for (j, scenario) in enumerate(scenarios[2:end])
            for x in 1:50
                lines!(ex, time_future, month_to_year_avg(basis_trends_future[:,x,i,j]), color=scenario_colors[scenario], alpha=0.1, linestyle=:dot)
            end
            mean_trend = month_to_year_avg(mean(basis_trends_future[:,:,i,j], dims=2))
            lines!(ex, time_future, mean_trend, color=scenario_colors[scenario], alpha=1, linestyle=:solid, linewidth=2, label=scenario)
        end
        axislegend(ex, position=:lb)

        ax2 = GeoAxis(fig[3,n],  title="Mode $i", ylabel=(i==1 ? "Precipitation" : ""))
        tmp = reshape(basis[M*N+1:end,i], (M, N))
        heatmap!(ax2, lonvec2, latvec, tmp, colormap=:plum)
    end
    # save("figs/basis_modes.png", fig)
    display(fig)
end


# generate visuals of all the modes

# for i in 200:-1:1
#     fig = Figure(resolution=(1000,750))
#     ax = GeoAxis(fig[1,1], title="Mode $i")
#     tmp = reshape(basis[1:M*N,i], (M, N))
#     heatmap!(ax, lonvec2, latvec, tmp, colormap=:thermometer)
#     save("figs/modes_og/temp_mode_$i.png", fig)
#     # display(fig)
# end

#modes 6 and 15 look vaguely enso-like in the non-dim version

### cumulative variance plot
parent_folder = "metrics"

hfile = h5open("data/$(parent_folder)/temp_precip_basis_2000d.hdf5", "r")
# hfile = h5open("data/$(parent_folder)/temp_basis_2000d.hdf5", "r")
basis = read(hfile, "basis")
S = read(hfile, "S")
close(hfile)

ss = sum(S.^2)
cum_var = cumsum(S.^2 ./ ss)

begin
    fig = Figure(resolution=(600, 400))
    ax = Axis(fig[1,1], xlabel="Mode", ylabel="Cumulative variance explained", title="Cumulative variance explained by the modes")
    lines!(ax, cum_var, color=:black, linestyle=:solid)
    display(fig)
    # save("figs/cum_var.png", fig)
end