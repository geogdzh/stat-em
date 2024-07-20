using CairoMakie, ProgressBars, HDF5
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl")

scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
scenario_colors = Dict("historical" => :red4, "ssp585" => :red, "ssp245" => :magenta3, "ssp119" => :indigo)
time_history = [x for x in 1850:2014]
time_future = [x for x in 2015:2100]
L1, L2 = 1980, 1032 #for CMIP6
l1, l2 = 165, 86

d=10
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

## load the emulator
hfile = h5open("data/$(parent_folder)/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
mean_coefs_1 = read(hfile, "mean_coefs_1")
mean_coefs_2 = read(hfile, "mean_coefs_2")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
num_ens_members = read(hfile, "num_ens_members")
ens_gmt = read(hfile, "ens_gmt")
close(hfile)

hfile = h5open("data/$(parent_folder)/training_data_ssp585_$(d)d_49ens.hdf5", "r")
ens_projts = read(hfile, "projts")[1:d, :, :]
# ens_gmt = read(hfile, "ens_gmt")
# num_ens_members = read(hfile, "num_ens_members")
close(hfile)

hfile = h5open("data/ssp119_gmts_50ens.hdf5", "r")
ens_gmt_119 = read(hfile, "ens_gmt")
close(hfile)
ens_gmt_119 = mean(ens_gmt_119, dims=1) # the MAXIMUM (turning point) is at index 34

hfile = h5open("data/$(parent_folder)/projts_ssp119_$(d)d_50ens.hdf5", "r")
ens_projts_119 = read(hfile, "projts")
close(hfile)

mean_coefs_119_1 = get_mean_coefs(ens_projts_119, ens_gmt_119, degree=1)
# mean_coefs_119_2 = get_mean_coefs(ens_projts_119, ens_gmt_119, degree=2)

#####
#averaging
year_ens_projts_119 = zeros((d, Int(L2/12), num_ens_members))
for i in 1:num_ens_members
    for mode in 1:d
        year_ens_projts_119[mode,:,i] = month_to_year_avg(ens_projts_119[mode,:,i])
    end
end


begin
    fig = Figure(resolution=(1500, 750))
    for mode in 1:10
        row_index = (mode - 1) ÷ 5 + 1
        col_index = (mode - 1) % 5 + 1
        ax = Axis(fig[row_index, col_index], title="Mode $mode")
        # for i in 1:50
            # scatter!(ax, ens_gmt_119[1:34], month_to_year_avg(ens_projts_119[mode,:,i])[1:34], color=:red, label="before")
            # scatter!(ax, ens_gmt_119[35:end], month_to_year_avg(ens_projts_119[mode,:,i])[35:end], color=:blue, label="after")
        # end
        hist!(ax, vec(year_ens_projts_119[mode,1:34,:]), bins=20, color=(:red, 0.5), label="before", normalization=:pdf)
        hist!(ax, vec(year_ens_projts_119[mode,35:end,:]), bins=20, color=(:blue, 0.5), label="after", normalization=:pdf)
        if mode == 5
            axislegend(ax, position=:rt)
        end
    end
    # save("figs/$parent_folder/before_after_hist.png", fig)
    
    display(fig)
end


#mean fits 
begin
    fig = Figure(resolution=(1500, 750))
    mon = 1
    for mode in 1:10
        row_index = (mode - 1) ÷ 5 + 1
        col_index = (mode - 1) % 5 + 1
        ax = Axis(fig[row_index, col_index], title="Mode $mode")
        for i in 1:num_ens_members
            scatter!(ax, ens_gmt[1:l1], ens_projts[mode, mon:12:L1, i], markersize=5, alpha=0.3, color=scenario_colors["historical"], label="historical") 
            scatter!(ax, ens_gmt[l1+1:end], ens_projts[mode, L1+mon:12:end, i], markersize=5, alpha=0.3, color=scenario_colors["ssp585"], label="ssp585") 
            scatter!(ax, ens_gmt_119[:], ens_projts_119[mode, mon:12:end, i], markersize=5, alpha=0.3, color=scenario_colors["ssp119"], label= "ssp119")
            # scatter!(ax, ens_gmt[1:l1], month_to_year_avg(ens_projts[mode, 1:L1, i]), markersize=5, alpha=0.3, color=scenario_colors["historical"]) 
            # scatter!(ax, ens_gmt[l1+1:end], month_to_year_avg(ens_projts[mode, L1+1:end, i]), markersize=5, alpha=0.3, color=scenario_colors["ssp585"]) 
            # scatter!(ax, ens_gmt_119[:], month_to_year_avg(ens_projts_119[mode, :, i]), markersize=5, alpha=0.3, color=scenario_colors["ssp119"]) 
            if mode == 1 && i == 1
                axislegend(ax, position=:rt)
            end
        end
        ens_mean = mean(ens_projts[mode, mon:12:end, :], dims=2)
        scatter!(ax, ens_gmt[:], ens_mean[:], color=:black, markersize=6, label="ensemble mean")
        lines!(ax, ens_gmt[:], [mean_coefs_1[mon, mode, 2].*x.+mean_coefs_1[mon, mode, 1] for x in ens_gmt[:]], color=:black, linewidth=3, label="linear fit")
        lines!(ax, ens_gmt[:], [mean_coefs_2[mon, mode, 3].*x.^2 .+ mean_coefs_2[mon, mode, 2].*x.+mean_coefs_2[mon, mode, 1] for x in ens_gmt[:]], color=:blue, linewidth=3, label="quadratic fit")

        if mode == 1
            elems = [MarkerElement(color=:black, marker=:marker, markersize=6), LineElement(color=:black, linestyle=:solid, linewidth=3), LineElement(color=:blue, linestyle=:solid, linewidth=3)]
            labels = ["ensemble\nmean", "linear fit", "quadratic fit"]
            axislegend(ax, elems, labels, position=:lb)
        end
    end
    # save("figs/$parent_folder/jan_mean_fits_10_modes_ssp_comparison.png", fig)
    display(fig)
end 

#var fits
#needs: vars, num_ens_members, var_coefs :: after an emulator has been trained

var_coefs, vars = get_var_coefs(ens_projts, ens_gmt, mean_coefs_1; return_vars=true)
var_coefs_119, vars_119 = get_var_coefs(ens_projts_119, ens_gmt_119, mean_coefs_119_1; return_vars=true)


begin
    fig = Figure(resolution=(1500, 750))
    mon = 1
    for mode in 1:10
        row_index = (mode - 1) ÷ 5 + 1
        col_index = (mode - 1) % 5 + 1
        ax = Axis(fig[row_index, col_index], title="Mode $mode")
        for i in 1:num_ens_members
            scatter!(ax, ens_gmt[1:l1], sqrt.(vars[mon, mode, i, 1:l1]), markersize=5, alpha=0.3, color=scenario_colors["historical"], label="historical") 
            scatter!(ax, ens_gmt[l1+1:end], sqrt.(vars[mon, mode, i, l1+1:end]), markersize=5, alpha=0.3, color=scenario_colors["ssp585"], label="ssp585") 
            scatter!(ax, ens_gmt_119[:], sqrt.(vars_119[mon, mode, i, :]), markersize=5, alpha=0.3, color=scenario_colors["ssp119"], label="ssp119")
            
            if mode == 1 && i == 1
                axislegend(ax, position=:rt)
            end
        end
        ens_mean = sqrt.(mean(dropdims(vars[mon:12:end, mode, :, :], dims=1), dims=1))
        scatter!(ax, ens_gmt[:], ens_mean[:], color=:black, markersize=7)
        lines!(ax, ens_gmt[:], sqrt.([var_coefs[mon, mode, 2].*x.+var_coefs[mon, mode, 1] for x in ens_gmt[:]]), color=:black, linewidth=3)
        # (actually plotting standard deviation)

        if mode == 2
            elems = [MarkerElement(color=:black, marker=:marker, markersize=6), LineElement(color=:black, linestyle=:solid, linewidth=3)]
            labels = ["ensemble\nmean", "linear fit"]
            axislegend(ax, elems, labels, position=:rt)
        end
    end
    # save("figs/$parent_folder/jan_var_fits_10_modes_ssp_comparison.png", fig)
    display(fig)
end 
