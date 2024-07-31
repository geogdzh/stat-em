d = 10

## load the emulator
hfile = h5open("data/$(parent_folder)/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
mean_coefs_1 = read(hfile, "mean_coefs_1")
mean_coefs_2 = read(hfile, "mean_coefs_2")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
ens_gmt = read(hfile, "ens_gmt")
close(hfile)

hfile = h5open("data/$(parent_folder)/training_data_ssp585_$(d)d_$(num_ens_members)ens.hdf5", "r")
ens_projts = read(hfile, "projts")[1:d, :, :]
close(hfile)

hfile = h5open("data/ssp119_gmts_$(num_ens_members)ens.hdf5", "r")
ens_gmt_119 = read(hfile, "ens_gmt")
close(hfile)
ens_gmt_119 = mean(ens_gmt_119, dims=1) # the MAXIMUM (turning point) is at index 34

hfile = h5open("data/$(parent_folder)/projts_ssp119_$(d)d_$(num_ens_members)ens.hdf5", "r")
ens_projts_119 = read(hfile, "projts")
close(hfile)

mean_coefs_119_1 = get_mean_coefs(ens_projts_119, ens_gmt_119, degree=1)

#####
#averaging
# year_ens_projts_119 = zeros((d, Int(L2/12), num_ens_members))
# for i in 1:num_ens_members
#     for mode in 1:d
#         year_ens_projts_119[mode,:,i] = month_to_year_avg(ens_projts_119[mode,:,i])
#     end
# end


begin
    fig = Figure(resolution=(1500, 750))
    for mode in 1:10
        row_index = (mode - 1) ÷ 5 + 1
        col_index = (mode - 1) % 5 + 1
        ax = Axis(fig[row_index, col_index], title="Mode $mode")
        # hist!(ax, vec(year_ens_projts_119[mode,1:34,:]), bins=20, color=(:red, 0.5), label="before", normalization=:pdf)
        # hist!(ax, vec(year_ens_projts_119[mode,35:end,:]), bins=20, color=(:blue, 0.5), label="after", normalization=:pdf)

        hist!(ax, vec(ens_projts_119[mode,1:12:34*12+1,:]), bins=20, color=(:red, 0.5), label="before", normalization=:pdf)
        hist!(ax, vec(ens_projts_119[mode,35*12+1:12:end,:]), bins=20, color=(:blue, 0.5), label="after", normalization=:pdf)
        if mode == 5
            axislegend(ax, position=:lt)
        end
    end
    save("figs/$parent_folder/before_after_hist.png", fig)
    
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
    save("figs/$parent_folder/jan_mean_fits_10_modes_ssp_comparison.png", fig)
    display(fig)
end 

#var fits
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
    save("figs/$parent_folder/jan_var_fits_10_modes_ssp_comparison.png", fig)
    display(fig)
end 


#### covariance fits
d = 100
scenario = "ssp585"

hfile = h5open("data/process/chols_$(scenario)_$(d)d.hdf5", "r")
chols = read(hfile, "chols")
chol_coefs = read(hfile, "chol_coefs")
ens_gmt = read(hfile, "ens_gmt")[:]
close(hfile)

hfile = h5open("data/process/covs_$(scenario)_$(d)d.hdf5")
covs = zeros(200, 200, 12, length(ens_gmt)) 
for i in 1:length(ens_gmt)
    covs[:,:,:,i] = read(hfile, "covs_$i")
end

begin
    for shift in [0, 95]
        fig = Figure(resolution=(1000, 1000))
        for i in 1:5
            for j in 1:5
                ax = Axis(fig[i, j], xticklabelrotation=45, title=(i==1 ? "Mode $(j+shift)" : ""), xlabel = (i==5 ? "GMT" : ""), ylabel = (j==1 ? "Mode $(i+shift)" : ""), ylabelfont=:bold)
                scatter!(ax, ens_gmt, covs[i+shift, j+shift, 1, :], color=:orange, alpha=0.5)
                lines!(ax, ens_gmt, covhats[i+shift, j+shift, 1, :], color=:blue)
            end
        end
        save("figs/$parent_folder/covariance_fits_$(shift==0 ? "top" : "bottom").png", fig)
        display(fig)
    end
end
