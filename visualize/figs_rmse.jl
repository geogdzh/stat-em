

function plot_rmse(ax, variable, measure, numbers; testing_k=false, rel_error=false)
    label = variable == "tas" ? variable : "two"

    linestyles = testing_k ? [  :dash, :solid] : [ :dot, :solid]
    for (j, scenario) in enumerate(scenarios[2:end])
        for (i, number) in enumerate(numbers)
            hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$(scenario).hdf5", "r")
            
            if testing_k
                ser = read(hfile, "rmse_$(measure)s_time_$(label)_k$(number)")
            elseif rel_error
                ser = read(hfile, "rmse_$(measure)s_time_$(label)_$(number)_rel")
            else    
                ser = read(hfile, "rmse_$(measure)s_time_$(label)_$(number)")
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


function generate_rmse_fig(variable)
    numbers = [10, 100]
    ks = [x for x in 1:2]
    label = variable == "tas" ? variable : "two"
    long_label = var_labels[variable]    

    begin 
        fig = Figure(resolution=(1500,1000)) #
        # lims = Dict("temp" => (0.15, 0.6), "pr" => (3e-6, 9e-6))

        rel_error = false # if true, remove ylims settings
        measure = "mean"
        ax = Axis(fig[1,1:4], title="a) Average RMSE of the ensemble $(measure) \n for $(long_label) for varied # of modes", xlabel="Year", ylabel="RMSE")
        # ylims!(ax, lims[variable])
        plot_rmse(ax, variable, measure, numbers; rel_error=rel_error)
        ax = Axis(fig[1,5:8], title="b) Average RMSE of the ensemble $(measure) \n for $(long_label) for varied degree of fit", xlabel="Year")
        # ylims!(ax, lims[variable])
        plot_rmse(ax, variable, measure, ks; rel_error=rel_error, testing_k=true)
        

        ax = GeoAxis(fig[2,1:5], title="d) RMSE of the ensemble $(measure) \n for $(long_label) on SSP119", xticklabelrotation=45.0)
        hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$("ssp119").hdf5", "r")
        begin
            data =  read(hfile, "rmse_$(measure)s_$(label)_100") 
            close(hfile)
            ext = (0., maximum(data))
            heatmap!(ax, lonvec2, latvec, data[:,:,1], colormap=:thermal, colorrange=ext)
            Colorbar(fig[2,6], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
        end 

        measure = "std"
        ax = Axis(fig[1,9:12], title="c) Average RMSE of the ensemble $(measure) \n for $(long_label) for varied # of modes", xlabel="Year")
        plot_rmse(ax, variable, measure, numbers; rel_error=rel_error)
        
        ax = GeoAxis(fig[2,7:11], title="e) RMSE of the ensemble $(measure) \n for $(long_label) on SSP119", xticklabelrotation=45.0)
        hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$("ssp119").hdf5", "r")
        begin
            data = read(hfile, "rmse_$(measure)s_$(label)_100") 
            close(hfile)
            ext = (0., maximum(data))
            heatmap!(ax, lonvec2, latvec, data[:,:,1], colormap=:thermal, colorrange=ext)
            Colorbar(fig[2,12], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
        end
        colsize!(fig.layout, 12, Relative(1/11))
        save("figs/$(parent_folder)/rmse_joint_$(variable).png", fig)
        display(fig)
    end 
end