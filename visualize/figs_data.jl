# fig 1: different scenarios

fig = Figure(resolution=(600, 500))
ax = Axis(fig[1,1], xlabel="Year", ylabel="Global Mean Temperature (K)", title="Global Mean Temperature in MPI-ESM1.2-LR", xticks=1850:50:2100)
for scenario in scenarios
    hfile = h5open("data/$(scenario)_gmts_$(num_ens_members)ens.hdf5", "r")
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
save("figs/$parent_folder/gmt_scenarios.png", fig)


# fig 2: show the basis
hfile = h5open("data/$(parent_folder)/basis_2000d.hdf5", "r")
basis = read(hfile, "basis")
close(hfile)

#get the lines for the first basis modes
basis_trends_history = zeros((L1, num_ens_members, 10, 1))
basis_trends_future = zeros((L2, num_ens_members, 10, 3))

hfile = h5open("data/$(parent_folder)/projts_historical_10d_$(num_ens_members)ens.hdf5", "r")
projts = read(hfile, "projts")
close(hfile)
# mean_projts = mean(projts, dims=3)[:,:,1]
for i in 1:10
    basis_trends_history[:,:,i,1] = projts[i,:,:]
end
for (j,scenario) in enumerate(scenarios[2:end]) #this except it doesn't really exist
    hfile = h5open("data/$(parent_folder)/projts_$(scenario)_10d_$(num_ens_members)ens.hdf5", "r")
    projts = read(hfile, "projts")
    close(hfile)
    for i in 1:10
        basis_trends_future[:,:,i,j] = projts[i,:,:]
    end
end



begin
    fig = Figure(resolution=(1500, 1000)) 
    for (n, i) in enumerate([x for x in 1:3]) #i in 1:3
        ax = GeoAxis(fig[1,n], title="Mode $i", ylabel=(i==1 ? "Temperature" : ""), xticklabelrotation=45.0)
        tmp = reshape(basis[1:M*N,i], (M, N))
        if n == 1
            heatmap!(ax, lonvec2, latvec,  tmp, colormap=:thermometer, colorrange=(minimum(tmp)*1.5, maximum(tmp)))
        else
            heatmap!(ax, lonvec2, latvec, tmp, colormap=:thermometer)
        end

        ex = Axis(fig[2,n], xlabel="Year", ylabel=(i==1 ? "Coefficient of mode (annual mean)" : ""), title="Mode $i")
        for x in 1:num_ens_members
            lines!(ex, time_history, (n==1 ? -1 : 1) .* month_to_year_avg(basis_trends_history[:,x,i,1]), color=scenario_colors["historical"], alpha=0.1, linestyle=:dot)
        end
        mean_trend = (n==1 ? -1 : 1) .* month_to_year_avg(mean(basis_trends_history[:,:,i,1], dims=2)) 
        lines!(ex, time_history, mean_trend, color=scenario_colors["historical"], alpha=1, linestyle=:solid, linewidth=2, label="historical")
        for (j, scenario) in enumerate(scenarios[2:end]) 
            for x in 1:num_ens_members
                lines!(ex, time_future, (n==1 ? -1 : 1) .* month_to_year_avg(basis_trends_future[:,x,i,j]), color=scenario_colors[scenario], alpha=0.1, linestyle=:dot)
            end
            mean_trend = (n==1 ? -1 : 1) .* month_to_year_avg(mean(basis_trends_future[:,:,i,j], dims=2))
            lines!(ex, time_future, mean_trend, color=scenario_colors[scenario], alpha=1, linestyle=:solid, linewidth=2, label=scenario)
        end
        if n==1
            axislegend(ex, position=:lt)
        end

        ax2 = GeoAxis(fig[3,n],  title="Mode $i", ylabel=(i==1 ? uppercasefirst(var_labels[second_var]) : ""), xticklabelrotation=45.0)
        tmp = reshape(basis[M*N+1:end,i], (M, N))
        heatmap!(ax2, lonvec2, latvec, tmp, colormap=:plum)
    end
    save("figs/$parent_folder/basis_modes.png", fig)
    display(fig)
end