function visualize_trajectories(variable)
    label = variable == "tas" ? variable : "two"
    d = 100

    hfile = h5open("data/ground_truth/location_samples_$(variable)_ssp585_$(num_ens_members)ens.hdf5", "r")
    sampling_indices = read(hfile, "sampling_indices")
    sampling_labels = read(hfile, "sampling_labels")
    end_points = read(hfile, "end_points")
    close(hfile)

    hfile = h5open("data/ssp585_gmts_$(num_ens_members)ens.hdf5", "r")
    ens_gmt = read(hfile, "ens_gmt")
    close(hfile)
    ens_gmt = mean(ens_gmt, dims=1)[:]
    sample_gmts = ens_gmt[end-9:end]

    #load the emulator
    hfile = h5open("data/$parent_folder/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
    mean_coefs = read(hfile, "mean_coefs_2")
    chol_coefs = read(hfile, "chol_coefs")
    basis = read(hfile, "basis")
    close(hfile)

    # sample a trajectory without transitions
    ensemble_size = 10
    sample_trajs = zeros(6, ensemble_size, 120) # hardcoded for 10 years
    for n in 1:ensemble_size
        traj = emulate_no_cov(sample_gmts, mean_coefs, chol_coefs)
        data = back_to_data(traj, basis)
        data2 = variable == "tas" ? shape_data(data[1:M*N,:], M, N, true) : shape_data(data[M*N+1:end,:], M, N, true)
        for i in 1:6
            sample_trajs[i, n, :] = mean(data2[sampling_indices[i,1]:sampling_indices[i,2], sampling_indices[i,3]:sampling_indices[i,4], :], dims=(1,2))[:]
        end
    end


    begin
        fig = Figure(resolution=(1500, 750))
        for i in 1:2
            for j in 1:3
                ind = (i-1)*3+j
                ax = Axis(fig[i,j], title=sampling_labels[ind], xlabel=(i == 2 ? "Time (months)" : ""), ylabel=(j==1 ? "$(uppercasefirst(var_labels[variable])) ($(unit_labels[variable]))" : ""))
                for n in 1:ensemble_size
                    lines!(ax, 1:120, end_points[ind, n, :], color=:red, alpha=0.3, label="MPI")
                    lines!(ax, 1:120, sample_trajs[ind, n, :], color=:blue, alpha=0.3, label="emulator")
                end
                if j == 1 && i == 1
                    axislegend(ax, position=:lb, merge=true)
                end
            end
        end
        save("figs/$parent_folder/$(variable)_sample_trajectories.png", fig)
        display(fig)
    end
end