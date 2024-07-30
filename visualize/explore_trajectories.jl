#load emulator
d = 100
hfile = h5open("data/$(parent_folder)/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
mean_coefs = read(hfile, "mean_coefs_1")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
if non_dim
    temp_factor = read(hfile, "temp_factor")
    if using_two
        pr_factor = read(hfile, "pr_factor")
    end
end
close(hfile)

#get gmt list
hfile = h5open("data/ssp245_gmts_$(num_ens_members)ens.hdf5")
ens_gmt = read(hfile, "ens_gmt")
close(hfile)

gmts = mean(ens_gmt, dims=1)[:]

ensemble_size = 100
trajectories = zeros(M, N, L2, ensemble_size)
trajectories_hurs = zeros(M, N, L2, ensemble_size)
for n in 1:ensemble_size
    traj = emulate(gmts, mean_coefs, chol_coefs) # it's not implemented to work with quadratic yet
    datatraj = back_to_data(traj, basis) 
    # trajectories[:,:,:,n] = shape_data(datatraj[1:M*N, :], M, N, true) #.* temp_factor
    trajectories_hurs[:,:,:,n] = shape_data(datatraj[M*N+1:end, :], M, N, true) #.* pr_factor
end

endtime = time_future[end-59:end]
starttime = time_future[1:60]

#load locations
hfile = h5open("data/ground_truth/location_samples_tas_ssp585_$(num_ens_members)ens.hdf5", "r")
sampling_indices = read(hfile, "sampling_indices")
sampling_labels = read(hfile, "sampling_labels")
close(hfile)

begin
    fig = Figure()
    for i in 1:3
        for j in 1:2
            ind = (i-1)*2 + j
            ax = Axis(fig[i,j], title=sampling_labels[ind])
            for n in 1:ensemble_size
                # shaped = trajectories[:,:,end-59:end,n]
                shaped = trajectories[:,:,1:60,n]
                subset = shaped[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], :]
                avg = mean(subset, dims=(1,2))[:]
                # lines!(ax, endtime, avg, color=:black, alpha=0.6)
                lines!(ax, starttime, avg, color=:black, alpha=0.6)
            end
        end
    end
    # save("data/$(parent_folder)/figs/$(variable)_start_tranjectories_ssp245.png", fig)
    display(fig)
end





#### dev
hfile = h5open("data/$(parent_folder)/training_data_ssp585_100d_28ens.hdf5", "r")
projts = read(hfile, "projts")