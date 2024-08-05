function calculate_rmse(numbers, variable, scenarios; for_k=false, which_k=1)
 
    label = variable == "tas" ? variable : "two"

    for scenario in scenarios[2:end]
        # once we have the emulator and have saved out the emulator ensemble vars/means etc; calculate the rmse compared to various scenarios

        hfile = h5open("data/ground_truth/vars_$(variable)_$(scenario)_$(num_ens_members)ens.hdf5", "r") # true CMIP vars
        true_var = read(hfile, "true_var")
        true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1]
        close(hfile)

        wfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$(scenario).hdf5", "cw") #let's make this also include means
        
        if for_k
            hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_k.hdf5", "r") #emulator vars - this includes means
  
            for number in numbers
                ens_means = read(hfile, "ens_means_$(label)_k$(number)")

                rmse_means = sqrt.(sum((true_ens_mean.-ens_means).^2, dims=3)[:,:,1]./size(true_ens_mean)[3]) #time average rmse (shaped as a map)
                rmse_means_time = sqrt.(weighted_avg((true_ens_mean.-ens_means).^2, latvec)) #spatial average rmse (shaped as a timeseries
        
                write(wfile, "rmse_means_$(label)_k$(number)", rmse_means)
                write(wfile, "rmse_means_time_$(label)_k$(number)", rmse_means_time)
            end

            close(hfile)
        else
            for number in numbers
                hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_$(number)d.hdf5", "r") #emulator vars - this includes means

                ens_means =  read(hfile, "ens_means_$(label)_$(number)_k$(which_k)")
                ens_vars = read(hfile, "ens_vars_$(label)_$(number)")

                rmse_stds = sqrt.(sum((sqrt.(true_var) .- sqrt.(ens_vars)).^2, dims=3)[:,:,1]./size(true_var)[3]) #time average rmse (shaped as a map)
                rmse_stds_time = sqrt.(weighted_avg((sqrt.(true_var) .- sqrt.(ens_vars)).^2, latvec)) #spatial average rmse (shaped as a timeseries)
                rmse_means = sqrt.(sum((true_ens_mean.-ens_means).^2, dims=3)[:,:,1]./size(true_ens_mean)[3]) #time average rmse (shaped as a map)
                rmse_means_time = sqrt.(weighted_avg((true_ens_mean.-ens_means).^2, latvec)) #spatial average rmse (shaped as a timeseries

                write(wfile, "rmse_stds_$(label)_$(number)", rmse_stds)
                write(wfile, "rmse_stds_time_$(label)_$(number)", rmse_stds_time)
                write(wfile, "rmse_means_$(label)_$(number)", rmse_means)
                write(wfile, "rmse_means_time_$(label)_$(number)", rmse_means_time)

                close(hfile)
            end
        end
        
        close(wfile)
    end
end