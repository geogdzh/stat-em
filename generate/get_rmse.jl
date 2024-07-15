using HDF5
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl") 

#################### ok let's test it out for real

#get a sample gmt list and the latvec to be used later on
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
file3 = file_head*"ssp585/tas/r1i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas")
lonvec, latvec = ts3.lonvec[:], ts3.latvec[:]
lonvec2 = lonvec .-180.

using_precip = true 
non_dim = true  
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


## generate statistics
scenarios = ["historical", "ssp585", "ssp245", "ssp119"]
time_history = [x for x in 1850:2014]
time_future = [x for x in 2015:2100]
L1, L2 = 1980, 1032 
l1, l2 = 165, 86


function calculate_rmse(numbers, variable, scenarios; rel_error=false, for_k=false)
    # println("working on $(variable) and d = $(numbers[1]) and for_k is $(for_k)")
    # flush(stdout) 
    if non_dim #is this gonna be necessary? maybe not
        hfile = h5open("data/$parent_folder/basis_1000d.hdf5", "r") 
        factor = read(hfile, "$(variable)_factor")
        close(hfile)
    end

    variable = variable == "temp" ? "tas" : "pr"

    for scenario in scenarios[2:end]
        # once we have the emulator and have saved out the emulator ensemble vars/means etc; calculate the rmse compared to various scenarios

        hfile = h5open("data/ground_truth/vars_$(variable)_$(scenario)_50ens.hdf5", "r") # true CMIP vars
        true_var = read(hfile, "true_var")
        true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1]
        close(hfile)

        wfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$(scenario).hdf5", "cw") #let's make this also include means
        
        if for_k
            hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_k.hdf5", "r") #emulator vars - this includes means

            for number in numbers
                ens_means = read(hfile, "ens_means_$(variable)_k$(number)")

                rmse_means = sqrt.(sum((true_ens_mean.-ens_means).^2, dims=3)[:,:,1]./size(true_ens_mean)[3]) #time average rmse (shaped as a map)
                rmse_means_time = sqrt.(weighted_avg((true_ens_mean.-ens_means).^2, latvec)) #spatial average rmse (shaped as a timeseries
        
                write(wfile, "rmse_means_$(variable)_k$(number)", rmse_means)
                write(wfile, "rmse_means_time_$(variable)_k$(number)", rmse_means_time)
            end

            close(hfile)
        else
            for number in numbers
                hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_$(number)d.hdf5", "r") #emulator vars - this includes means

                ens_means =  read(hfile, "ens_means_$(variable)_$(number)")
                ens_vars = read(hfile, "ens_vars_$(variable)_$(number)")

                # if non_dim #DOUBLE CHECK!
                #     ens_means = ens_means .* factor
                #     ens_vars = ens_vars .* factor.^2
                # end

                rmse_stds = sqrt.(sum((sqrt.(true_var) .- sqrt.(ens_vars)).^2, dims=3)[:,:,1]./size(true_var)[3]) #time average rmse (shaped as a map)
                rmse_stds_time = sqrt.(weighted_avg((sqrt.(true_var) .- sqrt.(ens_vars)).^2, latvec)) #spatial average rmse (shaped as a timeseries)
                rmse_means = sqrt.(sum((true_ens_mean.-ens_means).^2, dims=3)[:,:,1]./size(true_ens_mean)[3]) #time average rmse (shaped as a map)
                rmse_means_time = sqrt.(weighted_avg((true_ens_mean.-ens_means).^2, latvec)) #spatial average rmse (shaped as a timeseries

                write(wfile, "rmse_stds_$(variable)_$(number)", rmse_stds)
                write(wfile, "rmse_stds_time_$(variable)_$(number)", rmse_stds_time)
                write(wfile, "rmse_means_$(variable)_$(number)", rmse_means)
                write(wfile, "rmse_means_time_$(variable)_$(number)", rmse_means_time)

                close(hfile)
            end
        end
        
        close(wfile)
    end
end

for variable in ["temp", "pr"]
    rel_error = false

    numbers = [10, 100]
    calculate_rmse(numbers, variable, scenarios; rel_error=rel_error, for_k=false)
    
    ks = [x for x in 1:2]
    calculate_rmse(ks, variable, scenarios; rel_error=rel_error, for_k=true)
end
