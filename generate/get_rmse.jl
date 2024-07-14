using CairoMakie, ProgressBars, HDF5, GeoMakie
include("utils.jl")
include("eof_util.jl")
include("emulator_util.jl") 

#################### ok let's test it out for real

#get a sample gmt list and the latvec to be used later on
# file_head = "/Users/masha/urop_2022/cmip/CMIP6/interim/"
file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
file3 = file_head*"ssp585/tas/r1i1p1f1_ssp585_tas.nc"
ts3 = ncData(file3, "tas")
lonvec, latvec = ts3.lonvec[:], ts3.latvec[:]
lonvec2 = lonvec .-180.

using_precip = true 
non_dim = false  
use_metrics = true
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
    if non_dim
        hfile = h5open("data/$parent_folder/temp_precip_basis_1000d.hdf5", "r") 
        factor = read(hfile, "$(variable)_factor")
        close(hfile)
    end

    for scenario in scenarios[2:end]
        # once we have the emulator and have saved out the emulator ensemble vars/means etc; calculate the rmse compared to various scenarios

        hfile = h5open("data/temp_precip/vars_$(variable)_$(scenario)_50ens.hdf5", "r") # true CMIP vars
        true_var = read(hfile, "true_var")
        true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1]
        close(hfile)

        wfile = h5open("data/$(parent_folder)/ens_vars_withpr_rmse_$(scenario).hdf5", "cw") #let's make this also include means
        hfile = h5open("data/$(parent_folder)/ens_vars_$(scenario).hdf5", "r") #emulator vars - this includes means
        for number in numbers
            ens_means = variable == "temp" ? read(hfile, "ens_means_tas_$(for_k==true ? 100 : number)") : read(hfile, "ens_means_$(variable)_$(for_k==true ? 100 : number)")
            ens_vars = variable == "temp" ? read(hfile, "ens_vars_tas_$(for_k==true ? 100 : number)") : read(hfile, "ens_vars_$(variable)_$(for_k==true ? 100 : number)")

            if non_dim #DOUBLE CHECK!
                ens_means = ens_means .* factor
                ens_vars = ens_vars .* factor.^2
            end

            if rel_error #NOT updated
                rmse_stds = sum(abs.(sqrt.(true_var).-sqrt.(ens_vars)) ./ mean((true_var), dims=3), dims=3) ./size(true_var)[3] #relative to variance
                rmse_stds_time = weighted_avg(abs.(sqrt.(true_var).-sqrt.(ens_vars)) ./ mean((true_var), dims=3), latvec)
                # error relative to std for the means
                rmse_means = sum(abs.(true_ens_mean.-ens_means) ./ mean(sqrt.(true_var), dims=3), dims=3) ./size(true_ens_mean)[3]
                rmse_means_time = weighted_avg(abs.(true_ens_mean.-ens_means) ./ mean(sqrt.(true_var), dims=3), latvec)

                write(wfile, "rmse_stds_$(variable)_$(number)_rel", rmse_stds)
                write(wfile, "rmse_stds_time_$(variable)_$(number)_rel", rmse_stds_time)
                write(wfile, "rmse_means_$(variable)_$(number)_rel", rmse_means)
                write(wfile, "rmse_means_time_$(variable)_$(number)_rel", rmse_means_time)
            elseif for_k
                ens_means = variable == "temp" ? read(hfile, "ens_means_tas_k$(number)") : read(hfile, "ens_means_$(variable)_k$(number)")

                rmse_means = sqrt.(sum((true_ens_mean.-ens_means).^2, dims=3)[:,:,1]./size(true_ens_mean)[3]) #time average rmse (shaped as a map)
                rmse_means_time = sqrt.(weighted_avg((true_ens_mean.-ens_means).^2, latvec)) #spatial average rmse (shaped as a timeseries
        
                write(wfile, "rmse_means_$(variable)_k$(number)", rmse_means)
                write(wfile, "rmse_means_time_$(variable)_k$(number)", rmse_means_time)

            else
                rmse_stds = sqrt.(sum((sqrt.(true_var) .- sqrt.(ens_vars)).^2, dims=3)[:,:,1]./size(true_var)[3]) #time average rmse (shaped as a map)
                rmse_stds_time = sqrt.(weighted_avg((sqrt.(true_var) .- sqrt.(ens_vars)).^2, latvec))#; already_weighted=use_metrics)) #spatial average rmse (shaped as a timeseries)
                rmse_means = sqrt.(sum((true_ens_mean.-ens_means).^2, dims=3)[:,:,1]./size(true_ens_mean)[3]) #time average rmse (shaped as a map)
                rmse_means_time = sqrt.(weighted_avg((true_ens_mean.-ens_means).^2, latvec))#; already_weighted=use_metrics)) #spatial average rmse (shaped as a timeseries

                write(wfile, "rmse_stds_$(variable)_$(number)", rmse_stds)
                write(wfile, "rmse_stds_time_$(variable)_$(number)", rmse_stds_time)
                write(wfile, "rmse_means_$(variable)_$(number)", rmse_means)
                write(wfile, "rmse_means_time_$(variable)_$(number)", rmse_means_time)
            end

        end
        close(hfile)
        close(wfile)
    end
end

for variable in ["temp", "pr"]
    rel_error = false

    numbers = [10]
    calculate_rmse(numbers, variable, scenarios; rel_error=rel_error, for_k=false)
    
    ks = [x for x in 1:2]
    calculate_rmse(ks, variable, scenarios; rel_error=rel_error, for_k=true)
end
