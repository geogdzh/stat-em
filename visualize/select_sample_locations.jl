## let's identify the locations
function get_indices(lon1, lon2, lat1, lat2, ts)
    lons = ts(lon1, lon2, lat1, lat2).lonvec
    lats = ts(lon1, lon2, lat1, lat2).latvec
    return [findfirst(x -> x == lons[1], ts.lonvec), 
    findfirst(x -> x == lons[end], ts.lonvec), 
    findfirst(x -> x == lats[1], ts.latvec), 
    findfirst(x -> x == lats[end], ts.latvec)]
end

sampling_indices = zeros(Int64, (6, 4))
sampling_indices[1,:] = get_indices(272,273,41,42, ts3) #chicago
sampling_indices[2,:] = get_indices(300,301,-4,-3, ts3) #manaus, brazil
sampling_indices[3,:] = get_indices(7,9,11,12, ts3) #nigeria
sampling_indices[4,:] = get_indices(26,28,68,70, ts3) #finland
sampling_indices[5,:] = get_indices(90,92,23,25, ts3) #dhaka, bhangladesh
sampling_indices[6,:] = get_indices(174, 176, -38, -36, ts3) #auckland, nz

sampling_labels = ["Chicago, USA", "Manaus, Brazil", "Kano, Nigeria", "Inari, Finland", "Dhaka, Bhangladesh", "Auckland, NZ"]

#### start with: real data PDFs for last ten years of SSP585
for variable in ["tas", second_var]
    end_points = zeros(6, num_ens_members, 120)
    for n in 1:num_ens_members
        file = file_head*"ssp585/$(variable)/r$(n)i1p1f1_ssp585_$(variable).nc"
        ts = ncData(file, variable)
        data = apply_transform(ts.data, variable; hurs_option=hurs_option)
        for ind in 1:6
            subset = data[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], end-119:end]
            avg = mean(subset, dims=(1,2))[:]
            end_points[ind, n, :] = avg
        end
    end
    wfile = h5open("data/ground_truth/location_samples_$(variable)_ssp585_$(num_ens_members)ens.hdf5", "w")
    write(wfile, "end_points", end_points)
    write(wfile, "sampling_indices", sampling_indices)
    write(wfile, "sampling_labels", sampling_labels)

    # normal approximation to ground truth
    hfile = h5open("data/ground_truth/vars_$(variable)_ssp585_$(num_ens_members)ens.hdf5", "r") # true CMIP vars incl means
    true_var = read(hfile, "true_var")
    true_ens_mean = read(hfile, "true_ens_mean")
    close(hfile)
    sample_means = zeros(6, l2)
    sample_vars = zeros(6, l2)
    for ind in 1:6
        slice = true_ens_mean[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], 1:12:end] #only jan!!
        sample_means[ind, :] = mean(slice, dims=(1,2))[:]
        slice = true_var[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], 1:12:end] #only jan!!
        sample_vars[ind, :] = mean(slice, dims=(1,2))[:] 
    end
    write(wfile, "sample_means", sample_means)
    write(wfile, "sample_vars", sample_vars)
    close(wfile)
end

