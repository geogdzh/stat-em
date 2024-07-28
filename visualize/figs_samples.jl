

#####
label = variable == "tas" ? variable : "two"

# get the ground truth
hfile = h5open("data/ground_truth/location_samples_$(variable)_ssp585_$(num_ens_members)ens.hdf5", "r")
sampling_indices = read(hfile, "sampling_indices")
sampling_labels = read(hfile, "sampling_labels")
end_points = read(hfile, "end_points")
sample_means = read(hfile, "sample_means")
sample_vars = read(hfile, "sample_vars")
close(hfile)

## get emulator approximation
hfile = h5open("data/$parent_folder/ens_vars/ens_vars_ssp585_100d.hdf5", "r")
ens_means = read(hfile, "ens_means_$(label)_100")
ens_vars = read(hfile, "ens_vars_$(label)_100")
close(hfile)

sample_ens_means = zeros(6, l2)
sample_ens_vars = zeros(6, l2)
for ind in 1:6
    slice = ens_means[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], 1:12:end] #only jan
    sample_ens_means[ind, :] = mean(slice, dims=(1,2))[:] 
    slice = ens_vars[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], 1:12:end] #only jan
    sample_ens_vars[ind, :] = mean(slice, dims=(1,2))[:]  
end

##
begin
    data = end_points[:,:,1:12:end]# looking only at January!
    fig = Figure(resolution=(1000,750)) #title="Specific Humidity in January 2095-2100"
    for i in 1:2
        for j in 1:3
            ind = (i-1)*3+j
            ax = Axis(fig[i,j], title=sampling_labels[ind]) #xlabel="Temperature (K)", ylabel="Frequency",
            samples = reshape(data[ind,:,:], (num_ens_members*Int(size(end_points)[3]/12),))
            hist!(ax, samples, bins=15,  normalization=:pdf)
# 
            # add the distributions
            params = [(sample_means[ind, x], sqrt(sample_vars[ind, x])) for x in 76:86]
            dist = MixtureModel(Normal, params)
            plot!(ax, dist, label="True distribution", color=:red)

            # params2 = [(sample_ens_means[ind, x], sqrt(sample_ens_vars[ind, x])) for x in 76:86]
            # dist2 = MixtureModel(Normal, params2)
            # plot!(ax, dist2, label="Emulator distribution", color=:blue)
        end
    end
    # save("figs/$parent_folder/$(variable)_samples.png", fig)
    display(fig)
end

