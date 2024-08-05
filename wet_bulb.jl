
function wet_bulb(T, RH)
    Tw = (T .* atan.(0.151977 .* sqrt.(RH .+ 8.313659)) 
        .+ atan.(T .+ RH) .- atan.(RH .- 1.676331) 
        .+ 0.00391838 .* (RH .^ 1.5) .* atan.(0.023101 .* RH) .- 4.686035)
    return Tw
end

for (i, n) in enumerate(ensemble_members)
    nc =  ncData(file_head*"historical/tas/r$(n)i1p1f1_historical_tas.nc", "tas")
    tas = nc.data .- 273.15
    nc =  ncData(file_head*"historical/hurs/r$(n)i1p1f1_historical_hurs.nc", "hurs")
    hurs = nc.data

    Tw = wet_bulb(tas, hurs)
    wfile = h5open("data/ground_truth/wet_bulb/historical/r$(n)i1p1f1_historical_wet_bulb.hdf5", "w")
    write(wfile, "wet_bulb", Tw)
    close(wfile)
end

#get the no log version of hurs
L = L1

all_gmts = zeros((num_ens_members, Int(L/12))) 
all_data = zeros((M,N, L, num_ens_members))
for (i, n) in enumerate(ensemble_members)
    println("working on ensemble member $(n)")
    flush(stdout)
    file = file_head*"historical/hurs/r$(n)i1p1f1_historical_hurs.nc"
    ts = ncData(file, "hurs") 
    all_gmts[i,:] = get_gmt_list(ts)
    all_data[:,:,:,i] = apply_transform(ts.data[:,:,:], "hurs"; hurs_option="nolog")
end

true_var = var(all_data, dims=4)[:,:,:,1]
true_ens_gmt = mean(all_gmts, dims=1)[:]
true_ens_mean = mean(all_data, dims=4)[:,:,:,1]

hfile = h5open("data/ground_truth/vars_hurs_historical_$(num_ens_members)ens_NOLOG.hdf5", "w")
write(hfile, "true_var", true_var)
write(hfile, "true_ens_gmt", true_ens_gmt)
write(hfile, "true_ens_mean", true_ens_mean)
close(hfile)


#now get it from the mean
hfile = h5open("data/ground_truth/vars_tas_historical_28ens.hdf5", "r")
tas_mean = read(hfile, "true_ens_mean")
tas_mean = tas_mean .- 273.15
close(hfile)
hfile = h5open("data/ground_truth/vars_hurs_historical_28ens_NOLOG.hdf5", "r")
hurs_mean = read(hfile, "true_ens_mean")
close(hfile)

Tw_mean = wet_bulb(tas_mean, hurs_mean)
wfile = h5open("data/ground_truth/wet_bulb/historical/means.hdf5", "cw")
write(wfile, "Tw_from_means", Tw_mean)
close(wfile)

#now get it from the ensemble
alldata = zeros((M,N, L, num_ens_members))
for (i, n) in enumerate(ensemble_members)
    hfile = h5open("data/ground_truth/wet_bulb/historical/r$(n)i1p1f1_historical_wet_bulb.hdf5", "r")
    Tw = read(hfile, "wet_bulb")
    close(hfile)
    alldata[:,:,:,i] = Tw
end
Tw_mean = mean(alldata, dims=4)[:,:,:,1]
wfile = h5open("data/ground_truth/wet_bulb/historical/means.hdf5", "cw")
write(wfile, "Tw_from_ensemble", Tw_mean)
close(wfile)

#### compare
hfile = h5open("data/ground_truth/wet_bulb/historical/means.hdf5", "r")
Tw_from_mean = read(hfile, "Tw_from_means")
Tw_from_ensemble = read(hfile, "Tw_from_ensemble")
close(hfile)

avg_means = month_to_year_avg(weighted_avg(Tw_from_mean, latvec))
avg_ensemble = month_to_year_avg(weighted_avg(Tw_from_ensemble, latvec))

begin
    fig = Figure()
    ax = Axis(fig[1,1], ylabel="Global avg Wet Bulb Temperature (C)", xlabel="Year", title="Wet Bulb Temperature from Historical Ensemble")
    lines!(ax, time_history, avg_means, color=:red, linestyle=:solid, label="from means")
    lines!(ax, time_history, avg_ensemble, color=:blue, linestyle=:solid, label="from ensemble")
    axislegend(ax, position=:lt)
    fig
end

# see it locally?
hfile = h5open("data/ground_truth/location_samples_hurs_ssp585_28ens.hdf5", "r")
sampling_labels = read(hfile, "sampling_labels")
sampling_indices = read(hfile, "sampling_indices")
close(hfile)

end_points_means = zeros(6, 120)
end_points_ensemble = zeros(6, 120)
hfile = h5open("data/ground_truth/wet_bulb/historical/means.hdf5", "r")
Tw_from_ensemble = read(hfile, "Tw_from_ensemble")
Tw_from_mean  = read(hfile, "Tw_from_means")
close(hfile)
for ind in 1:6
    subset = Tw_from_mean[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], end-119:end]
    avg = mean(subset, dims=(1,2))[:]
    end_points_means[ind, :] = avg
    subset = Tw_from_ensemble[sampling_indices[ind,1]:sampling_indices[ind,2], sampling_indices[ind,3]:sampling_indices[ind,4], end-119:end]
    avg = mean(subset, dims=(1,2))[:]
    end_points_ensemble[ind, :] = avg
end

begin
    data = end_points[:,:,1:12:end]# looking only at January!
    fig = Figure(resolution=(1000,750)) #title="Specific Humidity in January 2095-2100"
    for i in 1:2
        for j in 1:3
            ind = (i-1)*3+j
            ax = Axis(fig[i,j], title=sampling_labels[ind],  ylabel=(j==1 ? "Frequency" : ""))
            samples = end_points_means[ind, :]
            hist!(ax, samples, bins=15,  normalization=:pdf, color=(:red,0.5), label="from mean")
            samples = end_points_ensemble[ind, :]
            hist!(ax, samples, bins=15,  normalization=:pdf, color=(:blue,0.5), label="from ensemble")
            if ind == 3
                axislegend(ax, position=:lt)
            end
        end
    end
    fig
end