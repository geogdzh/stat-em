function generate_pattern_scaling(variable, measure)
    hfile = h5open("data/ground_truth/vars_$(variable)_historical_28ens.hdf5", "r")
    true_ens_mean_hist = read(hfile, "true_ens_mean")
    ens_gmt_hist = read(hfile, "true_ens_gmt")
    close(hfile)

    hfile = h5open("data/ground_truth/vars_$(variable)_ssp585_28ens.hdf5", "r")
    true_ens_mean_ssp585 = read(hfile, "true_ens_mean")
    ens_gmt_ssp585 = read(hfile, "true_ens_gmt")
    close(hfile)

    GMT1850 = ens_gmt_hist[1]
    GMT2100 = ens_gmt_ssp585[end]

    function get_w(GMT)
        return (GMT - GMT1850)/(GMT2100 - GMT1850)
    end

    function pattern_scaling(GMT, month, true_ens_hist, true_ens_ssp585)
        # month as Int
        # Tbar doesn't actually have to be T
        w = get_w(GMT)
        Tbar_1850  = true_ens_hist[:,:,month]
        Tbar_2100 = true_ens_ssp585[:,:,end-12+month]
        return (1-w)*Tbar_1850 + w*Tbar_2100
    end

    ### let's generate it for all the GMTs
    for scenario in scenarios[2:end]
        hfile = h5open("data/$(scenario)_gmts_$(num_ens_members)ens.hdf5", "r") 
        ens_gmt = mean(read(hfile, "ens_gmt"), dims=1)[:]
        close(hfile)

        hfile = h5open("data/ground_truth/vars_$(variable)_historical_28ens.hdf5", "r")
        true_hist = measure == "mean" ? read(hfile, "true_ens_mean") : sqrt.(read(hfile, "true_var"))
        ens_gmt_hist = read(hfile, "true_ens_gmt")
        close(hfile)

        hfile = h5open("data/ground_truth/vars_$(variable)_ssp585_28ens.hdf5", "r")
        true_ssp585 = measure == "mean" ? read(hfile, "true_ens_mean") : sqrt.(read(hfile, "true_var"))
        ens_gmt_ssp585 = read(hfile, "true_ens_gmt")
        close(hfile)

        data = zeros(M, N, L2)
        for (i, T) in enumerate(ens_gmt)
            for month in 1:12
                data[:,:,12*(i-1)+month] = pattern_scaling(T, month, true_hist, true_ssp585) #need to double check indexing here
            end
        end

        wfile = h5open("data/pattern_scaling/projected_$(scenario)_$(variable)_$(num_ens_members)ens.hdf5", "cw")
        write(wfile, "$measure", data)
        close(wfile)
    end

    # # and get its rmse
    for scenario in scenarios[2:end]
        hfile = h5open("data/pattern_scaling/projected_$(scenario)_$(variable)_$(num_ens_members)ens.hdf5", "r")
        data = read(hfile, "$measure")
        close(hfile)

        hfile = h5open("data/ground_truth/vars_$(variable)_$(scenario)_$(num_ens_members)ens.hdf5", "r")
        true_data = measure == "mean" ? read(hfile, "true_ens_mean") : sqrt.(read(hfile, "true_var"))
        close(hfile)

        rmse = sqrt.(sum((true_data.-data).^2, dims=3)[:,:,1]./size(true_data)[3]) #time average rmse (shaped as a map)
        rmse_time = sqrt.(weighted_avg((true_data.-data).^2, latvec)) #spatial average rmse (shaped as a timeseries)

        wfile = h5open("data/pattern_scaling/rmse_$(scenario)_$(variable)_$(num_ens_members)ens.hdf5", "cw")
        write(wfile, "rmse_$measure", rmse)
        write(wfile, "rmse_time_$measure", rmse_time)
        close(wfile)
    end
end