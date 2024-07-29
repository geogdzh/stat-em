function get_projts_history(d::Int)
    if isfile("data/$(parent_folder)/projts_historical_$(d)d_$(num_ens_members)ens.hdf5") && isfile("data/$(parent_folder)/projts_ssp585_$(d)d_$(num_ens_members)ens.hdf5")
        println("projts already exists!")
        flush(stdout)
        return nothing
    end

    ############## load a basis
    hfile = h5open("data/$(parent_folder)/basis_2000d.hdf5", "r") #this basis is calculated from just one ens member
    basis = read(hfile, "basis")
    if non_dim
        temp_factor = read(hfile, "temp_factor")
        if using_two
            two_factor = read(hfile, "two_factor")
        end
    end
    if use_metrics
        metric = read(hfile, "metric")
    end
    close(hfile)
    basis = basis[:, 1:d]

    ############### 
    for scenario in scenarios
        println("working on $(scenario)")
        flush(stdout)
        projts = scenario == "historical" ? zeros((d, (L1), num_ens_members)) : zeros((d, (L2), num_ens_members))
        ens_gmt = scenario == "historical" ? zeros((num_ens_members, Int((L1)/12))) : zeros((num_ens_members, Int((L2)/12)))

        file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
        for (i, n) in enumerate(ensemble_members)
            println("working on ensemble member $(n)")
            flush(stdout)
            files = file_head*"$(scenario)/tas/r$(n)i1p1f1_$(scenario)_tas.nc"
            files_two = file_head*"$(scenario)/$(second_var)/r$(n)i1p1f1_$(scenario)_$(second_var).nc"
            tmps = ncData(files, "tas")
            
            if using_two
                data1 = use_metrics ? sqrt.(metric) .* tmps.data : tmps.data
                data1 = non_dim ? reshape_data(data1) ./ temp_factor : reshape_data(data1)

                twos = ncData(files_two, second_var)
                twodata = apply_transform(twos.data, second_var)
                data2 = use_metrics ? sqrt.(metric) .* twodata : twodata
                data2 = non_dim ? reshape_data(data2) ./ two_factor : reshape_data(data2)
                data = vcat(data1, data2)
                projts[:, :, i] =  project_timeseries(data, basis, reshaped=true)
            else
                projts[:, :, i] =  non_dim ? project_timeseries(tmps.data, basis) ./ temp_factor : project_timeseries(tmps.data, basis)
                if use_metrics #this is pure laziness
                    throw("metrics not implemented for single variable")
                end
            end
            ens_gmt[i, :] = get_gmt_list(tmps) 
        end

        hfile = h5open("data/$(parent_folder)/projts_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w")
        write(hfile, "projts", projts)
        write(hfile, "ens_gmt", ens_gmt)
        if non_dim
            write(hfile, "temp_factor", temp_factor)
            if using_two
                write(hfile, "two_factor", two_factor)
            end
        end
        if use_metrics
            write(hfile, "metric", metric)
        end
        close(hfile)
    end
end