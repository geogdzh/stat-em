function get_ens_vars(d, true_ens_gmt; get_means=false, k=2) # OR the means lol
    M = 192
    N = 96
    L = 1032
    hfile = h5open("data/$(parent_folder)/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
    mean_coefs = read(hfile, "mean_coefs_$(k)")
    chol_coefs = read(hfile, "chol_coefs")
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
    ens_vars_tas = zeros(M, N, L)
    ens_vars_two = zeros(M, N, L)
    for m in ProgressBar(1:Int(L/12))
        println("working on year $(m)")
        flush(stdout)
        if get_means
            for n in 1:12
                if k==1 
                    data = back_to_data(mean_coefs[n,:,2].*true_ens_gmt[m] .+ mean_coefs[n, :, 1], basis)
                elseif k==2
                    data = back_to_data(mean_coefs[n,:,2].*true_ens_gmt[m] .+ mean_coefs[n, :, 1] .+ mean_coefs[n,:,3].*true_ens_gmt[m].^2, basis)
                end
                newdata = non_dim ? shape_data(data[1:M*N,:], M, N) .* temp_factor : shape_data(data[1:M*N,:], M, N)
                newdata = use_metrics ? newdata ./ sqrt.(metric) : newdata
                ens_vars_tas[:,:,(m-1)*12+n] = newdata
                if using_two
                    new2data = non_dim ? shape_data(data[M*N+1:end,:], M, N) .* two_factor : shape_data(data[M*N+1:end,:], M, N)
                    new2data = use_metrics ? new2data ./ sqrt.(metric) : new2data
                    ens_vars_two[:,:,(m-1)*12+n] = new2data
                end
            end
        else
            co = get_cov(true_ens_gmt[m], chol_coefs) 
            for n in 1:12
                data = sum([co[:,:,n][i,j].*basis[:,i].*basis[:,j] for i in 1:d, j in 1:d]) 
                newdata = non_dim ? shape_data(data[1:M*N,:], M, N) .* temp_factor^2 : shape_data(data[1:M*N,:], M, N)
                newdata =  use_metrics ? newdata ./ sqrt.(metric) : newdata
                ens_vars_tas[:,:,(m-1)*12+n] = newdata
                if using_two
                    new2data = non_dim ? shape_data(data[M*N+1:end,:], M, N) .* two_factor^2 : shape_data(data[M*N+1:end,:], M, N)
                    new2data = use_metrics ? new2data ./ sqrt.(metric) : new2data
                    ens_vars_two[:,:,(m-1)*12+n] = new2data
                end
            end
        end
    end
    println("done")
    return using_two ? (ens_vars_tas, ens_vars_two) : ens_vars_tas
end


function run_ens_vars(param, d)
    if param == "d"
        if isfile("data/$(parent_folder)/ens_vars/ens_vars_historical_$(d)d.hdf5") #generalize!
            return nothing
        end
        # test the differnt values of n
        for scenario in scenarios[2:end] 
            println("working on $(scenario)")
            flush(stdout)
            #get the true gmt
            hfile = h5open("data/$(scenario)_gmts_$(num_ens_members)ens.hdf5", "r") 
            ens_gmt = read(hfile, "ens_gmt")
            true_ens_gmt = mean(ens_gmt, dims=1)[:]
            close(hfile)

            println("working on $(d)")
            flush(stdout)
            if using_two
                ens_vars_tas, ens_vars_two = get_ens_vars(d, true_ens_gmt)
                ens_means_tas, ens_means_two = get_ens_vars(d, true_ens_gmt; get_means=true)
                hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_$(d)d.hdf5", "w")
                write(hfile, "ens_vars_tas_$(d)", ens_vars_tas)
                write(hfile, "ens_vars_two_$(d)", ens_vars_two)
                write(hfile, "ens_means_tas_$(d)", ens_means_tas)
                write(hfile, "ens_means_two_$(d)", ens_means_two)
                close(hfile)
            else
                ens_vars_tas = get_ens_vars(d, true_ens_gmt)
                ens_means_tas = get_ens_vars(d, true_ens_gmt; get_means=true)
                hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_$(d)d.hdf5", "w")
                write(hfile, "ens_vars_tas_$(d)", ens_vars_tas)
                write(hfile, "ens_means_tas_$(d)", ens_means_tas)
                close(hfile)
            end
        end

    elseif param == "k"
        if isfile("data/$(parent_folder)/ens_vars/ens_vars_historical_k.hdf5") #generalize!
            return nothing
        end
        # # test the differnt values of k
        for scenario in scenarios[2:end] 
            println("working on $(scenario)")
            flush(stdout)
            hfile = h5open("data/$(scenario)_gmts_$(num_ens_members)ens.hdf5", "r") 
            ens_gmt = read(hfile, "ens_gmt")
            true_ens_gmt = mean(ens_gmt, dims=1)[:]
            close(hfile)

            hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_k.hdf5", "w")
            for k in 1:2
                ens_means_tas, ens_means_two = get_ens_vars(d, true_ens_gmt; get_means=true, k=k)
                write(hfile, "ens_means_tas_k$(k)", ens_means_tas)
                write(hfile, "ens_means_two_k$(k)", ens_means_two)
            end
            write(hfile, "d", d)
            close(hfile)
        end
    end
end