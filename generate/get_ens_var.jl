using ProgressBars, HDF5
include("../utils/data_util.jl")
include("../utils/eof_util.jl")
include("../utils/emulator_util.jl") 

scenarios = ["historical", "ssp585", "ssp245", "ssp119"]

param = ARGS[1]
d = parse(Int, ARGS[2])
# if param == "d"
#     d = parse(Int, ARGS[2])
# else #if testing k 
#     d = 10 #CHANGE DEFAULT
# end

using_two = (ARGS[3] == "true" )
second_var = ARGS[4] # "pr" or "huss"
non_dim = (ARGS[5] == "true" )  
use_metrics = (ARGS[6] == "true" )
if using_two
    if second_var == "pr"
        parent_folder = "temp_precip"
    else
        parent_folder = "temp_huss"
    end
else
    parent_folder = "temp"
end
if non_dim
    parent_folder = "nondim"
end
if use_metrics && using_two
    parent_folder = "metrics"
elseif use_metrics && !using_two
    parent_folder = "temp_metrics"
end


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
            pr_factor = read(hfile, "pr_factor")
        end
    end
    if use_metrics
        metric = read(hfile, "metric") 
    end
    close(hfile)
    ens_vars_tas = zeros(M, N, L)
    ens_vars_pr = zeros(M, N, L)
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
                    new2data = non_dim ? shape_data(data[M*N+1:end,:], M, N) .* pr_factor : shape_data(data[M*N+1:end,:], M, N)
                    new2data = use_metrics ? new2data ./ sqrt.(metric) : new2data
                    ens_vars_pr[:,:,(m-1)*12+n] = new2data
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
                    new2data = non_dim ? shape_data(data[M*N+1:end,:], M, N) .* pr_factor^2 : shape_data(data[M*N+1:end,:], M, N)
                    new2data = use_metrics ? new2data ./ sqrt.(metric) : new2data
                    ens_vars_pr[:,:,(m-1)*12+n] = new2data
                end
            end
        end
    end
    println("done")
    return using_two ? (ens_vars_tas, ens_vars_pr) : ens_vars_tas
end

if param == "d"
    # test the differnt values of n
    for scenario in scenarios[2:end] 
        println("working on $(scenario)")
        flush(stdout)
        #get the true gmt
        hfile = h5open("data/$(scenario)_gmts_50ens.hdf5", "r") 
        ens_gmt = read(hfile, "ens_gmt")
        true_ens_gmt = mean(ens_gmt, dims=1)[:]
        close(hfile)

        
        println("working on $(d)")
        flush(stdout)
        if using_two
            ens_vars_tas, ens_vars_pr = get_ens_vars(d, true_ens_gmt)
            ens_means_tas, ens_means_pr = get_ens_vars(d, true_ens_gmt; get_means=true)
            hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_$(d)d.hdf5", "w")
            write(hfile, "ens_vars_tas_$(d)", ens_vars_tas)
            write(hfile, "ens_vars_pr_$(d)", ens_vars_pr)
            write(hfile, "ens_means_tas_$(d)", ens_means_tas)
            write(hfile, "ens_means_pr_$(d)", ens_means_pr)
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
    # # test the differnt values of k
    for scenario in scenarios[2:end] 
        println("working on $(scenario)")
        flush(stdout)
        hfile = h5open("data/$(scenario)_gmts_50ens.hdf5", "r") 
        ens_gmt = read(hfile, "ens_gmt")
        true_ens_gmt = mean(ens_gmt, dims=1)[:]
        close(hfile)

        hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_k.hdf5", "w")
        for k in 1:2
            ens_means_tas, ens_means_pr = get_ens_vars(d, true_ens_gmt; get_means=true, k=k)
            write(hfile, "ens_means_tas_k$(k)", ens_means_tas)
            write(hfile, "ens_means_pr_k$(k)", ens_means_pr)
        end
        write(hfile, "d", d)
        close(hfile)
    end
end