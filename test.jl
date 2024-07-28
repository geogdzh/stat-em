using CairoMakie, ProgressBars, HDF5, GeoMakie
include("utils/data_util.jl")
include("utils/eof_util.jl")
include("utils/emulator_util.jl") 

parent_folder = "huss_only"

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"

file_tail = "historical/huss/r10i1p1f1_historical_huss.nc"
ts = ncData(file_head*file_tail, "huss")
M, N, L1 = size(ts.data)

X = reshape_data(ts.data)

file_tail = "ssp585/huss/r10i1p1f1_ssp585_huss.nc"
ts2 = ncData(file_head*file_tail, "huss")
X2 = reshape_data(ts2.data)
M, N, L2 = size(ts2.data)

fullX = hcat(X, X2)

U, S, V = svd(fullX)
d = 2000
basis = U[:, 1:d]

hfile = h5open("data/$(parent_folder)/basis_$(d)d.hdf5", "w") #but to use a smaller basis, can just take fewer modes
write(hfile, "basis", basis)
close(hfile)

scenarios = ["historical", "ssp585"]

d = 20
basis = basis[:, 1:d]

ss = sum(S.^2)
cum_var = cumsum(S.^2 ./ ss)

begin
    fig = Figure(resolution=(600, 400))
    ax = Axis(fig[1,1], xlabel="Mode", ylabel="Cumulative variance explained", title="Cumulative variance explained by the modes")
    lines!(ax, cum_var, color=:black, linestyle=:solid)
    display(fig)
    # save("figs/cum_var.png", fig)
end

for scenario in scenarios
    println("working on $(scenario)")
    flush(stdout)
    num_ens_members = 50 # number of model runs used to train the emulator
    projts = scenario == "historical" ? zeros((d, (L1), num_ens_members)) : zeros((d, (L2), num_ens_members))
    ens_gmt = scenario == "historical" ? zeros((num_ens_members, Int((L1)/12))) : zeros((num_ens_members, Int((L2)/12)))

    file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"
    errors = []
    for i in 1:num_ens_members
        try
            println("working on ensemble member $(i)")
            flush(stdout)
            files = file_head*"$(scenario)/huss/r$(i)i1p1f1_$(scenario)_huss.nc"
            tmps = ncData(files, "huss")
            projts[:, :, i] = project_timeseries(tmps.data, basis)
            ens_gmt[i, :] = get_gmt_list(tmps) 
        catch
            println("missing values in ensemble member $(i)")
            push!(errors, i)
            flush(stdout)
        end
    end

    if length(errors) > 0
        projts = projts[:,:,1:end .!= errors[:]] # HARDCODED FOR ONLY ONE ERROR - need to generalize. issue was scoping in script run # adn for some reason error when no errors?
        ens_gmt = ens_gmt[1:end .!= errors[:], :]
    end
    num_ens_members = size(ens_gmt)[1]

    hfile = h5open("data/$(parent_folder)/projts_$(scenario)_$(d)d_$(num_ens_members)ens.hdf5", "w")
    write(hfile, "projts", projts)
    write(hfile, "ens_gmt", ens_gmt)
    write(hfile, "num_ens_members", num_ens_members)
    close(hfile)
end

### assemble training data

hfile = h5open("data/$(parent_folder)/projts_historical_$(d)d_50ens.hdf5", "r")
projts_history = read(hfile, "projts")
ens_gmt_history = read(hfile, "ens_gmt")
close(hfile)

hfile = h5open("data/$(parent_folder)/projts_$("ssp585")_$(d)d_50ens.hdf5", "r")
projts_ssp585 = read(hfile, "projts")
ens_gmt_ssp585 = read(hfile, "ens_gmt")
close(hfile)

num_ens_members = 50 # number of model runs used to train the emulator
projts = hcat(projts_history, projts_ssp585)
ens_gmt = hcat(ens_gmt_history, ens_gmt_ssp585)


hfile = h5open("data/$(parent_folder)/training_data_$("ssp585")_$(d)d_$(num_ens_members)ens.hdf5", "w")
write(hfile, "projts", projts)
write(hfile, "ens_gmt", ens_gmt)
write(hfile, "num_ens_members", num_ens_members)
close(hfile)

scenario  = "ssp585"
ens_projts = projts
#### get the emulator itself
ens_gmt = mean(ens_gmt, dims=1)
mean_coefs_1 = get_mean_coefs(ens_projts, ens_gmt, degree=1)
mean_coefs_2 = get_mean_coefs(ens_projts, ens_gmt, degree=2)
gmt_cov(ens_projts, ens_gmt, "$(scenario)_$(d)d") #saves out covs
chol_coefs = get_chol_coefs(ens_gmt, "$(scenario)_$(d)d")

hfile = h5open("data/$(parent_folder)/gaussian_emulator_$(scenario)_$(d)d.hdf5", "w")
write(hfile, "mean_coefs_1", mean_coefs_1)
write(hfile, "mean_coefs_2", mean_coefs_2)
write(hfile, "chol_coefs", chol_coefs)
write(hfile, "basis", basis)
write(hfile, "num_ens_members", num_ens_members)
write(hfile, "ens_gmt", ens_gmt)
close(hfile)



###################

basis_shape = shape_data(basis, M, N, true)

begin
    fig = Figure(resolution = (2500, 1000))
    for i in 1:2
        for j in 1:5
            ind = (i-1)*5+j
            ax = Axis(fig[i,j], title="Mode $(ind)")
            heatmap!(ax, basis_shape[:,:,ind])
        end
    end
    # save("figs/$(parent_folder)/basis_$(d)d.png", fig)
    display(fig)
end

hfile = h5open("data/temp_huss/basis_2000d.hdf5", "r")
basis_th = read(hfile, "basis")
close(hfile)
basis_th = basis_th[:, 1:20]

basis_th_shape = shape_data(basis_th[M*N+1:end,:], M, N, true)

begin
    fig = Figure(resolution = (2500, 1000))
    for i in 1:2
        for j in 1:5
            ind = (i-1)*5+j
            ax = Axis(fig[i,j], title="Mode $(ind)")
            heatmap!(ax, basis_th_shape[:,:,ind])
        end
    end
    # save("figs/$(parent_folder)/basis_th_$(d)d.png", fig)
    display(fig)
end


####

function get_ens_vars(d, true_ens_gmt; get_means=false, k=2) # OR the means lol
    M = 192
    N = 96
    L = 1032
    hfile = h5open("data/$(parent_folder)/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
    mean_coefs = read(hfile, "mean_coefs_$(k)")
    chol_coefs = read(hfile, "chol_coefs")
    basis = read(hfile, "basis")
    close(hfile)
    ens_vars_tas = zeros(M, N, L)
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
                ens_vars_tas[:,:,(m-1)*12+n] = shape_data(data[1:M*N,:], M, N)
            end
        else
            co = get_cov(true_ens_gmt[m], chol_coefs) 
            for n in 1:12
                data = sum([co[:,:,n][i,j].*basis[:,i].*basis[:,j] for i in 1:d, j in 1:d]) 
                ens_vars_tas[:,:,(m-1)*12+n] =  shape_data(data[1:M*N,:], M, N)
            end
        end
    end
    println("done")
    return ens_vars_tas
end

param = "k"
d = 20

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

        ens_vars_tas = get_ens_vars(d, true_ens_gmt)
        ens_means_tas = get_ens_vars(d, true_ens_gmt; get_means=true)
        hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_$(scenario)_$(d)d.hdf5", "w")
        write(hfile, "ens_vars_tas_$(d)", ens_vars_tas)
        write(hfile, "ens_means_tas_$(d)", ens_means_tas)
        close(hfile)
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
            ens_means_tas = get_ens_vars(d, true_ens_gmt; get_means=true, k=k)
            write(hfile, "ens_means_tas_k$(k)", ens_means_tas)
        end
        write(hfile, "d", d)
        close(hfile)
    end
end


###



function calculate_rmse(numbers, true_variable, scenarios; rel_error=false, for_k=false)
    # println("working on $(variable) and d = $(numbers[1]) and for_k is $(for_k)")
    # flush(stdout) 
 
    variable = true_variable == "temp" ? "tas" : "pr" #pr is just the label at this point, could be whatever the second variable is

    for scenario in scenarios[2:end]
        # once we have the emulator and have saved out the emulator ensemble vars/means etc; calculate the rmse compared to various scenarios

        hfile = h5open("data/ground_truth/vars_$("huss")_$(scenario)_50ens.hdf5", "r") # true CMIP vars
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

for variable in ["temp"]
    rel_error = false

    numbers = [20]
    calculate_rmse(numbers, variable, scenarios; rel_error=rel_error, for_k=false)
    
    ks = [x for x in 1:2]
    calculate_rmse(ks, variable, scenarios; rel_error=rel_error, for_k=true)
end


## go define the rmse function then come back

numbers = [20]#[10, 100]
ks = [x for x in 1:2]

variable = "temp" #temp/pr
begin 
    fig = Figure(resolution=(1500,1000)) #
    # lims = Dict("temp" => (0.15, 0.6), "pr" => (3e-6, 9e-6))

    rel_error = false # if true, remove ylims settings
    measure = "mean"
    ax = Axis(fig[1,1:4], title="a) Average RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") for varied # of modes", xlabel="Year", ylabel="RMSE")
    # ylims!(ax, lims[variable])
    plot_rmse(ax, variable, measure, numbers; rel_error=rel_error)
    ax = Axis(fig[1,5:8], title="b) Average RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") for varied degree of fit", xlabel="Year")
    # ylims!(ax, lims[variable])
    plot_rmse(ax, variable, measure, ks; rel_error=rel_error, testing_k=true)
    

    ax = GeoAxis(fig[2,1:5], title="d) RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") on SSP119")
    hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$("ssp119").hdf5", "r") #toggle here
    begin
        data = rel_error ? read(hfile, "rmse_$(measure)s_$(variable)_100_rel") : read(hfile, "rmse_$(measure)s_$(variable == "temp" ? "tas" : variable)_20") #CHANGE BACK
        close(hfile)
        ext = (0., maximum(data))
        heatmap!(ax, lonvec2, latvec, data[:,:,1], colormap=:thermal, colorrange=ext)
        Colorbar(fig[2,6], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    end 

    measure = "std"
    ax = Axis(fig[1,9:12], title="c) Average RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") for varied # of modes", xlabel="Year")
    plot_rmse(ax, variable, measure, numbers; rel_error=rel_error)
    
    ax = GeoAxis(fig[2,7:11], title="e) RMSE of the ensemble $(measure) \n for $(variable == "temp" ? "temperature" : "specific humidity") on SSP119")
    hfile = h5open("data/$(parent_folder)/ens_vars/ens_vars_rmse_$("ssp119").hdf5", "r") #toggle here
    begin
        data = rel_error ? read(hfile, "rmse_$(measure)s_$(variable)_100_rel") : read(hfile, "rmse_$(measure)s_$(variable == "temp" ? "tas" : variable)_20") #CHANGE BACK
        close(hfile)
        ext = (0., maximum(data))
        heatmap!(ax, lonvec2, latvec, data[:,:,1], colormap=:thermal, colorrange=ext)
        Colorbar(fig[2,12], label="RMSE", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    end
    colsize!(fig.layout, 12, Relative(1/11))
    # save("figs/$(parent_folder)/rmse_joint_$(variable).png", fig)
    fig
end 

###########

hfile = h5open("data/ground_truth/vars_huss_ssp585_50ens.hdf5", "r")
true_var = read(hfile, "true_var")
true_ens_mean = read(hfile, "true_ens_mean")[:,:,:,1]
close(hfile)

avg_std = sum(sqrt.(true_var), dims=3)[:,:,1]./size(true_var)[3]
avg_mean = sum(true_ens_mean, dims=3)[:,:,1]./size(true_ens_mean)[3]

begin
    fig = Figure(resolution=(1500, 1000))
    ax = GeoAxis(fig[1,1])
    heatmap!(ax, lonvec2, latvec, avg_std, colormap=:thermal)
    ext = extrema(avg_std)
    Colorbar(fig[1,2], label="Standard deviation", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    ax = GeoAxis(fig[2,1])
    heatmap!(ax, lonvec2, latvec, avg_mean, colormap=:thermal)
    ext = extrema(avg_mean)
    Colorbar(fig[2,2], label="Mean", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    # save("figs/ground_truth/huss_std_mean.png", fig)
    display(fig)
end




hfile = h5open("data/$parent_folder/ens_vars/ens_vars_ssp585_20d.hdf5", "r")
ens_means_tas = read(hfile, "ens_means_tas_20")
ens_vars_tas = read(hfile, "ens_vars_tas_20")
close(hfile)

avg_std = sum(sqrt.(ens_vars_tas), dims=3)[:,:,1]./size(true_var)[3]
avg_mean = sum(ens_means_tas, dims=3)[:,:,1]./size(true_ens_mean)[3]

begin
    fig = Figure(resolution=(1500, 1000))
    ax = GeoAxis(fig[1,1])
    heatmap!(ax, lonvec2, latvec, avg_std, colormap=:thermal)
    ext = extrema(avg_std)
    Colorbar(fig[1,2], label="Standard deviation", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    ax = GeoAxis(fig[2,1])
    heatmap!(ax, lonvec2, latvec, avg_mean, colormap=:thermal)
    ext = extrema(avg_mean)
    Colorbar(fig[2,2], label="Mean", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    # save("figs/$parent_folder/ens_huss_std_mean.png", fig)
    display(fig)
end


#

hfile = h5open("data/$parent_folder/ens_vars/ens_vars_ssp585_100d.hdf5", "r")
ens_means_tas = read(hfile, "ens_means_pr_100")
ens_vars_tas = read(hfile, "ens_vars_pr_100")
close(hfile)

avg_std = sum(sqrt.(ens_vars_tas), dims=3)[:,:,1]./size(ens_vars_tas)[3]
avg_mean = sum(ens_means_tas, dims=3)[:,:,1]./size(ens_means_tas)[3]

begin
    fig = Figure(resolution=(1500, 1000))
    ax = GeoAxis(fig[1,1])
    heatmap!(ax, lonvec2, latvec, avg_std, colormap=:thermal)
    ext = extrema(avg_std)
    Colorbar(fig[1,2], label="Standard deviation", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    ax = GeoAxis(fig[2,1])
    heatmap!(ax, lonvec2, latvec, avg_mean, colormap=:thermal)
    ext = extrema(avg_mean)
    Colorbar(fig[2,2], label="Mean", colormap=:thermal, colorrange=ext, height = Relative(2/4))
    # save("figs/temp_huss/ens_huss_std_mean_100.png", fig)
    display(fig)
end


#####

hfile = h5open("data/$parent_folder/basis_2000d.hdf5", "r")
basis = read(hfile, "basis")
temp_factor = read(hfile, "temp_factor")
pr_factor = read(hfile, "pr_factor")
close(hfile)

true_ens_mean ./ pr_factor

hfile = h5open("data/ground_truth/vars_tas_ssp585_50ens.hdf5", "r")
true_var_tas = read(hfile, "true_var")
true_ens_mean_tas = read(hfile, "true_ens_mean")[:,:,:,1]
close(hfile)

true_ens_mean_tas/temp_factor
true_ens_mean

##
file = "/net/fs06/d3/mgeo/CMIP6/raw/$(scenario)/$("hurs")/r$(1)i1p1f1/$("hurs")_Amon_MPI-ESM1-2-LR_$(scenario)_r1i1p1f1_gn_$(1850)01-$(1869)12.nc" 
ncgen(file, "hurssample.jl")