using Statistics, Distributions, LinearAlgebra, Polynomials

### data functions
function get_gmt_list(ts::ncData; use_metrics=false)
    # applies weighting and returns monthly average over entire geographic field
    hist_mean_temp = use_metrics ? mean(ts.data, dims=(1,2))[:] : weighted_avg(ts) #this doesn't make sense!!
    hist_year_temp =  [mean(hist_mean_temp[i:min(i+12-1, end)]) for i in 1:12:length(hist_mean_temp)]
    return hist_year_temp
end

function month_to_year_avg(list)
    return [mean(list[i:min(i+12-1, end)]) for i in 1:12:length(list)]
end


### training functions
function get_mean_coefs(ens_projts, ens_gmt; degree=1)
    # takes in the ensemble of projts and the ensemble MEAN gmt
    d, ts, num_ens_members = size(ens_projts)
    mean_coefs = zeros((12, d, degree+1))
    A1 = fill(1., length(ens_gmt)*num_ens_members)
    A2 = repeat(ens_gmt', num_ens_members)
    A = hcat(A1, A2)
    for i in 1:12
        y = ens_projts[:,i:12:end,1] 
        for j in 2:num_ens_members
            y = hcat(y, ens_projts[:,i:12:end,j])
        end
        for j in 1:d   
            if degree == 1
                b = y[j,:]
                mean_coefs[i, j, :] = A \ b
            else #fix this later to not use Polynomials
                mean_coefs[i, j, :] =  Polynomials.fit(repeat(ens_gmt', num_ens_members)[:], y[j,:][:], degree)[:]
            end
        end
    end
    return mean_coefs
end

function gmt_cov(ens_projts, ens_gmt, offloading_tag; corrs=false, parent_folder=nothing)
    D = size(ens_projts)[1]*2
    num_years = length(ens_gmt)
    if isnothing(parent_folder)
        wfile = corrs == true ? h5open("data/process/corrs_$(offloading_tag).hdf5", "w") : h5open("data/process/covs_$(offloading_tag).hdf5", "w")
    else
        wfile = corrs == true ? h5open("data/process/corrs_$(offloading_tag).hdf5", "w") : h5open("data/process/$(parent_folder)/covs_$(offloading_tag).hdf5", "w")
    end
    write(wfile, "D", D)
    write(wfile, "num_years", num_years)
    for j in 1:num_years
        covs = zeros((D, D, 12))
        for i in 1:12
            proj1 = ens_projts[:,i:12:end,:] #this is all the jans, eg
            proj2 = i==12 ? ens_projts[:,1:12:end,:] : ens_projts[:,i+1:12:end,:] #all the febs
            cat = vcat(proj1, proj2)
            
            covs[:,:,i] = corrs==true ? cor(cat[:,j,:]; dims=2) : cov(cat[:,j,:]; dims=2)
        end    
        corrs == true ? write(wfile, "corrs_$j", covs) : write(wfile, "covs_$j", covs)
    end
    close(wfile)
end

# the ORIGINAL METHOD that does both -- ans is the one currently being used
function get_chol_coefs(ens_gmt, offloading_tag; return_chols=false, offload=false, parent_folder=nothing)
    hfile = isnothing(parent_folder) ? h5open("data/process/covs_$(offloading_tag).hdf5", "r") : h5open("data/process/$(parent_folder)/covs_$(offloading_tag).hdf5", "r")
    D = read(hfile, "D")
    num_years = read(hfile, "num_years")
    # NEED TO ADD IMPLEMENTATION FOR OFFLOAD: in which case call get_chols
    if !offload 
        chols = zeros((D, D, 12, num_years)) #raising errors for array size
        for i in 1:12
            for j in 1:num_years
                covs = read(hfile, "covs_$j") #this gets confusing, need to check that it's for the corresponding emulator version
                sc = covs[:,:,i]
                ll, vv = eigen(sc)
                sc = sc + sqrt(eps(ll[end])) .* I
                chols[:,:,i,j] = cholesky(sc).U
            end
        end
        chol_coefs = zeros((12, D, D, 2))
        for i in 1:12
            for j in 1:D 
                for k in 1:D
                    y = [chols[j,k,i,n] for n in 1:num_years]
                    A1 = fill(1., num_years)
                    A2 = ens_gmt'
                    A = hcat(A1, A2)
                    chol_coefs[i, j, k, :] = A \ y
                end
            end
        end
        if return_chols
            return chol_coefs, chols
        end
        return chol_coefs
    else
        get_chols(offloading_tag; parent_folder=parent_folder) #this gets just the chols and offloads them
        chol_coefs = fit_chol_coefs(ens_gmt, offloading_tag; parent_folder=parent_folder)
        return chol_coefs
    end
end

function get_chols(offloading_tag; parent_folder=nothing)
    hfile = isnothing(parent_folder) ? h5open("data/process/covs_$(offloading_tag).hdf5", "r") : h5open("data/process/$(parent_folder)/covs_$(offloading_tag).hdf5", "r")
    D = read(hfile, "D")
    num_years = read(hfile, "num_years")
    wfile = isnothing(parent_folder) ? h5open("data/process/chols_$(offloading_tag).hdf5", "w") : h5open("data/process/$(parent_folder)/chols_$(offloading_tag).hdf5", "w")
    for j in 1:num_years
        println("working on cholesky decomposition for year $j")
        flush(stdout)
        chols = zeros((D, D, 12))
        for i in 1:12 
            covs = read(hfile, "covs_$j")
            sc = covs[:,:,i]
            ll, vv = eigen(sc)
            sc = sc + sqrt(eps(ll[end])) .* I
            chols[:,:,i] = cholesky(sc).U
        end
        write(wfile, "chols_$j", chols)
    end
    # println("done with cholesky decompositions")
    write(wfile, "num_years", num_years)
    write(wfile, "D", D)
    close(wfile)
    close(hfile)
end

# function fit_chol_coefs(ens_gmt, offloading_tag)
#     hfile = h5open("data/process/chols_$(offloading_tag).hdf5", "r")
#     num_years = read(hfile, "num_years")
#     D = read(hfile, "D")
#     chol_coefs = zeros((12, D, D, 2)) 
#     for i in 1:12
#         println("working on chol_coefs for month $i")
#         flush(stdout)
#         for j in 1:D #(is already being done separately for each cell of matrix, only integrating over the years)
#             for k in 1:D
#                 y = [read(hfile, "chols_$n")[j,k,i] for n in 1:num_years] #maybe this isn't the most efficient?? 
#                 A1 = fill(1., num_years)
#                 A2 = ens_gmt'
#                 A = hcat(A1, A2)
#                 chol_coefs[i, j, k, :] = A \ y
#             end
#         end
#     end
#     println("done with chol_coefs")
#     close(hfile)
#     return chol_coefs
# end

# function fit_chol_coefs(ens_gmt, offloading_tag) # NEED TO MERGE THIS WITH ABOVE
#     hfile = h5open("data/process/chols_$(offloading_tag).hdf5", "r")
#     num_years = read(hfile, "num_years")
#     D = read(hfile, "D")
#     chol_coefs = zeros((12, D, D, 2)) 

#     A1 = fill(1., num_years)
#     A2 = ens_gmt'
#     A = hcat(A1, A2)
#     for i in 1:12
#         println("working on chol_coefs for month $i")
#         flush(stdout)
#         for j in 1:D #(is already being done separately for each cell of matrix, only integrating over the years)
#             for k in 1:D #ok well for one they're all upper/lower triangular so i don't need to run this for half the entries. it's gonna be zero
#                 println("we're at $j, $k")
#                 flush(stdout)
#                 y = zeros(num_years)
#                 for n in 1:num_years
#                     y[n] = read(hfile, "chols_$n")[j,k,i]
#                 end
#                 # y = [read(hfile, "chols_$n")[j,k,i] for n in 1:num_years] #I think this must be the issue
#                 chol_coefs[i, j, k, :] = A \ y
#             end
#         end
#     end
#     println("done with chol_coefs")
#     close(hfile)
#     return chol_coefs
# end



### testing functions

function get_cov(gmt, chol_coefs) 
    D = size(chol_coefs)[2]
    covs = zeros((D, D, 12))
    for i in 1:12
        L = chol_coefs[i, :, :, 2] .* gmt .+ chol_coefs[i, :, :, 1] #so technically this is U
        covs[:,:,i] = L'*L #but this is consistent with that 
    end
    return covs
end

function get_means(gmt, mean_coefs) 
    out = zeros((size(mean_coefs)[2],12))
    k = size(mean_coefs)[3] -1
    for i in 1:12
        if k == 1
            out[:,i] = mean_coefs[i,:,2].*gmt .+ mean_coefs[i, :, 1]
        elseif k == 2
            out[:,i] = mean_coefs[i,:,3].*gmt.^2 .+ mean_coefs[i,:,2].*gmt .+ mean_coefs[i, :, 1]
        end
    end
    return out
end

function emulate_no_cov(gmt_list, mean_coefs, chol_coefs)
    num_years = length(gmt_list)
    # gmt assumed constant for each year
    trajectory = zeros(size(mean_coefs)[2], 12*num_years) 
    for yr in 1:num_years
        gmt = gmt_list[yr]
        covs = get_cov(gmt, chol_coefs) #alternative to this (using cholesky) would be to have const corr
        for i in 1:12
            μ = get_means(gmt, mean_coefs)[:,i]
            Σ = Symmetric(covs[:,:,i][1:d,1:d]) #symmetric might not be necessary
            dist = MvNormal(μ, Σ)
            trajectory[:,(yr-1)*12 + i] = rand(dist)
        end
    end
    return trajectory
end


function emulate(gmt_list, mean_coefs, chol_coefs; no_cov=false)
    num_years = length(gmt_list)
    # gmt assumed constant for each year
    trajectory = zeros(size(mean_coefs)[2], 12*num_years) 

    gmt = gmt_list[1]
    covs = get_cov(gmt, chol_coefs)
    means = get_means(gmt, mean_coefs) #list of twelve
    Σ = covs[1:d,1:d,12]
    dec = rand(MvNormal(means[:,12], Σ))
    trajectory[:,1] = emulate_step(dec, 12, covs, means; new_means=means, no_cov=no_cov)

    for year in 1:num_years #year is an index
        for i in 1:11 #here i is the prev_month
            trajectory[:, (i+1)+(year-1)*12] = emulate_step(trajectory[:, i+(year-1)*12], i, covs, means)
        end
        if year != num_years #if that wasn't the last year
            #now set up for the next year
            gmt = gmt_list[year+1]
            new_means = get_means(gmt, mean_coefs)
            trajectory[:, (12+1)+(year-1)*12] = emulate_step(trajectory[:, 12+(year-1)*12], 12, covs, means; new_means=new_means) #do the december transition
            covs = get_cov(gmt, chol_coefs)
            means = new_means
        end
    end
    return trajectory
end


function emulate_step(prev_val, prev_month, covs, means; new_means=nothing, no_cov=false)
    # prev_val - value of preceding month
    # prev_month - integer of preceding month
    # covs - matrix of covariances for the GMT of the preceding month; size = (2d, 2d, 12) 
    # means - vector of means for the GMT of the preceding month
    # new_means - vector of means for the GMT of the new month, if different
    μ_1 = means[:,prev_month]
    μ_2 = isnothing(new_means) ? means[:,prev_month+1] : new_means[:,1] #hardcoded for switch to happen in january!
    Σ = covs[:,:,prev_month] 
    Σ += I .* sqrt(eps(maximum(Σ))) #+ sqrt(eps(maximum(covs[:,:,prev_month]))) .* I 
    d = Int(size(covs)[1]/2)
    Σ_11, Σ_12, Σ_21, Σ_22 = Σ[1:d,1:d], Σ[1:d,d+1:end], Σ[d+1:end,1:d], Σ[d+1:end,d+1:end] 
    # Σ_22 = prev_month == 12 ? covs[:,:,1][1:d,1:d] : covs[:,:,prev_month+1][1:d,1:d]
    μ_out = μ_2 .+ Σ_21 * inv(Σ_11) * (prev_val .- μ_1)
    if no_cov
        Σ_out = Σ_22
    else
        # Σ_out = Σ_22 .- Σ_21 * inv(Σ_11) * Σ_12 #this is the mathematical model
        Σ_out = Symmetric(Σ_22) - Symmetric(Σ_21 * (Σ_11 \ Σ_12))
        Σ_out = Σ_out + sqrt(eps(maximum(Σ_out))) .* I
    end
    dist = MvNormal(μ_out, Σ_out)
    return rand(dist)
end



###### for exploration/visualization purposes

function get_var_coefs(ens_projts, ens_gmt, mean_coefs; return_vars=false)
    d, ts, num_ens_members = size(ens_projts)
    t = Int(ts/12)
    var_coefs = zeros((12, d, 2))
    vars_array = zeros((12, d, num_ens_members, t))
    for i in 1:12
        y = ens_projts[:,i:12:end,1] 
        for j in 2:num_ens_members
            y = hcat(y, ens_projts[:,i:12:end,j])
        end
        for j in 1:d
            b, m = mean_coefs[i, j, :]
            fits = repeat([m*x+b for x in ens_gmt]',num_ens_members)
            vars = (y[j,:].-fits).^2                                #vars calculated as squared difference from the line of best fit
            A1 = fill(1., length(ens_gmt)*num_ens_members)
            A2 = repeat(ens_gmt', num_ens_members)
            A = hcat(A1, A2)
            var_coefs[i, j, :] = A \ vars
            vars_array[i,j,:,:] = reshape(vars, (t, num_ens_members))'
        end
    end
    if return_vars
        return var_coefs, vars_array
    end
    return var_coefs
end