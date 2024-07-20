using LinearAlgebra, HDF5
include("../utils/data_util.jl")
include("../utils/eof_util.jl")

#### get basis - based ONLY on one ens memeber (should this be substantiated?)
using_two = true 
second_var = "huss" # "pr" or "huss"
non_dim = false  
use_metrics = false
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

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"

#first ensemble member of historical run
file_tail = "historical/tas/r10i1p1f1_historical_tas.nc"
ts = ncData(file_head*file_tail, "tas")
M, N, L1 = size(ts.data)

if use_metrics
    latvec = ts.latvec
    Δϕ = reshape(2π / M * ones(M), (M, 1, 1))
    Δθ = reshape(π / N * ones(N) .* cos.(deg2rad.(latvec)), (1, N, 1))
    metric = (Δθ .* Δϕ) / (4π)
end

X = reshape_data(use_metrics ? sqrt.(metric) .* ts.data : ts.data)

#first ens member of the model run
file_tail = "ssp585/tas/r10i1p1f1_ssp585_tas.nc"
ts2 = ncData(file_head*file_tail, "tas")
X2 = reshape_data(ts2.data)
M, N, L2 = size(use_metrics ? sqrt.(metric) .* ts2.data : ts2.data)

fullX = hcat(X, X2)

if using_two
    phfile = file_head*"historical/$(second_var)/r10i1p1f1_historical_$(second_var).nc"
    pr = ncData(phfile, second_var)
    Xp = reshape_data(use_metrics ? sqrt.(metric) .* pr.data : pr.data) # use log of precip to normalize and convert to mm/day
    if second_var == "pr"
        Xp = log.(Xp .* 86400)
    end

    phfile2 = file_head*"ssp585/$(second_var)/r10i1p1f1_ssp585_$(second_var).nc"
    pr2 = ncData(phfile2, second_var)
    Xp2 = reshape_data(use_metrics ? sqrt.(metric) .* pr2.data : pr2.data) # same as above
    if second_var == "pr"
        Xp2 = log.(Xp2 .* 86400)
    end

    fullXp = hcat(Xp, Xp2)
end

if non_dim
    temp_factor = maximum(X)
    X = X ./ temp_factor
    X2 = X2 ./ temp_factor
    if using_two
        pr_factor = maximum(Xp)
        Xp = Xp ./ pr_factor
        Xp2 = Xp2 ./ pr_factor
    end
end


full = using_two ? vcat(fullX, fullXp) : fullX
U, S, V = svd(full)
d = 2000
basis = U[:,1:d]

hfile = h5open("data/$(parent_folder)/basis_$(d)d.hdf5", "w") #but to use a smaller basis, can just take fewer modes
write(hfile, "basis", basis)
if non_dim
    write(hfile, "temp_factor", temp_factor)
    write(hfile, "pr_factor", pr_factor)
end
if use_metrics
    write(hfile, "metric", metric)
end
close(hfile)