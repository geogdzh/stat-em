using CairoMakie, ProgressBars, HDF5, GeoMakie
include("utils/data_util.jl")
include("utils/eof_util.jl")
include("utils/emulator_util.jl") 

parent_folder = "hurs_only"

file_head = "/net/fs06/d3/mgeo/CMIP6/interim/"

file_tail = "historical/hurs/r10i1p1f1_historical_hurs.nc"
ts = ncData(file_head*file_tail, "hurs")
M, N, L1 = size(ts.data)

X = reshape_data(ts.data)

file_tail = "ssp585/hurs/r10i1p1f1_ssp585_hurs.nc"
ts2 = ncData(file_head*file_tail, "hurs")
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

###