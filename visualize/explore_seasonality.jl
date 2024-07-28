d = 100

hfile = using_two ? h5open("data/$(parent_folder)/temp_precip_basis_2000d.hdf5", "r") : h5open("data/$(parent_folder)/temp_basis_2000d.hdf5", "r") #this basis is calculated from just one ens member
basis = read(hfile, "basis")
if non_dim ##none of this is actually used here
    temp_factor = read(hfile, "temp_factor")
    pr_factor = read(hfile, "pr_factor")
end
if use_metrics
    metric = read(hfile, "metric")
end
close(hfile)
basis = basis[:, 1:d]

hfile = using_two ? h5open("data/$(parent_folder)/gaussian_emulator_withpr_ssp585_$(d)d.hdf5", "r") : h5open("data/$(parent_folder)/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
mean_coefs = read(hfile, "mean_coefs_2") #lets default to quadratic
# mean_coefs = read(hfile, "mean_coefs")
chol_coefs = read(hfile, "chol_coefs")
basis = read(hfile, "basis")
if use_metrics
    metric = read(hfile, "metric") 
end
close(hfile)

hfile = h5open("data/ssp585_gmts_50ens.hdf5")
ens_gmt = read(hfile, "ens_gmt")
close(hfile)
ens_gmt = mean(ens_gmt, dims=1)[:]


tempmeans = zeros(12, length(ens_gmt))
prmeans = zeros(12, length(ens_gmt))
tempstds = zeros(12, length(ens_gmt))
prstds = zeros(12, length(ens_gmt))
for (i, t) in enumerate(ens_gmt)
    println("working on year $(i)")
    flush(stdout)
    co = get_cov(t, chol_coefs)
    for n in 1:12
        data = back_to_data(mean_coefs[n,:,3].*t^2 .+ mean_coefs[n,:,2].*t .+ mean_coefs[n, :, 1], basis)
        tempdata = shape_data(data[1:M*N,:], M, N)
        prdata = shape_data(data[M*N+1:end,:], M, N)
        tempmeans[n, i] = tempdata[18, 83]
        prmeans[n, i] = prdata[18, 83]
        vardata = sum([co[:,:,n][i,j].*basis[:,i].*basis[:,j] for i in 1:d, j in 1:d])
        tempvardata = shape_data(vardata[1:M*N], M, N)
        prvardata = shape_data(vardata[M*N+1:end], M, N)
        tempstds[n, i] = tempvardata[18, 83]
        prstds[n, i] = prvardata[18, 83]
    end
end

hfile = h5open("data/$(parent_folder)/anchorage_$(d)d.hdf5", "cw")
write(hfile, "tempmeans", tempmeans)
write(hfile, "prmeans", prmeans)
write(hfile, "tempvars", tempstds)
write(hfile, "prvars", prstds)
close(hfile)

hfile = h5open("data/$(parent_folder)/anchorage_$(d)d.hdf5", "r")
tempmeans = read(hfile, "tempmeans")
prmeans = read(hfile, "prmeans")
tempstds = sqrt.(read(hfile, "tempvars"))
prstds = sqrt.(read(hfile, "prvars"))
close(hfile)

tempdiffs = zeros(length(ens_gmt), 3)
prdiffs = zeros(length(ens_gmt), 3)
for i in 1:length(ens_gmt)
    max, maxi = findmax(tempmeans[:,i])
    min, mini = findmin(tempmeans[:,i])
    tempdiffs[i, 2] = max - min
    tempdiffs[i, 1] = (max + tempstds[maxi, i]) - (min - tempstds[mini, i])
    tempdiffs[i, 3] = (max - tempstds[maxi, i]) - (min + tempstds[mini, i])
#
    max, maxi = findmax(prmeans[:,i])
    min, mini = findmin(prmeans[:,i])
    prdiffs[i, 2] = max - min
    prdiffs[i, 1] = (max + prstds[maxi, i]) - (min - prstds[mini, i])
    prdiffs[i, 3] = (max - prstds[maxi, i]) - (min + prstds[mini, i])
end

####plot
begin
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1:3, 1], xlabel = "Month", ylabel = "Temp", title = "Monthly Temp in Fairbanks, AK")
    lines!(ax, 1:12, tempmeans[:,1] .+ tempstds[:,1], color = :blue, linewidth = 1, alpha=0.5)
    lines!(ax, 1:12, tempmeans[:,1], color = :blue, linewidth = 2, label = "2014")
    lines!(ax, 1:12, tempmeans[:,1] .- tempstds[:,1], color = :blue, linewidth = 1, alpha=0.5)
    lines!(ax, 1:12, tempmeans[:,end] .+ tempstds[:,end], color = :red, linewidth = 1, alpha=0.5)
    lines!(ax, 1:12, tempmeans[:, end], color = :red, linewidth = 2, label = "2100")
    lines!(ax, 1:12, tempmeans[:,end] .- tempstds[:,end], color = :red, linewidth = 1, alpha=0.5)


    axislegend(ax, position = :lt)
    ax2 = Axis(fig[4,1], xlabel = "Year", ylabel = "Max Difference")
    lines!(ax2, time_future, tempdiffs[:,1], color = :black, linewidth = 1, linestyle=:dash)
    lines!(ax2, time_future, tempdiffs[:,2], color = :black, linewidth = 2)
    lines!(ax2, time_future, tempdiffs[:,3], color = :black, linewidth = 1, linestyle=:dash)
    # save("figs/fairbanks_temp.png", fig)
    fig
end

begin
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1:3, 1], xlabel = "Month", ylabel = "Precip", title = "Monthly Precip in Fairbanks, AK")
    lines!(ax, 1:12, prmeans[:,1] .+ prstds[:,1], color = :blue, linewidth = 1, alpha=0.5)
    lines!(ax, 1:12, prmeans[:,1], color = :blue, linewidth = 2, label = "2014")
    lines!(ax, 1:12, prmeans[:,1] .- prstds[:,1], color = :blue, linewidth = 1, alpha=0.5)
    lines!(ax, 1:12, prmeans[:,end] .+ prstds[:,end], color = :red, linewidth = 1, alpha=0.5)
    lines!(ax, 1:12, prmeans[:, end], color = :red, linewidth = 2, label = "2100")
    lines!(ax, 1:12, prmeans[:,end] .- prstds[:,end], color = :red, linewidth = 1, alpha=0.5)


    axislegend(ax, position = :lt)
    ax2 = Axis(fig[4,1], xlabel = "Year", ylabel = "Max Difference")
    lines!(ax2, time_future, prdiffs[:,1], color = :black, linewidth = 1, linestyle=:dash)
    lines!(ax2, time_future, prdiffs[:,2], color = :black, linewidth = 2)
    lines!(ax2, time_future, prdiffs[:,3], color = :black, linewidth = 1, linestyle=:dash)
    save("figs/fairbanks_precip.png", fig)
    fig
end
