########### covariance fits
d = 100
scenario = "ssp585"

hfile = h5open("data/process/chols_$(scenario)_$(d)d.hdf5", "r")
chols = read(hfile, "chols")
chol_coefs = read(hfile, "chol_coefs")
ens_gmt = read(hfile, "ens_gmt")[:]
close(hfile)

hfile = h5open("data/process/covs_$(scenario)_$(d)d.hdf5")
covs = zeros(200, 200, 12, length(ens_gmt)) 
for i in 1:length(ens_gmt)
    covs[:,:,:,i] = read(hfile, "covs_$i")
end

covhats = zeros(200, 200,12, length(ens_gmt))
for j in 1:12
    for (i,x) in enumerate(ens_gmt)
        U = chol_coefs[j, :, :, 2].*x .+ chol_coefs[j, :, :, 1]
        covhats[:,:,j,i] = U' * U
    end
end

begin
    for shift in [0, 95]
        fig = Figure(resolution=(1000, 1000))
        for i in 1:5
            for j in 1:5
                ax = Axis(fig[i, j], xticklabelrotation=45, title=(i==1 ? "Mode $(j+shift)" : ""), xlabel = (i==5 ? "GMT" : ""), ylabel = (j==1 ? "Mode $(i+shift)" : ""), ylabelfont=:bold)
                scatter!(ax, ens_gmt, covs[i+shift, j+shift, 1, :], color=:orange, alpha=0.5)
                lines!(ax, ens_gmt, covhats[i+shift, j+shift, 1, :], color=:blue)
            end
        end
        save("figs/$parent_folder/covariance_fits_$(shift==0 ? "top" : "bottom").png", fig)
        display(fig)
    end
end


##################

# get the EOF coefficients
hfile = h5open("data/$parent_folder/training_data_ssp585_$(d)d_$(num_ens_members)ens.hdf5", "r")
ens_gmt = read(hfile, "ens_gmt")
ens_projts = read(hfile, "projts")
close(hfile)


# get means and covs and marginalize them
hfile = h5open("data/$parent_folder/gaussian_emulator_ssp585_$(d)d.hdf5", "r")
means_coefs = read(hfile, "mean_coefs_2")
close(hfile)

muhats = zeros(100, 12, length(ens_gmt))
for j in 1:12
    for (i,x) in enumerate(ens_gmt)
        muhats[:,j,i] = means_coefs[j, :, 3].*x^2 + means_coefs[j, :, 2].*x .+ means_coefs[j, :, 1]
    end
end

begin   
    sample_length = 10
    month = 1
    fig = Figure(resolution=(1000, 1000))

    for (a, i) in enumerate(1:5)
        for (b, j) in enumerate(6:10)
            ax = Axis(fig[a,b], xlabel="Mode $i", ylabel="Mode $j") #title=(a==1 ? "Mode $j" : "")

            samplesA = reshape(ens_projts[i, month:12:12*sample_length, :], (sample_length*num_ens_members))
            samplesB = reshape(ens_projts[j, month:12:12*sample_length, :], (sample_length*num_ens_members))
            scatter!(ax, samplesA, samplesB, color=:red)

            params = []
            for n in 1:sample_length
                mu  = [muhats[i, month, n], muhats[j, month, n]]
                covh = [covs[i, i, month, n] covs[i, j, month, n]; covs[j, i, month, n] covs[j, j, month, n]]
                push!(params, (mu, covh))
            end
            dist = MixtureModel(MvNormal, params)
            samples = rand(dist, 100)
            scatter!(ax, samples[1,:], samples[2,:], color=:black)

            hidedecorations!(ax, label=false)
        end
    end
    save("figs/$parent_folder/eof_gaussianity.png", fig)
    display(fig)
end

