using LinearAlgebra, Statistics

function reshape_data(data)
    M, N, L = size(data)
    return reshape(data, (M * N, L))
end

reshape_data(d::ncData) = reshape_data(d.data)

function shape_data(data, M, N, full_array=false)
    if !full_array
        return reshape(data, (M, N))
    else
        L = size(data)[2]
        return reshape(data, (M,N, L))
    end
end


# function get_eof(d::ncData, full=false)
#     X = reshape_data(d.data)
#     U, S, V = svd(X) 
#     return full ? U : U,S,V
# end


function projection(v, U)
    #for an orthogonal basis
    dim = size(U,2)
    proj = zeros(dim)
    for i in 1:dim
        proj[i] = dot(v, U[:, i])/dot(U[:,i], U[:,i])
    end
    return proj
end

function project_timeseries(data, U; reshaped=false)
    flattened = !reshaped ? reshape_data(data) : data
    timelen = size(data)[end] #!reshaped ? size(data)[3] : size(data)[2] 
    out = Array{Float64}(undef, (size(U)[2], timelen))
    for i in 1:timelen
        out[:,i] = projection(flattened[:,i], U) 
    end
    return out
end

function back_to_data(projts, basisU)
    return basisU*projts
end

function weighted_avg(data, latvec; already_weighted=false)
    #already_weighted is a flag to indicate if the data is already weighted by the area of the grid cell
    if already_weighted
        return mean(data, dims = (1,2))[:]
    end
    M, N, L = size(data)
    Δϕ = reshape(2π / M * ones(M), (M, 1, 1))
    Δθ = reshape(π / N * ones(N) .* cos.(deg2rad.(latvec)), (1, N, 1))
    metric = (Δθ .* Δϕ) / (4π)
    return sum(metric .* data, dims = (1,2))[:]
end

weighted_avg(ts::ncData) = weighted_avg(ts.data, ts.latvec; already_weighted=false)