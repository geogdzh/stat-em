using Dates, NCDatasets


function find_closest(vec, value)
    #finds index of closest value in vector
    geq = searchsortedfirst(vec, value)
    leq = searchsortedlast(vec, value)
    if abs(vec[leq]-value) < abs(vec[geq]-value)
        return leq
    else 
        return geq
    end
end

function convert_lon(lon)
    # in: lon on a +/- scale
    # out: lon on 360
    return lon > 0 ? lon : 360 + lon
end

struct ncData{D, NV, TV, T}
    data::D
    lonvec::NV
    latvec::TV
    timevec::T
end

function ncData(file, varname)
    ds = Dataset(file)
    lonvec = ds["lon"][:]
    latvec = ds["lat"][:]
    timevec = ds["time"][:]
    data = ds[varname][:,:,:]
    return ncData(data, lonvec, latvec, timevec)
end

function slice_map(d::ncData, lon1, lon2, lat1, lat2)
    lon1 = find_closest(d.lonvec, convert_lon(lon1))
    lon2 = find_closest(d.lonvec, convert_lon(lon2))
    lat1 = find_closest(d.latvec, lat1)
    lat2 = find_closest(d.latvec, lat2)
    if lon1 > lon2
        #we're passing over the prime meridian
        dataslice = vcat(d.data[lon1:end, lat1:lat2, :], d.data[1:lon2, lat1:lat2, :])
        newlonvec = vcat(d.lonvec[lon1:end], d.lonvec[1:lon2])
    else
        dataslice = d.data[lon1:lon2, lat1:lat2, :]
        newlonvec = d.lonvec[lon1:lon2]
    end
    newlatvec = d.latvec[lat1:lat2]
    return ncData(dataslice, newlonvec, newlatvec, d.timevec)
end

(d::ncData)(lon1, lon2, lat1, lat2) = slice_map(d, lon1, lon2, lat1, lat2)

function slice_time(d::ncData, start_year, end_year)
    #inclusive
    first = findfirst(x -> x == DateTime(start_year, 1, 1), d.timevec)
    last = findfirst(x -> x == DateTime(end_year, 12, 1), d.timevec)

    dataslice = d.data[:,:,first:last]
    newtimevec = d.timevec[first:last]
    return ncData(dataslice, d.lonvec, d.latvec, newtimevec)

end

(d::ncData)(start_year, end_year) = slice_time(d, start_year, end_year)

