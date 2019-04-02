struct RadarDataJPL{T}
    object::String
    utcepoch::DateTime
    delay::T
    delay_sigma::T
    delay_units::String
    doppler::T
    doppler_sigma::T
    doppler_units::String
    freq::T
    rcvr::Int
    xmit::Int
    bouncepoint::String
    #Inner constructor
    function RadarDataJPL{T}(
        object::String,
        utcepoch::DateTime,
        delay::T,
        delay_sigma::T,
        delay_units::String,
        doppler::T,
        doppler_sigma::T,
        doppler_units::String,
        freq::T,
        rcvr::Int,
        xmit::Int,
        bouncepoint::String
    ) where {T<:Number}
        return new{T}(object, utcepoch, delay, delay_sigma, delay_units,
            doppler, doppler_sigma, doppler_units, freq, rcvr, xmit, bouncepoint)
    end
end

function RadarDataJPL(object::String,
    utcepoch::DateTime,
    delay::T,
    delay_sigma::T,
    delay_units::String,
    doppler::T,
    doppler_sigma::T,
    doppler_units::String,
    freq::T,
    rcvr::Int,
    xmit::Int,
    bouncepoint::String
) where {T<:Number}
    return RadarDataJPL{T}(object, utcepoch, delay, delay_sigma, delay_units,
        doppler, doppler_sigma, doppler_units, freq, rcvr, xmit, bouncepoint)
end

function process_radar_data_jpl(radar_data_jpl_file_path::String)
    radar_data_jpl = readdlm(radar_data_jpl_file_path, '\t')
    dataformat_jpl = "y-m-d H:M:S"
    jpl_datetime = DateTime.(radar_data_jpl[:,2], dataformat_jpl)
    del_ind = findall(x->x=="us", radar_data_jpl[:,5])
    dop_ind = findall(x->x=="Hz", radar_data_jpl[:,5])
    jpl_datetime_unrep = sort(
    union(
        DateTime.(radar_data_jpl[del_ind,2], dataformat_jpl),
        DateTime.(radar_data_jpl[dop_ind,2], dataformat_jpl)
        )
    )
    radar_data_jpl_processed = Array{Any}(undef, length(jpl_datetime_unrep), 12)
    for i in eachindex(jpl_datetime_unrep)
        inds = findall(x->x==jpl_datetime_unrep[i], jpl_datetime)
        if length(inds) == 1
            if radar_data_jpl[inds[1],5] == "us"
                x = vcat(radar_data_jpl[inds[1],1], jpl_datetime_unrep[i], radar_data_jpl[inds[1],3:5], [NaN, NaN, "X"], radar_data_jpl[inds[1], 6:end])
            elseif radar_data_jpl[inds[1],5] == "Hz"
                x = vcat(radar_data_jpl[inds[1],1], jpl_datetime_unrep[i], [NaN, NaN, "X"], radar_data_jpl[inds[1],3:5], radar_data_jpl[inds[1], 6:end])
            else
                @warn "Measurement units not us nor Hz"
            end
            radar_data_jpl_processed[i,:] .= x
        elseif length(inds) == 2
            if radar_data_jpl[inds[1],5] == "us"
                y = vcat(radar_data_jpl[inds[1],1], jpl_datetime_unrep[i], radar_data_jpl[inds[1],3:5], radar_data_jpl[inds[2],3:5], radar_data_jpl[inds[1], 6:end])
            elseif radar_data_jpl[inds[1],5] == "Hz"
                y = vcat(radar_data_jpl[inds[1],1], jpl_datetime_unrep[i], radar_data_jpl[inds[2],3:5], radar_data_jpl[inds[1],3:5], radar_data_jpl[inds[1], 6:end])
            else
                @warn "Measurement units not us nor Hz"
            end
            radar_data_jpl_processed[i,:] .= y
        else
            @warn "length(inds) != 1 && length(inds) != 2"
        end
    end
    radar_data_vec = RadarDataJPL.(
        String.(radar_data_jpl_processed[:,1]),
        DateTime.(radar_data_jpl_processed[:,2]),
        Float64.(radar_data_jpl_processed[:,3]),
        Float64.(radar_data_jpl_processed[:,4]),
        String.(radar_data_jpl_processed[:,5]),
        Float64.(radar_data_jpl_processed[:,6]),
        Float64.(radar_data_jpl_processed[:,7]),
        String.(radar_data_jpl_processed[:,8]),
        Float64.(radar_data_jpl_processed[:,9]),
        Int.(radar_data_jpl_processed[:,10]),
        Int.(radar_data_jpl_processed[:,11]),
        String.(radar_data_jpl_processed[:,12])
    )
    return radar_data_vec
end