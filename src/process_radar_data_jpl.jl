@doc raw"""
    RadarDataJPL{T}

Radar measurements.

# Fields

- `object::String`: object name.
- `utcepoch::DateTime`: UTC time. 
- `delay::T`: time-delay.
- `delay_sigma`: time-delay uncertainty.
- `delay_units::String`: units of time delay. 
- `doppler`: Doppler shift.
- `doppler_sigma`: Doppler shift uncertainty.
- `doppler_units`: units of Doppler shift.
- `freq`: frequency of the measurement.
- `rcvr`: ID of reciever antenna.
- `xmit`: ID of emission antenna. 
- `bouncepoint`: bounce point. 
"""
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
    # Inner constructor
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

# Outer constructors
RadarDataJPL(r::RadarDataJPL{T}) where {T<:Number} = r
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

# Functions to get specific fields ofe a RadarDataJPL object 

utcepoch(rdata::RadarDataJPL) = rdata.utcepoch
delay(rdata::RadarDataJPL) = rdata.delay
delay_sigma(rdata::RadarDataJPL) = rdata.delay_sigma
delay_units(rdata::RadarDataJPL) = rdata.delay_units
doppler(rdata::RadarDataJPL) = rdata.doppler
doppler_sigma(rdata::RadarDataJPL) = rdata.doppler_sigma
doppler_units(rdata::RadarDataJPL) = rdata.doppler_units
freq(rdata::RadarDataJPL) = rdata.freq
rcvr(rdata::RadarDataJPL) = rdata.rcvr
xmit(rdata::RadarDataJPL) = rdata.xmit
bouncepoint(rdata::RadarDataJPL) = rdata.bouncepoint

@doc raw"""
    ismonostatic(rdata::RadarDataJPL)

Checks whether the setup is monostatic, i.e., receiver and transmitter are the same. 
"""
ismonostatic(rdata::RadarDataJPL) = rdata.rcvr == rdata.xmit

@doc raw"""
    issband(f_MHz::T) where {T<:Real}

Checks whether the transmission frequency `f_MHz` (in MHz) belongs to the S band 
(IEEE nomenclature). 
"""
issband(f_MHz::T) where {T<:Real} = 2000.0 ≤ f_MHz ≤ 4000.0

@doc raw"""
    isxband(f_MHz::T) where {T<:Real}

Checks whether the transmission frequency `f_MHz` (in MHz) belongs to the X band 
(IEEE nomenclature). 
"""
isxband(f_MHz::T) where {T<:Real} = 8000.0 ≤ f_MHz ≤ 12000.0

@doc raw"""
    process_radar_data_jpl(radar_data_jpl_file_path::String)

Parses the file of radar data `radar_data_jpl` and returns a vector of 
`RadarDataJPL` objects, each one corresponding to one measurement in `radar_data_jpl`. 
"""
function process_radar_data_jpl(radar_data_jpl_file_path::String)
    # Read radar data file 
    radar_data_jpl = readdlm(radar_data_jpl_file_path, '\t')
    # Desired format of date
    dataformat_jpl = "y-m-d H:M:S"
    # Transform date (second column) to desired format 
    jpl_datetime = DateTime.(radar_data_jpl[:,2], dataformat_jpl)
    # Find time delays (rows in units of us) 
    del_ind = findall(x->x=="us", radar_data_jpl[:,5])
    # Find Doppler shifts (rows in units of Hz)
    dop_ind = findall(x->x=="Hz", radar_data_jpl[:,5])
    # Sort dates
    jpl_datetime_unrep = sort(
    union(
        DateTime.(radar_data_jpl[del_ind,2], dataformat_jpl),
        DateTime.(radar_data_jpl[dop_ind,2], dataformat_jpl)
        )
    )
    # Allocate memmory for rada data 
    radar_data_jpl_processed = Array{Any}(undef, length(jpl_datetime_unrep), 12)

    # Iterate over each row 
    for i in eachindex(jpl_datetime_unrep)
        # Searh the i-th row (in sorted)
        inds = findall(x->x==jpl_datetime_unrep[i], jpl_datetime)
        # If there is only one measurement in the corresponding date...
        if length(inds) == 1
            # If the measurement is a time delay 
            if radar_data_jpl[inds[1],5] == "us"
                # i-th row 
                x = vcat(radar_data_jpl[inds[1],1], jpl_datetime_unrep[i], radar_data_jpl[inds[1],3:5], [NaN, NaN, "X"], radar_data_jpl[inds[1], 6:end])
            # Elseif the measurement is a Doppler shift
            elseif radar_data_jpl[inds[1],5] == "Hz"
                # i-th row 
                x = vcat(radar_data_jpl[inds[1],1], jpl_datetime_unrep[i], [NaN, NaN, "X"], radar_data_jpl[inds[1],3:5], radar_data_jpl[inds[1], 6:end])
            # Else the measurement is neither a time delay or a Doppler shift 
            else
                @warn "Measurement units not us nor Hz"
            end
            # Add i-th row 
            radar_data_jpl_processed[i,:] .= x
        # Else if there are two measurements in the corresponding date...
        elseif length(inds) == 2
            # First measurement is a time delay, second is a Doppler shift
            if radar_data_jpl[inds[1],5] == "us"
                # i-th row 
                y = vcat(radar_data_jpl[inds[1],1], jpl_datetime_unrep[i], radar_data_jpl[inds[1],3:5], radar_data_jpl[inds[2],3:5], radar_data_jpl[inds[1], 6:end])
            # First measurement is a Doppler shift, second is a time delay
            elseif radar_data_jpl[inds[1],5] == "Hz"
                # i-th row 
                y = vcat(radar_data_jpl[inds[1],1], jpl_datetime_unrep[i], radar_data_jpl[inds[2],3:5], radar_data_jpl[inds[1],3:5], radar_data_jpl[inds[1], 6:end])
            # neither of the two former cases 
            else
                @warn "Measurement units not us nor Hz"
            end
            # Add i-th row 
            radar_data_jpl_processed[i,:] .= y
        # Neither one or two measurements in the corresponding date 
        else
            @warn "length(inds) != 1 && length(inds) != 2"
        end
    end
    # Return vector of RadarDataJPL objects
    radar_data_vec = RadarDataJPL.(
        String.(radar_data_jpl_processed[:,1]),     # Name of the object
        DateTime.(radar_data_jpl_processed[:,2]),   # UTC time 
        Float64.(radar_data_jpl_processed[:,3]),    # time delay 
        Float64.(radar_data_jpl_processed[:,4]),    # time delay uncertainty
        String.(radar_data_jpl_processed[:,5]),     # time delay units
        Float64.(radar_data_jpl_processed[:,6]),    # Doppler shift
        Float64.(radar_data_jpl_processed[:,7]),    # Doppler shift uncertainty
        String.(radar_data_jpl_processed[:,8]),     # Doppler shift units 
        Float64.(radar_data_jpl_processed[:,9]),    # Frequency
        Int.(radar_data_jpl_processed[:,10]),       # ID of reciever antenna
        Int.(radar_data_jpl_processed[:,11]),       # ID of transmitter antenna
        String.(radar_data_jpl_processed[:,12])     # Bounce point 
    )
    return radar_data_vec
end