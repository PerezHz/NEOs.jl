#eastlong must be in degrees; result is in degrees
function localsiderealtime(jd::Real, east_long)
    J0 = floor(jd)-0.5 # Julian date at 0 h UT
    T0 = (J0-2.451545e6)/36525 # Julian centuries since J2000.0 = 2.451545e6
    # Greenwich sidereal time at T0 (degrees)
    θG0 = mod(100.4606184 + 36000.77004T0 + 0.000387933T0^2 - 2.583E-8T0^3, 360.0)
    UT = (jd-J0) # UT (in units of days)
    θG = θG0 + 360.98564724UT # Greenwich sidereal time at jd (degrees)
    return mod(θG + east_long, 360.0)
end

function localsiderealtime(jd::Taylor1{T}, east_long) where {T<:Real}
    J0 = Taylor1([floor(constant_term(jd))-0.5, one(eltype(jd))], jd.order) # Julian date at 0 h UT
    T0 = (J0-2.451545e6)/36525 # Julian centuries since J2000.0 = 2.451545e6
    # Greenwich sidereal time at T0 (degrees)
    θG0 = mod(100.4606184 + 36000.77004T0 + 0.000387933T0^2 - 2.583E-8T0^3, 360.0)
    UT = (jd-J0) # UT (in units of days)
    θG = θG0 + 360.98564724UT # Greenwich sidereal time at jd (degrees)
    return mod(θG + east_long, 360.0)
end

function observer_position(station_code, t_utc::T) where {T<:Real}
    # east_long: East longitude (deg)
    # r_cos_phi: distance from spin axis (km), taken from Yeomans et al. (1992)
    # r_sin_phi: height above equatorial plane (km), taken from Yeomans et al. (1992)
    # East long data is consistent with MPC database (2019-Mar-7)
    # Goldstone site numbers 13, 14 seem to be interchanged in Yeomans et al. (1992) paper
    # TODO: add more radar stations
    if station_code == 251 # Arecibo
        east_long = 293.24692 #deg
        r_cos_phi = 6056.525 #km
        r_sin_phi = 1994.665 #km
    elseif station_code == 252 # Goldstone DSS 13 (Venus site), Fort Irwin
        east_long = 243.20512 #deg
        r_cos_phi = 5215.484 #km
        r_sin_phi = 3660.957 #km
    elseif station_code == 253 # Goldstone DSS 14 (Mars site), Fort Irwin
        east_long = 243.11047 #deg
        r_cos_phi = 5203.997 #km
        r_sin_phi = 3677.052 #km
    elseif station_code == 254 # Haystack, Westford, MA
        east_long = 288.51128 #deg
        r_cos_phi = 4700.514 #km
        r_sin_phi = 4296.900 #km
    else
        @error "Unknown station."
    end
    # local sidereal time of observation
    lmst = deg2rad( localsiderealtime(t_utc, east_long) )
    # cartesian components of topocentric position of observer
    x_gc = r_cos_phi*cos(lmst)
    y_gc = r_cos_phi*sin(lmst)
    z_gc = r_sin_phi*one(lmst)

    # rotate observer topocentric position from Earth frame to inertial frame
    t_tdb = julian(TDBEpoch(UTCEpoch(t_utc, origin=:julian))).Δt
    W = inv(earth_pole_rotation(t_tdb-J2000))

    vec_pos_pod = [x_gc, y_gc, z_gc]
    vec_pos_frame = W*vec_pos_pod
    @show vec_pos_pod
    @show vec_pos_frame

    return vec_pos_frame
end

function observer_position(station_code, t_utc::Taylor1{T}) where {T<:Number}
    # east_long: East longitude (deg)
    # r_cos_phi: distance from spin axis (km), taken from Yeomans et al. (1992)
    # r_sin_phi: height above equatorial plane (km), taken from Yeomans et al. (1992)
    # East long data is consistent with MPC database (2019-Mar-7)
    # Goldstone site numbers 13, 14 seem to be interchanged in Yeomans et al. (1992) paper
    # TODO: add more radar stations
    if station_code == 251 # Arecibo
        east_long = 293.24692 #deg
        r_cos_phi = 6056.525 #km
        r_sin_phi = 1994.665 #km
    elseif station_code == 252 # Goldstone DSS 13 (Venus site), Fort Irwin
        east_long = 243.20512 #deg
        r_cos_phi = 5215.484 #km
        r_sin_phi = 3660.957 #km
    elseif station_code == 253 # Goldstone DSS 14 (Mars site), Fort Irwin
        east_long = 243.11047 #deg
        r_cos_phi = 5203.997 #km
        r_sin_phi = 3677.052 #km
    elseif station_code == 254 # Haystack, Westford, MA
        east_long = 288.51128 #deg
        r_cos_phi = 4700.514 #km
        r_sin_phi = 4296.900 #km
    else
        @error "Unknown station."
    end
    # local sidereal time of observation
    lmst = deg2rad( localsiderealtime(t_utc, east_long) )
    # cartesian components of topocentric position of observer
    x_gc = r_cos_phi*cos(lmst)
    y_gc = r_cos_phi*sin(lmst)
    z_gc = r_sin_phi*one(lmst)

    # rotate observer topocentric position from Earth frame to inertial frame
    t_tdb = julian(TDBEpoch(UTCEpoch(constant_term(t_utc), origin=:julian))).Δt
    W = inv( earth_pole_rotation(Taylor1([t_tdb, one(t_tdb)], t_utc.order)-J2000) )

    vec_pos_pod = [x_gc, y_gc, z_gc]
    vec_pos_frame = W*vec_pos_pod
    # @show vec_pos_pod
    # @show vec_pos_frame

    return vec_pos_frame
end
