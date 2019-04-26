module Apophis

__precompile__(false)

export main, au, yr, observer_position, apophisdofs, sundofs, earthdofs,
    ssdofs, c_au_per_day, μ, range_ae, radvel_ae, delay_doppler,
    delay_doppler_jpleph, mas2rad, t2c_rotation_iau_00_06,
    process_radar_data_jpl, RadarDataJPL, ismonostatic,
    RNp1BP_pN_A_J234E_J2S_ng!, RNp1BP_pN_A_J234E_J2S_ng_d225!,
    RNp1BP_pN_A_J234E_J2S_ng_d225_srp!

using Reexport
@reexport using TaylorIntegration, LinearAlgebra # so that JLD may interpret previously saved Taylor1 objects saved in .jld files
@reexport using Printf
using Dates: DateTime, julian2datetime, datetime2julian
using DelimitedFiles
using Test
using JLD
using AstroTime, EarthOrientation, SOFA, SPICE
# using Statistics: mean

# integration parameters
const varorder = 10
const order = 30
const abstol = 1.0E-30

const ea = 4 #Earth's index

# vector of G*m values
const μ = [0.0002959122082855911, 4.912248045036476e-11, 7.24345233264412e-10,
8.887692445125634e-10, 1.093189450742374e-11, 9.54954869555077e-11,
2.82534584083387e-7, 8.45970607324503e-8, 1.29202482578296e-8,
1.52435734788511e-8, 1.4004765625463623e-13, 0.0]

# vector of J2*R^2 values
const Λ2 = [4.5685187392703475e-12, 0.0, 0.0, 1.9679542578489185e-12,
2.7428745500623694e-14, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

const N = length(μ)

const au = 1.495978707E8 # astronomical unit value in km
const yr = 365.25 # days in a Julian year
const daysec = 86_400 # number of seconds in a day
const c_au_per_day = daysec*(299_792.458/au) # speed of light in au per day
const c_au_per_sec = 299_792.458/au # speed of light in au per sec
const c_cm_per_sec = 100_000*299_792.458 # speed of light in cm per sec

const apophisdofs = union(34:36, 70:72)
const sundofs = union(1:3, 37:39)
const earthdofs = union(3ea-2:3ea, 3(N+ea)-2:3(N+ea))
const ssdofs = setdiff(1:72, apophisdofs)

const J2000 = 2.451545e6

# standard value of nominal mean angular velocity of Earth (rad/day), ESAA 2014 Sec 7.4.3.3 p. 296
const ω = daysec*7.292115e-5
# The relationship of the angular velocity of the earth Omega with LOD is (https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html)
# where `omega` is in units of rad/day, and `lod` is in units of milliseconds
omega(lod) = (1e-12daysec)*(72921151.467064 - 0.843994809lod)

const R_sun = 696000.0/au # Solar radius in au, value taken from DE430 docs

const A_sun = 1.06e8 # Solar corona parameter A [cm^-3] (ESAA 2014, Table 8.5 p. 329)
const a_sun = 4.89e5 # Solar corona parameter a [cm^-3] (ESAA 2014, Table 8.5 p. 329)
const b_sun = 3.91e5 # Solar corona parameter b [cm^-3] (ESAA 2014, Table 8.5 p. 329)

const S0_sun = 63.15E6 # Sun radiated power intensity at photosphere surface, Watt/meter^2
const m2_s3_to_au2_day3 = 1e-6daysec^3/au^2 # conversion factor from m^2/sec^3 to au^2/day^3

# const A_sun = 1.67e8 # Solar corona parameter A [cm^-3] (Anderson 1978)
# const a_sun = 0.750e6 # Solar corona parameter a [cm^-3] (Anderson 1978)
# const b_sun = 0.1836e6 # Solar corona parameter b [cm^-3] (Anderson 1978)

function __init__()
    @show length(methods(RNp1BP_pN_A_J234E_J2S_ng!))
    @show length(methods(TaylorIntegration.jetcoeffs!))
    @show methods(RNp1BP_pN_A_J234E_J2S_ng!)
    # load JPL ephemerides
    loadjpleph()
end

include("process_radar_data_jpl.jl")
include("jpl-de-430-431-earth-orientation-model.jl")
include("topocentric.jl")
include("asteroid_dynamical_models.jl")
include("initial_conditions.jl")
include("delay_doppler.jl")

#root-finding functions (simplified versions of range_ae, radvel_ae)
g(t,x,dx) = (x[3N-2]-x[3ea-2])^2+(x[3N-1]-x[3ea-1])^2+(x[3N]-x[3ea])^2
g2(t,x,dx) = (x[3N-2]-x[3ea-2])*(x[6N-2]-x[3(N+ea)-2])+(x[3N-1]-x[3ea-1])*(x[6N-1]-x[3(N+ea)-1])+(x[3N]-x[3ea])*(x[6N]-x[3(N+ea)])

# instantaneous geocentric range, radial velocity
range_ae(x) = sqrt( (x[3N-2]-x[3ea-2])^2+(x[3N-1]-x[3ea-1])^2+(x[3N]-x[3ea])^2 )
radvel_ae(x) = ( (x[3N-2]-x[3ea-2])*(x[6N-2]-x[3(N+ea)-2])+(x[3N-1]-x[3ea-1])*(x[6N-1]-x[3(N+ea)-1])+(x[3N]-x[3ea])*(x[6N]-x[3(N+ea)]) )/range_ae(x)

function main(objname::String, datafile::String, dynamics::Function, maxsteps::Int,
    newtoniter::Int, t0::T, tspan::T; output::Bool=true, radarobs::Bool=true,
    jt::Bool=true) where {T<:Real}

    # get initial conditions
    q0 = initialcond(length(μ))

    if jt
        #construct jet transport initial condition as Vector{Taylor1{Float64}} from `q0`
        q0T1 = Taylor1.(q0,varorder)
        q0T1[1:end-1] = Taylor1.(q0[1:end-1],varorder)
        q0T1[end] = Taylor1([q0[end],1e-14],varorder) #note the 1e-14!!!
        __q0 = q0T1
    else
        __q0 = q0
    end

    @show tmax = t0+tspan*yr #final time of integration

    # do integration
    if radarobs
        # read object radar astrometry from JPL date
        radar_data_jpl = process_radar_data_jpl(datafile)
        #construct vector of observation times (UTC) > t0
        tv_jpl_utc = UTCEpoch.([x.utcepoch for x in radar_data_jpl])
        # convert to TDB
        tv_jpl_tdb = TDBEpoch.(tv_jpl_utc)
        # date/time to Julian date
        tv_jpl_tdb_julian =  map(x->x.Δt, julian.(tv_jpl_tdb))
        # construct time range variable with t0 and observation times > t0, removing repeated values
        tv = union(t0, tv_jpl_tdb_julian[tv_jpl_tdb_julian .> t0])
        @show all(diff(tv) .> 0)
        @time sol_objs = taylorinteg(dynamics, g2, __q0, tv, order, abstol; maxsteps=maxsteps, newtoniter=newtoniter);
        tup_names = (:xv1, :tvS1, :xvS1, :gvS1)
        sol = NamedTuple{tup_names}(sol_objs)
    else
        @time sol_objs = taylorinteg(dynamics, g2, __q0, t0, tmax, order, abstol; maxsteps=maxsteps, newtoniter=newtoniter);
        tup_names = (:tv1, :xv1, :tvS1, :xvS1, :gvS1)
        sol = NamedTuple{tup_names}(sol_objs)
    end

    #write solution to .jld files
    if output
        filename = string(objname, "_jt.jld")
        println("Saving solution to file: $filename")
        #first, deal with `tv_jpl_integ`
        jldopen(filename, "w") do file
            if radarobs
                println("Saving variable: tv1")
                write(file, "tv1", tv)
            end
            #loop over variables
            for ind in eachindex(sol)
                varname = string(ind)
                println("Saving variable: ", varname)
                write(file, varname, sol[ind])
            end
        end
        #check that tv_jpl_integ was recovered succesfully
        println("Checking that all variables were saved correctly...")
        if radarobs
            recovered_tv_jpl_integ = load(filename, "tv1")
            @show recovered_tv_jpl_integ == tv
        end
        #loop over rest of variables
        for ind in eachindex(sol)
            varname = string(ind)
            #read varname from files and assign recovered variable to recovered_sol_i
            recovered_sol_i = load(filename, varname)
            #check that varname was recovered succesfully
            @show recovered_sol_i == sol[ind]
        end
        println("Saved solution")
    end
    return nothing
end

end