###julia --machine-file <host-file> distributed.jl
###julia -p <number-of-processors> distributed.jl
#using Distributed
@everywhere begin
    import Pkg
    Pkg.activate("../")
    # Pkg.instantiate()
end
@everywhere begin
    using Apophis
    using Dates
    using TaylorSeries
    using SPICE: furnsh

    #script parameters (TODO: use ArgParse.jl instead)
    const objname = "Apophis"
    const maxsteps = 10000
    const nyears = 5.0
    const dense = true#false
    const apophisjlpath = dirname(pathof(Apophis))
    const radarobsfile = joinpath(apophisjlpath, "../Apophis_JPL_data_2012_2013.dat")
    const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph!
    const t0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
    const tmax = t0+365.25nyears #final time of integration
    @show t0 == 2454733.5
    @show tmax

    function earth_et(et)
        return ss16asteph( Apophis.etsecs2julian(et) )[union(3*4-2:3*4,3*(27+4)-2:3*(27+4))]
    end
    function sun_et(et)
        return ss16asteph( Apophis.etsecs2julian(et) )[union(3*1-2:3*1,3*(27+1)-2:3*(27+1))]
    end
end

# path to local Solar System ephemeris file
# ast_eph_file = joinpath(dirname(pathof(Apophis)), "../jpleph", "ss16ast343_eph_24yr_tx.jld")
ast_eph_file = joinpath(dirname(pathof(Apophis)), "../jpleph", "ss16ast343_eph_5yr_tx.jld")

ss16asteph, acc_eph, newtonianNb_Potential = Apophis.loadeph(ast_eph_file)
# Load TT-TDB from DE430 ephemeris TODO: use AstroTime.jl instead
# furnsh( joinpath(Apophis.jplephpath, "TTmTDB.de430.19feb2015.bsp") )

# aux = (1,2,3)
aux = (ss16asteph, acc_eph, newtonianNb_Potential, earth_et, sun_et)
for i in 1:nworkers()
    # @spawnat i+1 aux = (1,2,3)
    @spawnat i+1 aux = (ss16asteph, acc_eph, newtonianNb_Potential, earth_et, sun_et)
end

#warmup (compilation) short run on all processes
parallel_run(objname, dynamics, 1, t0, tmax, ast_eph_file, aux, output=true)
println("*** Finished warmup")

#Full jet transport integration until ~2038: about 8,000 steps
#parallel_run(objname, dynamics, maxsteps, t0, tmax, ast_eph_file, aux, radarobsfile=radarobsfile)
#println("*** Finished full jet transport integration")
