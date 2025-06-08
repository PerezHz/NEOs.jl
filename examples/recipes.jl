# Example: plot NEOs' custom objects leveraging the recipes in
# ext/NEOsRecipesBaseExt.jl
using NEOs, Plots

# Download optical astrometry of asteroid 2016 TU93
radec = fetch_radec_mpc("2016 TU93")

# Orbit determination
params = Parameters(
     coeffstol = Inf, bwdoffset = 0.042, fwdoffset = 0.042,
     gaussorder = 2, safegauss = length(od.tracklets) == 3 ? true : false,
     tsaorder = 2, adamiter = 500, adamQtol = 1e-5, jtlsorder = 2,
     jtlsmask = false, jtlsiter = 20, lsiter = 10, significance = 0.99,
     outrej = true, χ2_rec = sqrt(9.21), χ2_rej = sqrt(10),
     fudge = 100.0, max_per = 34.0,
)
od = ODProblem(newtonian!, radec)
orbit = initialorbitdetermination(od, params)

# Plot optical residuals
plot(orbit.res, label = "", xlabel = "αcos(δ) [arcsec]", ylabel = "δ [arcsec]",
     title = "2016 TU93 optical O-C residuals")

# Plot asteroid's orbit (in barycentric cartesian coordinates)
# N: number of points between t0 and tf (100 by default)
# projection: which coordinates to plot, options are: :x, :y, :z, :xy,
#             :xz, :yz and :xyz (default)
t0, tf = dtutc2days.(minmaxdates(orbit))
plot(orbit, t0, tf, N = 10_000, projection = :xyz, label = "2016 TU93",
     xlabel = "x [au]", ylabel = "y [au]", zlabel = "z [au]", color = :blue)

# We can also visualize the orbits of the Sun and the Earth
scatter!(params.eph_su, t0, tf, N = 10_000, projection = :xyz,
         label = "Sun", color = :yellow)
plot!(params.eph_ea, t0, tf, N = 10_000, projection = :xyz,
      label = "Earth", color = :green)

# Plot the boundary of an admissible region
# N: number of points along the boundary (100 by default)
# ρscale: horizontal axis scale, either :linear (default) or :log
using NEOs: AdmissibleRegion
A = AdmissibleRegion(orbit.tracklets[1], params)
plot(A, N = 10_000, ρscale = :log, label = "", xlabel = "log10(ρ) [au]",
     ylabel = "v_ρ [au/day]", title = "2016 TU93 first tracklet AR",
     c = :red, linewidth = 2)