# Example: plot NEOs' custom objects leveraging the recipes in
# ext/NEOsRecipesBaseExt.jl
using NEOs, Plots

# Download optical astrometry of asteroid 2023 DW
radec = fetch_radec_mpc("2023 DW")

# Orbit determination
params = NEOParameters(bwdoffset = 0.007, fwdoffset = 0.007)
sol = orbitdetermination(radec, params)

# Plot optical residuals
plot(sol.res, label = "", xlabel = "αcos(δ) [arcsec]", ylabel = "δ [arcsec]",
     title = "2023 DW optical O-C residuals")

# Plot asteroid's orbit (in barycentric cartesian coordinates)
# N: number of points between t0 and tf (100 by default)
# projection: which coordinates to plot, options are: :x, :y, :z, :xy,
#             :xz, :yz and :xyz (default)
t0, tf = sol.bwd.t0 + sol.bwd.t[end], sol.fwd.t0 + sol.fwd.t[end]
plot(sol, t0, tf, N = 10_000, projection = :xyz, label = "2023 DW",
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
A = AdmissibleRegion(sol.tracklets[1], params)
plot(A, N = 10_000, ρscale = :log, label = "", xlabel = "log10(ρ) [au]",
     ylabel = "v_ρ [au/day]", title = "2023 DW first tracklet AR",
     c = :red, linewidth = 2)

