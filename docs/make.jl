using Documenter
using DocumenterCitations
using NEOs

bib = CitationBibliography(joinpath(@__DIR__, "neos.bib"), style = :numeric)

makedocs(
    sitename = "NEOs.jl",
    format = Documenter.HTML(
        assets = String["assets/citations.css"],
    ),
    modules = [NEOs],
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Astrometry" => "astrometry.md",
        "Propagation" => "propagation.md",
        "Orbit determination" => "orbitdetermination.md",
        "Impact monitoring" => "impactmonitoring.md",
        "Advanced examples" => [
            "NEOCP" => "neocp.md",
        ],
        "References" => "references.md",
    ],
    plugins = [bib]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
