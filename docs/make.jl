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

deploydocs(
    repo = "github.com/PerezHz/NEOs.jl.git",
    push_preview= true,
)
