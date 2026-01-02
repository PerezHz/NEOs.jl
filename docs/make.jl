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
    ],
    plugins = [bib]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
