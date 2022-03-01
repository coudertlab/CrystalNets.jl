using Documenter, Literate
using CrystalNets

# from LiveServer documentation
const src = joinpath(@__DIR__, "src")
const lit = joinpath(@__DIR__, "lit")
for (root, _, files) ∈ walkdir(lit), file ∈ files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>src))[1]
    Literate.markdown(ipath, opath)
end

DocMeta.setdocmeta!(CrystalNets, :DocTestSetup, quote
    using CrystalNets
    using PeriodicGraphs
    import CrystalNets: Options
end; recursive=true)

makedocs(
    sitename = "CrystalNets",
    format = Documenter.HTML(),
    modules = [CrystalNets],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Visualization"      => "man/visualization.md",
            "Database studies"   => "man/database.md",
            "Clustering options" => "man/clustering.md",
            "Custom archive"     => "man/archive.md",
        ],
        "Common issues" => "faq.md",
        "Library" => [
            "Public"    => "lib/public.md",
            "Internals" => "lib/internals.md",
        ]
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
