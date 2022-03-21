using Documenter
using CrystalNets

# from LiveServer + Literate documentation
# using Literate
# const src = joinpath(@__DIR__, "src")
# const lit = joinpath(@__DIR__, "lit")
# for (root, _, files) ∈ walkdir(lit), file ∈ files
#     splitext(file)[2] == ".jl" || continue
#     ipath = joinpath(root, file)
#     opath = splitdir(replace(ipath, lit=>src))[1]
#     Literate.markdown(ipath, opath)
# end

DocMeta.setdocmeta!(CrystalNets, :DocTestSetup, quote
    using CrystalNets
    import CrystalNets: Options, Clustering, Bonding, StructureType
    const PeriodicGraphs = CrystalNets.PeriodicGraphs
    using .PeriodicGraphs

    CrystalNets.toggle_export(false)
    CrystalNets.toggle_warning(false)
end; recursive=true)

makedocs(
    sitename = "CrystalNets.jl",
    format = Documenter.HTML(),
    modules = [CrystalNets],
    pages = [
        "Home" => "index.md",
        "Visualization" => "man/visualization.md",
        "FAQ" => "faq.md",
        "Library" => [
            "Public"    => "lib/public.md",
            "Internals" => "lib/internals.md",
        ]
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/coudertlab/CrystalNets.jl.git"
)
