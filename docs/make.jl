using Documenter
using CrystalNets

DocMeta.setdocmeta!(CrystalNets, :DocTestSetup, quote
    using CrystalNets
    import CrystalNets: Options
end; recursive=true)

makedocs(
    sitename = "CrystalNets",
    format = Documenter.HTML(),
    modules = [CrystalNets]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
