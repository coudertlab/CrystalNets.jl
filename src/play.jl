function play()
    topo_paths = [
        "/Users/joshgoldman/Documents/Research/AI4ChemS/Tilings/CrystalNets.jl/inputs/acs_sym_7_mc_4__L_3.cif",
        "/Users/joshgoldman/Documents/Research/AI4ChemS/Tilings/CrystalNets.jl/inputs/csq_sym_5_on_13_sym_8_mc_8_L_12.cif",
    ]
    opts = Options(
        bonding=Bonding.Input,
        structure=StructureType.MOF,
        clusterings=[Clustering.SingleNodes],
        export_input=false,
        export_trimmed=false,
        export_attributions=false,
        export_clusters=false,
        export_net=false,
        export_subnets=false,
        skip_minimize=false,
        dimensions=Set(3),
        throw_error=true,
        track_mapping=Vector{Int}(),
        _pos=SVector{3,Float64}[],
    )

    determine_topology(topo_paths[2], opts)
    @show opts.track_mapping

    println("Done play")
end