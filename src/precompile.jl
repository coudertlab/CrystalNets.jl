using CrystalNets, PeriodicGraphs, PeriodicGraphEmbeddings, StaticArrays, Graphs
using SnoopPrecompile

@precompile_setup begin
    tmpdir = tempname()
    mkdir(tmpdir)
    cifs = joinpath(dirname(dirname(pathof(CrystalNets))), "test", "cif")
    path_to_im19 = joinpath(cifs, "IM-19.cif");
    path_to_rro = joinpath(cifs, "RRO.cif");
    path_to_mil100 = joinpath(cifs, "MOFs", "MIL-100.cif");
    path_to_calfig = joinpath(cifs, "CALFIG.cif");
    path_to_abw = joinpath(cifs, "ABW.cif")

    old_stdout = stdout
    old_stderr = stderr
    redirect_stdout(devnull)
    redirect_stderr(devnull)

    @precompile_all_calls begin
        # Rational{Int32}
        im19 = determine_topology(path_to_im19; bonding=Bonding.Guess, structure=StructureType.MOF,
                    clusterings=[Clustering.SingleNodes,Clustering.AllNodes,Clustering.PEM,Clustering.PE,Clustering.Standard,Clustering.Auto,Clustering.EachVertex],
                    export_input=tmpdir, export_trimmed=tmpdir, export_subnets=tmpdir, export_attributions=tmpdir, export_clusters=tmpdir)
        # Rational{Int64}
        rro = determine_topology(path_to_rro; structure=StructureType.Zeolite,
                    clusterings=[Clustering.EachVertex,Clustering.PEM,Clustering.PE,Clustering.Standard,Clustering.Auto],
                    export_input=tmpdir, export_trimmed=tmpdir, export_subnets=tmpdir, export_attributions=tmpdir, export_clusters=tmpdir)
        # Rational{Int128}
        mil100 = determine_topology(path_to_mil100; structure=StructureType.Guess,
                    clusterings=[Clustering.EachVertex,Clustering.PEM,Clustering.PE,Clustering.Standard,Clustering.Auto],
                    export_input=tmpdir, export_trimmed=tmpdir, export_subnets=tmpdir, export_attributions=tmpdir, export_clusters=tmpdir)
        # Rational{BigInt}
        sfv = topological_genome(CrystalNet(PeriodicGraph(REVERSE_CRYSTALNETS_ARCHIVE["*SFV"])))
        # 2D
        hcb = topological_genome(CrystalNet(PeriodicGraph(REVERSE_CRYSTALNETS_ARCHIVE["hcb"])))
        # non-periodic
        calfig = determine_topology(path_to_calfig; structure=StructureType.MOF, clusterings=[Clustering.Auto])
        # other calls to determine_topology
        determine_topology(path_to_abw; structure=StructureType.Zeolite)
        determine_topology(path_to_abw; structure=StructureType.Cluster, bonding=Bonding.Guess)
        determine_topology(path_to_abw; structure=StructureType.Auto, clusterings=[Clustering.SingleNodes])
        determine_topology(path_to_abw; structure=StructureType.Guess, bonding=Bonding.Auto, clusterings=[Clustering.SingleNodes])
        determine_topology(path_to_abw; clusterings=[Clustering.AllNodes,Clustering.Auto])
        determine_topology(path_to_calfig; bonding=Bonding.Input)
        determine_topology(path_to_calfig; bonding=Bonding.Input, clusterings=[Clustering.PEM])
        CrystalNets.toggle_export()
        CrystalNets.toggle_warning()
        CrystalNets.toggle_export(false)
        CrystalNets.toggle_warning(false)
    end

    rm(tmpdir; recursive=true)
    redirect_stdout(old_stdout)
    redirect_stderr(old_stderr)
end
