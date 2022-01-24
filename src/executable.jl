## Entry point of the executable

using ArgParse

function invalid_input_error(msg, e, stacktrace)
    println(stderr, msg)
    Base.display_error(stderr, e, stacktrace)
    return 2
end

function parse_error(msg)
    println(stderr, msg)
    println(stderr, "See --help for usage.")
    return 3
end

function internal_error(msg, e, stacktrace)
    println(stderr, msg)
    Base.display_error(stderr, e, stacktrace)
    return 4
end

function unhandled_error(msg, e, stacktrace)
    println(stderr, msg)
    Base.display_error(stderr, e, stacktrace)
    return 5
end

function parse_commandline(args)
    s = ArgParseSettings(prog = "CrystalNets" * (Sys.iswindows() ? ".exe" : ""),
                         description = "Automatic recognition of crystal net topologies.",
                         epilog = """\n\n\n\nCLUSTERING_MODE options:\n\n
                         \ua0\ua0\ua0* input: use the input residues as clusters. Fail if some atom does not belong to a residue.\n\n
                         \ua0\ua0\ua0* atom: each atom is its own cluster.\n\n
                         \ua0\ua0\ua0* mof: discard the input residues and consider the input as a MOF. Identify organic and inorganic clusters using a simple heuristic based on the atom types.\n\n
                         \ua0\ua0\ua0* guess: discard the input residues and try clustering mode "mof". If it fails, fall back to "atom".\n\n
                         \ua0\ua0\ua0* auto: if the input assigns each atom to a residue, equivalent to "input". Otherwise, fall back to "guess". Default option.\n\n
                         \n\n
                         BONDING_MODE options:\n\n
                         \ua0\ua0\ua0* input: use the bonds explicitly given by the input file. Fail if bonds are not provided by the input.\n\n
                         \ua0\ua0\ua0* guess: guess bonds using a variant of chemfiles / VMD algorithm.\n\n
                         \ua0\ua0\ua0* auto: if the input possesses explicit bonds, use them unless they are suspicious. Otherwise, fall back to "guess". Default option.\n\n
                         \n\n
                         CREATE_MODE options:\n\n
                         \ua0\ua0\ua0* empty: empty archive, unable to recognize any topological structure.\n\n
                         \ua0\ua0\ua0* rcsr: RCSR Systre archive.\n\n
                         \ua0\ua0\ua0* zeolites: zeolites topologies from the database of zeolite structures.\n\n
                         \ua0\ua0\ua0* full: combined rcsr and zeolites archives. Default option.\n\n
                         """,
                         preformatted_epilog = false,
                         autofix_names = true,
                         add_help = false,
                         usage = """
                   usage: CrystalNets [-a ARCHIVE_PATH [-u NAME [-f] | -r [-f]]]
                                      [[-c CLUSTERING_MODE] [-b BONDING_MODE] | -g]
                                      [--no-export | -e DIR_PATH] CRYSTAL_FILE   (Form A)
                          CrystalNets -a ARCHIVE_PATH -n [CREATE_MODE] [-f]      (Form B)
                          CrystalNets -a ARCHIVE_PATH -d [-f]                    (Form C)
                   """)

    add_arg_group!(s, "Options common to all forms")
    @add_arg_table! s begin
        "--help", "-h"
            help = "Show this help message and exit."
            action = :store_true

        "--no-warn"
            help = "Discard all warnings and information messages."
            action = :store_true

        "--archive", "-a"
            help = """Specify the path to an archive used to recognize topologies.
            If unspecified while using Form A, defaults to a combination of the RCSR
            Systre archive (available at http://rcsr.net/systre) and the known zeolite
            topologies (registered at http://www.iza-structure.org/).
            """
            metavar = "ARCHIVE_PATH"

        "--force", "-f"
            help = """Force the modification of the archive. Can only be used with
            --update-archive, --new-archive, --remove-from-archive or --delete-archive.
            """
            action = :store_true
    end

    add_arg_group!(s, "Form A: give the name of the topology of a net")
    @add_arg_table! s begin
        "input"
            help = "Path to the crystal to be analyzed"
            required = false
            metavar = "CRYSTAL_FILE"

        "--clustering", "-c"
            help = """Clustering mode, to be chosen between input, atom, mof, guess and auto.
            See bottom for more details.
            """
            metavar = "CLUSTERING_MODE"

        "--bond-detect", "-b"
            help = """Bond detection mode, to be chosen between input, guess and auto.
            See bottom for more details.
            """
            metavar = "BONDING_MODE"

        "--no-export"
            help = "Do not automatically export the parsed input."
            action = :store_true

        "--export-to", "-e"
            help = """Automatically export the parsed input to the directory at DIR_PATH.
            By default this option is enabled with DIR_PATH=$(tempdir())"""
            metavar = "DIR_PATH"

        "--genome", "-g"
            help = """If set, consider the CRYSTAL_FILE parameter as a topological
            genome rather than the path to a crystal.
            """
            action = :store_true

        "--update-archive", "-u"
            help = """Give a name to the topology of the crystal and add it to the
            archive at ARCHIVE_PATH. If the topology is already known to the archive,
            this will do nothing unless --force is passed on. The returned name is
            the previous binding, if present, or "UNKNOWN" followed by the topological
            genome otherwise.
            """
            metavar = "NAME"

        "--remove-from-archive", "-r"
            help = """Remove the topology from the archive at ARCHIVE_PATH. If the
            topology is absent, this will error unless --force is passed on. The
            returned name is the previous binding, if present, or "UNKNOWN" followed
            by the topological genome otherwise.
            """
            action = :store_true
    end

    add_arg_group!(s, "Form B: create a new archive")
    @add_arg_table! s begin
        "--new-archive", "-n"
            help = """Create an archive at ARCHIVE_PATH. CREATE_MODE can be
            either empty, full, rcsr or zeolites.
            See bottom for more details.
            """
            metavar = "CREATE_MODE"
            nargs = '?'
            constant = "full"
    end

    add_arg_group!(s, "Form C: delete an archive")
    @add_arg_table! s begin
        "--delete-archive", "-d"
            help = """Delete the archive at ARCHIVE_PATH."""
            action = :store_true
    end

    ret = try
        parse_args(args, s; as_symbols=true)
    catch e
        return internal_error("An error happened while parsing the input argument and options:",
                              e, catch_backtrace())
    end
    if ret[:help]
        ArgParse.show_help(s; exit_when_done=false)
        return 0
    end

    return ret
end


function parse_to_str_or_nothing(@nospecialize(x))::Union{Nothing,String,Int}
    if isnothing(x)
        return nothing
    end
    str = try
        string(x)
    catch e
        return parse_error("Unrecognized argument format: $x.")
    end
    return String(strip(str))
end
macro parse_to_str_or_nothing(x, name=x)
    return quote
        $(esc(name)) = parse_to_str_or_nothing($(esc(:parsed_args))[$(esc(QuoteNode(x)))])
        $(esc(name)) isa $(esc(Int)) && return $(esc(name))
        $(esc(name))::$(esc(Union{Nothing,String}))
    end
end

main(x::String) = main(split(x))

function main(args)
    try
        parsed_args = parse_commandline(args)
        if parsed_args isa Int
            return parsed_args
        end
        parsed_args::Dict{Symbol,Any}

        @toggleassert !parsed_args[:help]

        force::Bool = parsed_args[:force]

        if parsed_args[:no_warn]
            toggle_warning(false)
        end
        if parsed_args[:no_export]
            if !isnothing(parsed_args[:export_to])
                return parse_error("""Cannot use both arguments "--export-to" and "--no-export".""")
            end
            toggle_export(false)
        end

        @parse_to_str_or_nothing export_to
        export_path::String = if export_to isa Nothing
            DOEXPORT[] ? tempdir() : ""
        else
            toggle_export(true)
            export_to
        end

        @parse_to_str_or_nothing archive

        @parse_to_str_or_nothing new_archive new_archive_mode

        if new_archive_mode isa String
            if parsed_args[:delete_archive]
                return parse_error("Cannot execute both forms B and C.")
            elseif !isnothing(parsed_args[:input])
                return parse_error("Cannot execute both forms A and B.")
            end
            if archive isa Nothing
                return parse_error("Cannot create an archive without specifying its location: use the --archive option to provide a path for the archive.")
            else
                if force && isfile(archive)
                    parse_error("The path specified by --archive already exists. Use --force to remove the existing file and replace it.")
                end
                if new_archive_mode == "full"
                    export_arc(archive, false)
                elseif new_archive_mode == "empty"
                    export_arc(archive, true)
                elseif new_archive_mode*".arc" ∈ readdir(arc_location)
                    flag, arc = parse_arc(arc_location * new_archive_mode *".arc")
                    if !flag
                        internal_error("""CrystalNets.jl appears to have a broken installation (the archive version is older than that package's).
                        Please rebuild CrystalNets.jl with `import Pkg; Pkg.build("CrystalNets")`.
                        """, AssertionError("flag"), catch_backtrace())
                    end
                    export_arc(archive, false, arc)
                else
                    return parse_error("""Unknown archive: $new_archive_mode. Choose between "full", "empty", "rcsr" or "zeolites".""")
                end
            end
            return 0
        end

        if parsed_args[:delete_archive]
            if !isnothing(parsed_args[:input])
                return parse_error("Cannot execute both forms A and C.")
            end
            if archive isa Nothing
                return parse_error("""Cannot delete an archive without specifying its location: use the --archive option to provide a path for the archive.""")
            else
                if force
                    if !isfile(archive)
                        ifwarn("The specified archive does not exist.")
                    end
                    try
                        run(`rm -f $archive`)
                    catch e
                        return internal_error("Error encountered while deleting the archive:",
                                              e, catch_backtrace())
                    end
                else
                    if !isfile(archive)
                        return parse_error("The specified archive does not exist.")
                    end
                    try
                        run(`rm $archive`)
                    catch e
                        return internal_error("""The following error was encountered while deleting the archive. Try with option --force.""",
                                              e, catch_backtrace())
                    end
                end
            end
            return 0
        end

        @parse_to_str_or_nothing clustering cluster_mode
        clustering::ClusteringMode._ClusteringMode = begin
            if cluster_mode isa Nothing
                ClusteringMode.Auto
            else
                if cluster_mode == "input"
                    ClusteringMode.Input
                elseif cluster_mode == "atom"
                    ClusteringMode.EachVertex
                elseif cluster_mode == "mof"
                    ClusteringMode.MOF
                elseif cluster_mode == "zeolite"
                    ClusteringMode.Zeolite
                elseif cluster_mode == "cluster"
                    ClusteringMode.Cluster
                elseif cluster_mode == "guess"
                    ClusteringMode.Guess
                elseif cluster_mode == "auto"
                    ClusteringMode.Auto
                else
                    return parse_error("""Unknown clustering mode: $cluster_mode. Choose between "input", "atom", "mof", "guess" and "auto".""")
                end
            end
        end

        @parse_to_str_or_nothing bond_detect
        bondingmode::BondingMode._BondingMode = begin
            if bond_detect isa Nothing
                BondingMode.Auto
            else
                if bond_detect == "input"
                    BondingMode.Input
                elseif bond_detect == "guess"
                    BondingMode.Guess
                elseif bond_detect == "auto"
                    BondingMode.Auto
                else
                    return parse_error("""Unknown bond detection mode: $bond_detect. Choose between "input", "guess" and "auto".""")
                end
            end
        end

        @parse_to_str_or_nothing update_archive new_topology_name
        if new_topology_name isa String && isnothing(archive)
            return parse_error("""Cannot update an archive without specifying its location: use the --archive option to provide a path for the archive.""")
        end

        remove_from_archive::Bool = parsed_args[:remove_from_archive]
        if remove_from_archive && isnothing(archive)
            return parse_error("""Cannot modify an archive without specifying its location: use the --archive option to provide a path for the archive.""")
        end
        if remove_from_archive && new_topology_name isa String
            return parse_error("""Cannot both add and remove a topology from the archive. Choose only one of the options --update-archive and --remove-from-archive.""")
        end

        isgenome::Bool = parsed_args[:genome]
        if isgenome && cluster_mode isa String
            return parse_error("""Cannot consider the input as a topological genome while also using a specified clustering method on its vertices because genomes miss atom type information.""")
        end

        @parse_to_str_or_nothing input input_file
        if input_file isa Nothing
            return parse_error("Missing a CRYSTAL_FILE.")
        end
        input_file::String

        if archive isa String
            try
                change_current_archive!(archive)
            catch e
                invalid_input_error("""Cannot use the specified archive because of the following error:""",
                                    e, catch_backtrace())
            end
        end

        genome::String = if isgenome
            g = try
                PeriodicGraph(input_file)
            catch e
                return invalid_input_error("""Impossible to parse the given topological genome because of the following error:""",
                                           e, catch_backtrace())
            end
            try
                string(topological_genome(g))
            catch e
                return internal_error("""Internal error encountered while computing the topological genome:""",
                                      e, catch_backtrace())
            end
        else
            crystal::Crystal = try
                parse_chemfile(input_file, Options(; export_input=export_path,
                                                     export_clusters=export_path,
                                                     bonding_mode, clustering))
            catch e
                return invalid_input_error("""The input file could not be correctly parsed as as a crystal because of the following error:""",
                                           e, catch_backtrace())
            end
            net::CrystalNet = try
                clusters, _net = do_clustering(crystal)
                if !isempty(clusters) && !isempty(export_path)
                    try
                        Crystal{Clusters}(crystal, clusters)
                    catch e
                        return internal_error("""Internal error while exporting the clustering of vertices""",
                                                e, catch_backtrace())
                    end
                end
                _net
            catch e
                return invalid_input_error("""The given crystal cannot be analyzed because of the following error:""",
                                           e, catch_backtrace())
            end
            if isempty(net.pos)
                println(stderr, """Error: cannot interpret the input as a crystal or a crystalline framework. Make sure the input is connected and periodic.""")
                return 4
            end
            try
                string(topological_genome(net))
            catch e
                return internal_error("""Internal error encountered while computing the topological genome:""",
                                      e, catch_backtrace())
            end
        end

        id = recognize_topology(genome)

        if new_topology_name isa String
            if force
                _update_archive!(new_topology_name, genome)
                export_arc(archive, false)
            else
                flag = true
                try
                    add_to_current_archive!(new_topology_name, genome)
                catch
                    @ifwarn @info """The archive was not updated because either the name or the genome is already present."""
                    flag = false
                end
                if flag
                    export_arc(archive, false)
                end
            end
        end

        if remove_from_archive
            if force || id isa String
                delete!(CRYSTAL_NETS_VERSION, genome)
                export_arc(archive, false)
            else
                return parse_error("""The genome "$genome" was unknown to the archive and thus could not be deleted. Use option --force to discard this error.""")
            end
        end

        println(id)
        if startswith(id, "UNKNOWN")
            return 1
        end
        return 0

    catch e
        return unhandled_error("CrystalNets encountered an unhandled exception:",
                               e, catch_backtrace())
    end

end


Base.@ccallable function julia_main()::Cint
    main(ARGS)
end
