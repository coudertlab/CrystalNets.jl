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

const usable_structuretypes = "mof, zeolite, cluster, guess and auto"
const usable_clusterings = "auto, input, atom, singlenodes, allnodes, standard, pe and pem"
const usable_bondings = "input, guess and auto"

function parse_commandline(args)
    s = ArgParseSettings(prog = "CrystalNets" * (Sys.iswindows() ? ".exe" : ""),
                         description = "Automatic recognition of crystal net topologies.",
                         epilog = """\n\n\n\nSTRUCTURE options:\n\n
                         \ua0\ua0\ua0* mof: consider the input as a MOF. Identify organic and inorganic clusters using heuristics.\n\n
                         \ua0\ua0\ua0* zeolite: attempt to force each O atom two have exactly 2 non-O neighbours.\n\n
                         \ua0\ua0\ua0* cluster: like mof but metallic atoms are not given a larger radius for bond detection.\n\n
                         \ua0\ua0\ua0* guess: discard the input residues and try structure mode "cluster". If it fails, fall back to "auto".\n\n
                         \ua0\ua0\ua0* auto: no specific structure. Default option.\n\n
                         \n\n
                         BONDING options:\n\n
                         \ua0\ua0\ua0* input: use the bonds explicitly given by the input file. Fail if bonds are not provided by the input.\n\n
                         \ua0\ua0\ua0* guess: guess bonds using a variant of chemfiles / VMD algorithm.\n\n
                         \ua0\ua0\ua0* auto: if the input possesses explicit bonds, use them unless they are suspicious. Otherwise, fall back to "guess". Default option.\n\n
                         \n\n
                         CLUSTERING options:\n\n
                         \ua0\ua0\ua0* auto: cluster according to the structure.\n\n
                         \ua0\ua0\ua0* input: use the input residues as clusters. Fail if some atom does not belong to a residue.\n\n
                         \ua0\ua0\ua0* atom: each atom is its own cluster.\n\n
                         \ua0\ua0\ua0* singlenodes: each SBU is clustered into a single vertex.\n\n
                         \ua0\ua0\ua0* allnodes: within each SBU, aromatic cycles are collapsed into vertices.\n\n
                         \ua0\ua0\ua0* standard: each metal is its own vertex.\n\n
                         \ua0\ua0\ua0* pe: each metallic SBU is implicitly represented by its points of extensions.\n\n
                         \ua0\ua0\ua0* pem: each metallic SBU is represented by its points of extension and its metal.\n\n
                         """,
                        #  \n\n
                        #  CREATE_MODE options:\n\n
                        #  \ua0\ua0\ua0* empty: empty archive, unable to recognize any topological structure.\n\n
                        #  \ua0\ua0\ua0* rcsr: RCSR Systre archive.\n\n
                        #  \ua0\ua0\ua0* zeolites: zeolites topologies from the database of zeolite structures.\n\n
                        #  \ua0\ua0\ua0* epinet: systre nets from the EPINET database.\n\n
                        #  \ua0\ua0\ua0* full: combined rcsr, zeolites and epinet archives. Default option.\n\n
                        #  """,
                         preformatted_epilog = false,
                         autofix_names = true,
                         add_help = false,
                         usage = """
                    usage: CrystalNets [-a ARCHIVE_PATH [-u NAME [-f] | -r [-f]]]
                                       [[-s STRUCTURE] [-b BONDING] [-c CLUSTERING] | -k]
                                       [--no-export | -e DIR_PATH] CRYSTAL_FILE
                    """)
                    #       CrystalNets -a ARCHIVE_PATH -n [CREATE_MODE] [-f]              (Form B)
                    #       CrystalNets -a ARCHIVE_PATH -d [-f]                            (Form C)
                    # """)

    # add_arg_group!(s, "Options common to all forms")
    @add_arg_table! s begin
        "--help", "-h"
            help = "Show this help message and exit."
            action = :store_true

        "--no-warn"
            help = "Discard all warnings and information messages."
            action = :store_true

        "--no-export"
            help = "Do not automatically export the parsed input."
            action = :store_true

        "--archive", "-a"
            help = """Specify the path to an archive used to recognize topologies.
            If unspecified while using Form A, defaults to a combination of the RCSR
            Systre archive (available at http://rcsr.net/systre), the known zeolite
            topologies (registered at http://www.iza-structure.org/) and EPINET s-nets
            (available at http://epinet.anu.edu.au).
            """
            metavar = "ARCHIVE_PATH"

        # "--force", "-f"
        #     help = """Force the modification of the archive. Can only be used with
        #     --update-archive, --new-archive, --remove-from-archive or --delete-archive.
        #     """
        #     action = :store_true
    end

    # add_arg_group!(s, "Form A: give the name of the topology of a net")
    @add_arg_table! s begin
        "input"
            help = "Path to the crystal to be analyzed"
            required = false
            metavar = "CRYSTAL_FILE"

        "--structure", "-s"
            help = """Structure mode, to be chosen between $usable_structuretypes.
            See bottom for more details.
            """
            metavar = "STRUCTURE"

        "--clustering", "-c"
            help = """Clustering algorithms, to be chosen between $usable_clusterings.
            See bottom for more details.
            """
            metavar = "CLUSTERING"

        "--bond-detect", "-b"
            help = """Bond detection mode, to be chosen between $usable_bondings.
            See bottom for more details.
            """
            metavar = "BONDING"

        "--export-to", "-e"
            help = """Automatically export the parsed input to the directory at DIR_PATH.
            By default this option is enabled with DIR_PATH=$(tempdir())"""
            metavar = "DIR_PATH"

        "--key", "-k"
            help = """If set, consider the CRYSTAL_FILE parameter as a topological key
            rather than the path to a crystal.
            """
            action = :store_true

        # "--update-archive", "-u"
        #     help = """Give a name to the topology of the crystal and add it to the
        #     archive at ARCHIVE_PATH. If the topology is already known to the archive,
        #     this will do nothing unless --force is passed on. The returned name is
        #     the previous binding, if present, or "UNKNOWN" followed by the topological
        #     genome otherwise.
        #     """
        #     metavar = "NAME"

        # "--remove-from-archive", "-r"
        #     help = """Remove the topology from the archive at ARCHIVE_PATH. If the
        #     topology is absent, this will error unless --force is passed on. The
        #     returned name is the previous binding, if present, or "UNKNOWN" followed
        #     by the topological genome otherwise.
        #     """
        #     action = :store_true
    end

    # add_arg_group!(s, "Form B: create a new archive")
    # @add_arg_table! s begin
    #     "--new-archive", "-n"
    #         help = """Create an archive at ARCHIVE_PATH. CREATE_MODE can be
    #         either empty, full, rcsr, zeolites or epinet.
    #         See bottom for more details.
    #         """
    #         metavar = "CREATE_MODE"
    #         nargs = '?'
    #         constant = "full"
    # end

    # add_arg_group!(s, "Form C: delete an archive")
    # @add_arg_table! s begin
    #     "--delete-archive", "-d"
    #         help = """Delete the archive at ARCHIVE_PATH."""
    #         action = :store_true
    # end

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
        return parse_error(lazy"Unrecognized argument format: $x.")
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

function split_clusterings(s)
    if s == "auto"
        Clustering.Auto
    elseif s == "input"
        Clustering.Input
    elseif s == "atom"
        Clustering.EachVertex
    elseif s == "singlenodes"
        Clustering.SingleNodes
    elseif s == "allnodes"
        Clustering.AllNodes
    elseif s == "standard"
        Clustering.Standard
    elseif s == "pe"
        Clustering.PE
    elseif s == "pem"
        Clustering.PEM
    else
        return parse_error(lazy"Unknown clustering mode: $s. Choose between $usable_clusterings.")
    end
end

"""
    main(ARGS)

Function called when using the module as an executable.

Return code can be:
* 0: no error
* 1: the chemical bond system has no periodicity
* 2: invalid input
* 3: parsing error
* 4: internal CrystalNets.jl error
* 5: unhandled CrystalNets.jl error, please report
"""
function main(args)
    try
        _parsed_args = parse_commandline(args)
        if _parsed_args isa Int
            return _parsed_args
        end
        parsed_args::Dict{Symbol,Any} = _parsed_args

        @toggleassert !parsed_args[:help]

        # force::Bool = parsed_args[:force]

        if parsed_args[:no_warn]::Bool
            toggle_warning(false)
        end
        if parsed_args[:no_export]::Bool
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

        # @parse_to_str_or_nothing new_archive new_archive_mode

        # if new_archive_mode isa String
        #     if parsed_args[:delete_archive]::Bool
        #         return parse_error("Cannot execute both forms B and C.")
        #     elseif !isnothing(parsed_args[:input])
        #         return parse_error("Cannot execute both forms A and B.")
        #     end
        #     if archive isa Nothing
        #         return parse_error("Cannot create an archive without specifying its location: use the --archive option to provide a path for the archive.")
        #     else
        #         if force && isfile(archive)
        #             parse_error("The path specified by --archive already exists. Use --force to remove the existing file and replace it.")
        #         end
        #         if new_archive_mode == "full"
        #             export_arc(archive, false)
        #         elseif new_archive_mode == "empty"
        #             export_arc(archive, true)
        #         elseif new_archive_mode*".arc" ∈ readdir(arc_location)
        #             flag, arc = parse_arc(arc_location * new_archive_mode *".arc")
        #             if !flag
        #                 internal_error("""CrystalNets.jl appears to have a broken installation (the archive version is older than that package's).
        #                 Please rebuild CrystalNets.jl with `import Pkg; Pkg.build("CrystalNets")`.
        #                 """, AssertionError("flag"), catch_backtrace())
        #             end
        #             export_arc(archive, false, arc)
        #         else
        #             return parse_error("""Unknown archive: $new_archive_mode. Choose between "full", "empty", "rcsr", "zeolites" or "epinet".""")
        #         end
        #     end
        #     return 0
        # end

        # if parsed_args[:delete_archive]::Bool
        #     if !isnothing(parsed_args[:input])
        #         return parse_error("Cannot execute both forms A and C.")
        #     end
        #     if archive isa Nothing
        #         return parse_error("""Cannot delete an archive without specifying its location: use the --archive option to provide a path for the archive.""")
        #     else
        #         if force
        #             if !isfile(archive)
        #                 ifwarn("The specified archive does not exist.")
        #             end
        #             try
        #                 run(`rm -f $archive`)
        #             catch e
        #                 return internal_error("Error encountered while deleting the archive:",
        #                                       e, catch_backtrace())
        #             end
        #         else
        #             if !isfile(archive)
        #                 return parse_error("The specified archive does not exist.")
        #             end
        #             try
        #                 run(`rm $archive`)
        #             catch e
        #                 return internal_error("""The following error was encountered while deleting the archive. Try with option --force.""",
        #                                       e, catch_backtrace())
        #             end
        #         end
        #     end
        #     return 0
        # end

        @parse_to_str_or_nothing structure structure_mode
        structure::StructureType._StructureType = begin
            if structure_mode isa Nothing
                StructureType.Auto
            else
                if structure_mode == "mof"
                    StructureType.MOF
                elseif structure_mode == "zeolite"
                    StructureType.Zeolite
                elseif structure_mode == "cluster"
                    StructureType.Cluster
                elseif structure_mode == "guess"
                    StructureType.Guess
                elseif structure_mode == "auto"
                    StructureType.Auto
                else
                    return parse_error(lazy"Unknown structure type: $structure_mode. Choose between $usable_structuretypes.")
                end
            end
        end

        @parse_to_str_or_nothing clustering clustering_mode
        clusterings::Vector{Clustering._Clustering} = begin
            if clustering_mode isa Nothing
                [Clustering.Auto]
            else
                clustering_splits = split(clustering_mode, ',')
                _clusterings = Vector{Clustering._Clustering}(undef, length(clustering_splits))
                for (i,s) in enumerate(clustering_splits)
                    _clust = split_clusterings(s)
                    if _clust isa Clustering._Clustering
                        _clusterings[i] = _clust
                    else
                        return _clust
                    end
                end
                _clusterings
            end
        end

        @parse_to_str_or_nothing bond_detect
        bonding::Bonding._Bonding = begin
            if bond_detect isa Nothing
                Bonding.Auto
            else
                if bond_detect == "input"
                    Bonding.Input
                elseif bond_detect == "guess"
                    Bonding.Guess
                elseif bond_detect == "auto"
                    Bonding.Auto
                else
                    return parse_error(lazy"Unknown bond detection mode: $bond_detect. Choose between $usable_bondings.")
                end
            end
        end

        # @parse_to_str_or_nothing update_archive new_topology_name
        # if new_topology_name isa String && isnothing(archive)
        #     return parse_error("""Cannot update an archive without specifying its location: use the --archive option to provide a path for the archive.""")
        # end

        # remove_from_archive::Bool = parsed_args[:remove_from_archive]
        # if remove_from_archive && isnothing(archive)
        #     return parse_error("""Cannot modify an archive without specifying its location: use the --archive option to provide a path for the archive.""")
        # end
        # if remove_from_archive && new_topology_name isa String
        #     return parse_error("""Cannot both add and remove a topology from the archive. Choose only one of the options --update-archive and --remove-from-archive.""")
        # end

        iskey::Bool = parsed_args[:key]
        if iskey && (structure_mode isa String || clustering_mode isa String)
            msg = structure_mode isa String ? "structure type" : "clustering mode"
            return parse_error(lazy"Cannot consider the input as a topological key while also using a specified $(msg) because keys miss atom type information.")
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

        unets::UnderlyingNets = try
            if iskey
                g = try
                    PeriodicGraph(input_file)
                catch e
                    return invalid_input_error("""Impossible to parse the given topological key because of the following error:""",
                                            e, catch_backtrace())
                end
                UnderlyingNets(g)
            else
                crystal::Crystal = try
                    parse_chemfile(input_file, Options(;export_input=export_path,
                                                        export_net=export_path,
                                                        structure,
                                                        bonding,
                                                        clusterings,
                                                        throw_error=true,
                                                      ))
                catch e
                    return invalid_input_error("""The input file could not be correctly parsed as as a crystal because of the following error:""",
                                            e, catch_backtrace())
                end
                UnderlyingNets(crystal)
            end
        catch e
            return invalid_input_error("""The input cannot be analyzed because of the following error:""",
                                        e, catch_backtrace())
        end
        genomes::InterpenetratedTopologyResult = try
            topological_genome(unets)
        catch e
            return internal_error("""Internal error encountered while computing the topological genome:""",
                                    e, catch_backtrace())
        end

        #=
        if new_topology_name isa String
            if force
                _update_archive!(new_topology_name, genomes)
                export_arc(archive, false)
            else
                flag = true
                try
                    add_to_current_archive!(new_topology_name, genomes)
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
                delete!(CRYSTALNETS_ARCHIVE_VERSION, genome)
                export_arc(archive, false)
            else
                return parse_error("""The genome "$genome" was unknown to the archive and thus could not be deleted. Use option --force to discard this error.""")
            end
        end
        =#

        if length(genomes) == 0
            println(genomes)
            return 1
        end

        println(genomes)
        return 0
    catch e
        return unhandled_error("CrystalNets encountered an unhandled exception:",
                               e, catch_backtrace())
    end

end


function julia_main()::Cint
    main(ARGS)
end
