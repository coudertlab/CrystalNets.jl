## Handling of the topological archive internally used to recognize topologies.
using Pkg.Artifacts

const CRYSTALNETS_ARCHIVE_VERSION = "2.2"
const arc_location = artifact"archives"

"""
    const CRYSTALNETS_ARCHIVE::Dict{String,String}

The archive used to recognize known topologies.

You probably don't need to access it directly: rely on [`recognize_topology`](@ref) to read
and the various archive functions like [`add_to_current_archive!`](@ref) to write.
"""
const CRYSTALNETS_ARCHIVE = if isdir(arc_location) && !isempty(readdir(arc_location))
    flag, parsed = parse_arcs(arc_location)
    if !flag
        error("""CrystalNets.jl appears to have a broken installation (incompatible archive version).
        Please rebuild CrystalNets.jl with `import Pkg; Pkg.build("CrystalNets")`.
        """)
    end
    parsed
else
    error("""CrystalNets.jl appears to have a broken installation (missing default archive).
          Please rebuild CrystalNets.jl with `import Pkg; Pkg.build("CrystalNets")`.

          This issue may be due to an error while using `CrystalNets.clean_default_archive!(...)`.
          If you have not encountered any such error and you did not modify or erase the archive located at $arc_location please open an issue at https://github.com/coudertlab/CrystalNets.jl/issues/new
          """)
end

"""
    const REVERSE_CRYSTALNETS_ARCHIVE::Dict{String,String}

Reverse of [`CRYSTALNETS_ARCHIVE`](@ref).

Can be used to query the topological genome of known nets, as in:
```jldoctest
julia> REVERSE_CRYSTALNETS_ARCHIVE["dia"]
"3 1 2 0 0 0 1 2 0 0 1 1 2 0 1 0 1 2 1 0 0"

julia> topological_genome(CrystalNet(PeriodicGraph(ans)))
dia
```

!!! note
    It is also possible to directly access the topological genome as a `PeriodicGraph`
    by parsing the name as a [`TopologicalGenome`](@ref):
    ```jldoctest
    julia> PeriodicGraph(parse(TopologicalGenome, "pcu"))
    PeriodicGraph3D(1, PeriodicEdge3D[(1, 1, (0,0,1)), (1, 1, (0,1,0)), (1, 1, (1,0,0))])

    julia> string(PeriodicGraph(parse(TopologicalGenome, "nbo"))) == REVERSE_CRYSTALNETS_ARCHIVE["nbo"]
    true
    ```
"""
const REVERSE_CRYSTALNETS_ARCHIVE = Dict{String,String}(id => (startswith(key, "unstable") ? key[10:end] : key) for (key, id) in CRYSTALNETS_ARCHIVE)

export REVERSE_CRYSTALNETS_ARCHIVE

export clean_default_archive!,
       set_default_archive!,
       empty_default_archive!,
       change_current_archive!,
       refresh_current_archive!,
       add_to_current_archive!

# export make_archive

function _reset_archive!()
    global CRYSTALNETS_ARCHIVE
    global REVERSE_CRYSTALNETS_ARCHIVE
    global arc_location
    empty!(CRYSTALNETS_ARCHIVE)
    merge!(CRYSTALNETS_ARCHIVE, last(parse_arcs(arc_location)))
    empty!(REVERSE_CRYSTALNETS_ARCHIVE)
    merge!(REVERSE_CRYSTALNETS_ARCHIVE, Dict{String,String}(last(x) => first(x) for x in CRYSTALNETS_ARCHIVE))
    nothing
end

function validate_archive(custom_arc, avoid_recompute=true)::Dict{String,String}
    arc = try
        Serialization.deserialize(custom_arc)
    catch
        nothing
    end
    if arc isa Dict{String,String} && !isempty(arc) && isnumeric(first(first(keys(arc))))
        @ifwarn @info "Processing input as a serialized CrystalNets.jl archive"
    else
        try
            flag, parsed = parse_arc(custom_arc)
            if flag && avoid_recompute
                arc = parsed
            else
                @ifwarn begin
                    @info "Processing input as a topological .arc file"
                    if occursin("Made by CrystalNets.jl", readline(custom_arc))
                        @info "This archive was generated by an older version of CrystalNets."
                    end
                    @info "Keys will be converted to the topological genome used by CrystalNets. This may take a while."
                end
                arc_per_thread = [Pair{String,String}[] for _ in 1:nthreads()]
                Threads.@threads for (key, id) in collect(parsed)
                    _g = topological_genome(CrystalNet(PeriodicGraph(key)))
                    if _g.unstable || !isempty(_g.error)
                        if _g.unstable
                            println(stderr, "Net ", id, " is unstable.")
                        else
                            println(stderr, "Failed for net ", id, " with error:", g.error)
                        end
                        continue
                    end
                    genome = string(_g.genome)
                    if !isempty(genome)
                        push!(arc_per_thread[threadid()], (genome => id))
                    end
                end
                arc = Dict{String,String}(pop!(arc_per_thread))
                for _arc in arc_per_thread
                    merge!(arc, Dict{String,String}(_arc))
                end
            end
        catch e
            print(stderr, """
            Impossible to parse input as a topological archive. Please use a format
            similar to that of the RCSR Systre .arc file, with for each entry at
            least the "id" and "key" fields.

            This error may also occur if the given archive was empty. If you
            wish to set an empty archive, use `CrystalNets.empty_default_archive!()`
            """, '\n', "Encountered error: ")
            rethrow()
            # showerror(stderr, e)
            # Base.display_error(stderr, e, catch_backtrace())
        end
    end
    return arc
end

"""
    clean_default_archive!(custom_arc; validate=true, refresh=true)

Erase the default archive used by CrystalNets.jl to recognize known topologies
and replace it with a new one from the file located at `custom_arc`.

The `validate` parameter controls whether the new file is checked and converted
to a format usable by CrystalNets.jl. If unsure, leave it set.

The `refresh` optional parameter controls whether the current archive should be
replaced by the new default one.

!!! warning
    This archive will be kept and used for subsequent runs of CrystalNets.jl, even
    if you restart your Julia session.

    To only change the archive for the current session, use [`change_current_archive!(custom_arc)`](@ref change_current_archive!).

    See also [`refresh_current_archive!`](@ref) for similar uses.

!!! warning
    The previous default archive cannot be recovered afterwards, so make sure to
    keep a copy if necessary. The default archive is the set of ".arc" files located
    at `joinpath(dirname(dirname(pathof(CrystalNets))), "archives")`.
"""
function clean_default_archive!(custom_arc=nothing; validate=true, refresh=true, name="new")
    rm(arc_location; recursive=true)
    mkdir(arc_location)
    if validate
        arc = validate_archive(custom_arc)
        export_arc(joinpath(arc_location, name*".arc"), false, arc)
    else
        cp(custom_arc, arc_location)
    end
    if refresh
        refresh_current_archive!()
    end
    nothing
end

"""
    set_default_archive!()

Set the current archive as the new default archive.

!!! warning
    This archive will be kept and used for subsequent runs of CrystalNets.jl, even
    if you restart your Julia session.
"""
function set_default_archive!(name="new")
    global CRYSTALNETS_ARCHIVE
    export_arc(joinpath(arc_location, name*".arc"))
end

"""
    empty_default_archive!(; refresh=true)

Empty the default archive. This will prevent CrystalNets from recognizing any
topology before they are explicitly added.

The `refresh` optional parameter controls whether the current archive should also
be emptied.

!!! warning
    This empty archive will be kept and used for subsequent runs of CrystalNets.jl, even
    if you restart your Julia session. If you only want to empty the current archive,
    do `empty!(CrystalNets.CRYSTALNETS_ARCHIVE)`.
"""
function empty_default_archive!(; refresh=true)
    global CRYSTALNETS_ARCHIVE
    global REVERSE_CRYSTALNETS_ARCHIVE
    export_arc(arc_location, true)
    if refresh
        empty!(CRYSTALNETS_ARCHIVE)
        empty!(REVERSE_CRYSTALNETS_ARCHIVE)
    end
    nothing
end


function _change_current_archive!(newarc)
    global CRYSTALNETS_ARCHIVE
    global REVERSE_CRYSTALNETS_ARCHIVE
    empty!(CRYSTALNETS_ARCHIVE)
    empty!(REVERSE_CRYSTALNETS_ARCHIVE)
    merge!(CRYSTALNETS_ARCHIVE, newarc)
    merge!(REVERSE_CRYSTALNETS_ARCHIVE,
           Dict{String,String}(last(x) => first(x) for x in CRYSTALNETS_ARCHIVE))
    nothing
end

"""
    change_current_archive!(custom_arc; validate=true)

Erase the current archive used by CrystalNets.jl to recognize known topologies and
replace it with the archive stored in the file located at `custom_arc`.

The `validate` optional parameter controls whether the new file is checked and converted
to a format usable by CrystalNets.jl. If unsure, leave it set.

!!! note
    This modification will only last for the duration of this Julia session.

    If you wish to change the default archive and use it for subsequent runs, use
    [`clean_default_archive!`](@ref).

!!! warning
    Using an invalid archive will make CrystalNets.jl unusable. If this happens,
    simply run [`refresh_current_archive!()`](@ref) to revert to the
    default archive.
"""
function change_current_archive!(custom_arc; validate=true)
    arc::Dict{String,String} = if validate
         validate_archive(custom_arc)
    else
        last(parse_arc(custom_arc))
    end
    _change_current_archive!(arc)
end

"""
    refresh_current_archive!()

Revert the current topological archive to the default one.
"""
function refresh_current_archive!()
    _change_current_archive!(last(parse_arcs(arc_location)))
end


function _update_archive!(id, genome)
    global CRYSTALNETS_ARCHIVE
    global REVERSE_CRYSTALNETS_ARCHIVE
    CRYSTALNETS_ARCHIVE[genome] = id
    REVERSE_CRYSTALNETS_ARCHIVE[id] = genome
    nothing
end

"""
    add_to_current_archive!(id, genome)

Mark `genome` as the topological genome associated with the name `id` in the
current archive.

The input `id` and `genome` are not modified by this operation.

!!! note
    This modification will only last for the duration of this Julia session.

    If you wish to save the archive and use it for subsequent runs, use
    [`set_default_archive!`](@ref) after calling this function.
"""
function add_to_current_archive!(id::AbstractString, genome::AbstractString)
    if !isnumeric(first(genome))
        throw(ArgumentError(lazy"""
            This genome ("$genome") does not look like a genome. Are you sure you did not mix `id` and `genome`?

            If you really want to associate this id with this genome, use `CrystalNets._update_archive!(id, genome)`
            """))
    end
    global CRYSTALNETS_ARCHIVE
    for (x,y) in CRYSTALNETS_ARCHIVE
        if x == genome
            y == id && return
            throw(ArgumentError(lazy"""
                This genome is already registered under the name "$y".

                If you really want to change the name associated with it, use `CrystalNets._update_archive!(id, genome)`
                """))
        end
        if y == id
            throw(ArgumentError(lazy"""
                The name $id already corresponds to a different genome: "$x"

                If you really want to store another genome with the same name, use `CrystalNets._update_archive!(id, genome)`
                """))
        end
    end
    _update_archive!(id, genome)
end

"""
    make_archive(path, destination=nothing, verbose=false)

Make an archive from the files located in the directory given by `path` and export
it to `destination`, if specified. Each file of the directory should correspond
to a unique topology: if a topology is encountered multiple times, it will be assigned
the name of the latest file that bore it.

The archive can then be used with [`change_current_archive!(destination; validate=false)`](@ref change_current_archive!)
for instance.
"""
function make_archive(path, destination, verbose=false)
    arc = Dict{String,String}()
    Threads.@threads for f in readdir(path; join=false)
        name = splitext(f)[1]
        verbose && print("Handling "*name*"... ")
        flag = false
        flagerror = Ref{Any}(Tuple{Vector{Int},String}[])
        results::InterpenetratedTopologyResult = try
            x = topological_genome(UnderlyingNets(parse_chemfile(path*f)))
            verbose && println(name*" done.")
            x
        catch e
            flag = true
            flagerror[] = e
            InterpenetratedTopologyResult()
        end
        for (i, (topology, nfold)) in enumerate(results)
            genome = string(topology)
            if startswith(genome, "unstable") || genome == "non-periodic"
                flag = true
                push!(flagerror[]::Vector{Tuple{Vector{Int},String}}, (vmap, genome))
                continue
            end
            verbose && nfold != 1 && println(nfold, "-fold catenated net found for ", name)
            arc[genome] = length(results) == 1 ? name : (name * '_' * string(i))
        end
        if flag
            e = flagerror[]
            isinterrupt(e) && rethrow()
            if e isa Vector{Tuple{Vector{Int},String}}
                for (vmap, instability) in e
                    println(stderr, "The component of file ", f, " containing atoms ",
                                    vmap, " was found to be ", instability, '.')
                end
            else
                println(stderr, "Failed for file ", f, " with error:")
                showerror(stderr, e)
                Base.display_error(stderr, e, catch_backtrace())
                println(stderr)
            end
        end
    end
    if !(destination isa Nothing)
        export_arc(destination, false, arc)
    end
    return arc
end
