## CrystalNets main API

"""
    topological_genome(net::CrystalNet{D,T})::String where {D,T}

Compute the topological genome of a net. The topological genome is an invariant
if the net, meaning that it does not depend on its representation. It is also
the string representation of a D-periodic graph such that
`PeriodicGraph{D}(topological_genome(net))` is isomorphic to `net.pge.g` (except
possibly if the `ignore_types` option is unset).

!!! info
    Options must be passed directly within `net`.
"""
function topological_genome(net::CrystalNet{D,T})::TopologicalGenome where {D,T}
    isempty(net.types) && return TopologicalGenome(net.options.error)
    if net.options.ignore_types
        net = CrystalNet{D,T}(net.pge, fill(Symbol(""), length(net.types)), net.options)
    end
    shrunk_net, _collisions = collision_nodes(net)
    _collisions isa Nothing && return TopologicalGenome(net.pge.g, nothing, true)
    collisions::Vector{CollisionNode} = _collisions

    if !net.options.skip_minimize
        flag = true
        try
            shrunk_net, _collisions = minimize(shrunk_net, collisions)
            flag = false
        catch e
            isinterrupt(e) && rethrow()
            if T == Rational{BigInt} || !isoverfloworinexact(e)
                net.options.throw_error && rethrow()
                return TopologicalGenome(string(e)::String)
            end
        end
        if flag # not in the catch to avoid a StackOverflow of errors in case something goes wrong
            newnet = CrystalNet{D,widen(T)}(net; ignore_types=false)
            return topological_genome(newnet)
        end
        _collisions isa Nothing && return TopologicalGenome(shrunk_net.pge.g, nothing, true)
        collisions = _collisions
    end

    export_net = isempty(net.options.export_net) ? isempty(net.options._pos) ?
                 net.options.export_subnets : "" : net.options.export_net
    export_default(net, "net", net.options.name, export_net)

    return topological_genome(shrunk_net, collisions)
end

topological_genome(net::CrystalNet{0,T}) where {T} = TopologicalGenome(net.options.error)

function topological_genome(net::CrystalNet{D,T}, collisions::Vector{CollisionNode})::TopologicalGenome where {D,T}
    try
        g::PeriodicGraph{D} = topological_key(net, collisions)
        unstable = g.width[] == -2
        unstable && (g.width[] = -1)
        ne(g) == 0 && return TopologicalGenome(net.options.error)
        return TopologicalGenome(g, recognize_topology(g), unstable)
    catch e
        isinterrupt(e) && rethrow()
        if T == Rational{BigInt} || !isoverfloworinexact(e)
            net.options.throw_error && rethrow()
            return TopologicalGenome(string(e)::String)
        end
    end
    return topological_genome(CrystalNet{D,widen(T)}(net), collisions)
end

"""
    topological_genome(g::Union{String,PeriodicGraph}, options::Options=Options())
    topological_genome(g::Union{String,PeriodicGraph}; kwargs...)

Compute the topological genome of a periodic graph.
If given a topological key (as a string), it is converted to a `PeriodicGraph` first.
"""
function topological_genome(g::PeriodicGraph, options::Options)
    nets = UnderlyingNets(g, options)
    return topological_genome(nets)
end
topological_genome(s::String, options::Options) = topological_genome(PeriodicGraph(s), options)
topological_genome(g::Union{String,PeriodicGraph}; kwargs...) = topological_genome(g, Options(; throw_error=true, kwargs...))


function _loop_group!(ex::Expr, id::Symbol, net::Symbol, group::Expr)
    for j in 1:length(ex.args)
        arg = ex.args[j]
        if arg isa Symbol
            if arg === :id
                ex.args[j] = id
            elseif arg === :net
                ex.args[j] = net
            elseif arg === :group
                ex.args[j] = group
            end
        elseif arg isa Expr
            _loop_group!(arg, id, net, group)
        end
    end
    nothing
end
macro loop_group(ex)
    @toggleassert ex.head === :for
    ret = Expr(:block)
    for i in 1:3
        id = Symbol("id", i)
        net = Symbol("net", i)
        D = Symbol("D", i)
        group = :(group.$D)
        newex = quote
            if !isempty(group)
                currallowed = first(group)[2][1].options.dimensions
                if isempty(currallowed) || $i in currallowed
                    $(deepcopy(ex))
                end
            end
        end
        _loop_group!(newex, id, net, group)
        push!(ret.args, esc(newex))
    end
    return ret
end

"""
    topological_genome(group::UnderlyingNets)::Vector{Tuple{Vector{Int},String}}

Compute the topological genome of each subnet stored in `group`.

!!! info
    Options must be passed directly within the subnets.
"""
function topological_genome(group::UnderlyingNets)
    ret = Tuple{Vector{Int},TopologyResult}[]
    @loop_group for (id, net) in group
        encountered = Dict{PeriodicGraph,_Clustering}()
        subret = Vector{Tuple{_Clustering,Union{_Clustering,TopologicalGenome}}}(undef, length(net))
        for (j, subnet) in enumerate(net)
            clust = only(subnet.options.clusterings)
            refclust = get!(encountered, subnet.pge.g, clust)
            subret[j] = (clust, refclust == clust ? topological_genome(subnet) : refclust)
        end
        push!(ret, (id, TopologyResult(subret)))
    end
    return ret
end

"""
    recognize_topology(g::PeriodicGraph, arc=CRYSTAL_NETS_ARCHIVE)
    recognize_topology(genome::AbstractString, arc=CRYSTAL_NETS_ARCHIVE)

Attempt to recognize a topological genome from an archive of known genomes.

!!! warning
    This function does a simple lookup, not any kind of topology computation.
    To identify the topology of a `PeriodicGraph` or a `CrystalNet` `x`, query
    `topological_genome(x)` instead.
"""
function recognize_topology(genome::AbstractString, arc=CRYSTAL_NETS_ARCHIVE)
    get(arc, genome, nothing)
end
function recognize_topology(genome::PeriodicGraph, arc=CRYSTAL_NETS_ARCHIVE)
    recognize_topology(string(genome), arc)
end

"""
    determine_topology(path, options::Options)
    determine_topology(path; kwargs...)

Compute the topology of the structure described in the file located at `path`.
This is essentially equivalent to calling
`topological_genome(UnderlyingNets(parse_chemfile(path, options)))`.

In the case where the structure is not made of interpenetrating nets, return the topology
of the only net.
"""
function determine_topology(path, options::Options)
    genomes::Vector{Tuple{Vector{Int},TopologyResult}} =
        topological_genome(UnderlyingNets(parse_chemfile(path, options)))
    if length(genomes) == 1
        return genomes[1][2]
    end
    length(genomes) == 0 && return TopologyResult("")
    return genomes
end
determine_topology(path; kwargs...) = determine_topology(path, Options(; kwargs...))


macro ifvalidgenomereturn(opts, msg, skipcrystal=false)
    crystaldef = skipcrystal ? nothing : :(crystal = parse_chemfile(path, $opts))
    msg = lazy"(found by $msg)"
    ifprintinfo = isempty(msg) ? nothing : :(CrystalNets.@ifwarn @info $msg)

    return esc(quote
    $crystaldef
    group = UnderlyingNets(crystal)
    dim, subnets = isempty(group.D3) ? isempty(group.D2) ? (1, group.D1) : (2, group.D2) : (3, group.D3)
    for (_, nets) in subnets
        for net in nets
            net isa CrystalNet || continue
            sig = get(encountered_graphs, net.pge.g, net.pge.g)
            if haskey(encountered_genomes, net.pge.g)
                unstable, counter = encountered_genomes[sig]
                encountered_genomes[sig] = (unstable, counter + 1)
            else
                genome = topological_genome(net)::TopologicalGenome
                if dim == maxdim && genome.name !== nothing
                    $ifprintinfo
                    return genome
                end
                encountered_graphs[sig] = genome.genome
                encountered_genomes[genome.genome] = (genome.unstable, 1)
            end
        end
    end
    end)
end


"""
    guess_topology(path, options::Options)
    guess_topology(path; kwargs...)

Tries to determine the topology of the file at `path` by passing various options
(starting from the provided options if any) until finding a known topology.
If none is found, return the topological genome encountered most often through
the variation of options.
"""
function guess_topology(path, defopts::Options)
    maxdim = maximum(defopts.dimensions; init=0)
    encountered_graphs = Dict{PeriodicGraph,PeriodicGraph}()
    encountered_genomes = Dict{PeriodicGraph,Tuple{Bool,Int}}()
    if defopts.name == "unnamed"
        defopts = Options(defopts; name=splitext(splitdir(path)[2])[1])
    end

    @ifvalidgenomereturn defopts "using structure=StructureType.Guess" false

    # dryrun = Dict{Symbol,Union{Nothing,Set{Symbol}}}()
    # crystal = parse_chemfile(path, Options(defopts; dryrun))
    atoms = Set{Symbol}(crystal.types)

    if any(x -> (at = atomic_numbers[x]; (ismetal[at] | ismetalloid[at])), atoms)
        @ifvalidgenomereturn Options(defopts; structure=StructureType.MOF) "using MOF clusters"
    end

    defopts = Options(; structure=StructureType.Auto)
    # if defopts.cutoff_coeff == 0.75 # if default cutoff was used, try larger
    #     @ifvalidgenomereturn Options(defopts; cutoff_coeff=0.85) "using longer cutoff"
    # end
    # invalidatoms = union(get(dryrun, :invalidatoms, Set{Symbol}()),
    #                      get(dryrun, :suspect_smallbonds, Set{Symbol}()))
    # hashydrogens = :H ∈ atoms || :D ∈ atoms
    # if hashydrogens
    #     @ifvalidgenomereturn Options(defopts; ignore_atoms=(:H,:D)) "ignoring H"
    #     # delete!(invalidatoms, :H)
    #     # delete!(invalidatoms, :D)
    # end
    # for a in invalidatoms
    #     @ifvalidgenomereturn Options(defopts; ignore_atoms=tuple(a)) "ignoring $a"
    #     if hashydrogens
    #         @ifvalidgenomereturn Options(defopts; ignore_atoms=(a, :H, :D)) "ignoring H and $a"
    #     end
    # end
    # if any(==(:C), atoms)# && :C ∉ invalidatoms # organic solvent
    #     @ifvalidgenomereturn Options(defopts; ignore_atoms=(:C,)) "ignoring C"
    # end
    # if haskey(dryrun, :try_Input_bonds)
    #     @ifvalidgenomereturn Options(defopts; bonding=Bonding.Input) "using input bonds"
    # end
    # @ifvalidgenomereturn Options(defopts; ignore_low_occupancy=true) "removing atoms with occupancy < 0.5"
    if :Al ∈ atoms && (:P ∈ atoms || :Pc ∈ atoms) # ALPO datastructure
        # flag = false
        # for i in vertices(crystal.pge.g)
        #     crystal.types[i] === :O || continue
        #     currt = :O
        #     for x in neighbors(crystal.pge.g, i)
        #         t = crystal.types[x.v]
        #         if (t === :P) | (t === :Al)
        #             if currt === :O
        #                 currt = t
        #             elseif currt === t
        #                 flag = true
        #                 break
        #             end
        #         end
        #     end
        # end
        # if flag
            @ifvalidgenomereturn Options(defopts; ignore_homoatomic_bonds=(:Al,:P)) "removing Al-O-Al and P-O-P bonds"
        # end
        # if :Sn ∈ atoms # observed in some cases for ALPO databases
        #     @ifvalidgenomereturn Options(defopts; ignore_atoms=(:Sn,)) "ignoring Sn"
        # end
    end
    # if haskey(dryrun, :try_no_Auto_bonds)
    #     @ifvalidgenomereturn Options(defopts; bonding=Bonding.Guess) "discarding input bonds and guessing them"
    # end
    # if haskey(dryrun, :collisions)
    #     @ifvalidgenomereturn Options(defopts; authorize_pruning=false) "retaining all colliding atoms"
    # end

    # Finally, if everything fails, return the one encountered most
    maxencounter = 0
    most_plausible_genome::PeriodicGraph = PeriodicGraph{0}()
    final_unstable = false
    for (genome, (unstable, i)) in encountered_genomes
        if i > maxencounter
            maxencounter = i
            most_plausible_genome = genome
            final_unstable = unstable
        end
    end
    return most_plausible_genome == PeriodicGraph{0}() ? begin
            issetequal(defopts.dimensions, [1,2,3]) ? TopologicalGenome() :
                TopologicalGenome("could not guess (perhaps because no net with suitable dimension in $(defopts.dimensions) was found?)")
            end : TopologicalGenome(most_plausible_genome, nothing, final_unstable)
end
guess_topology(path; kwargs...) = guess_topology(path, Options(structure=StructureType.Guess; kwargs...))


"""
    determine_topology_dataset(path, save, autoclean, options::Options)
    determine_topology_dataset(path, save=true, autoclean=true; kwargs...)

Given a path to a directory containing structure input files, compute the
topology of each structure within the directory.
Return a dictionary linking each file name to the result.
The result is a [`TopologyResult`](@ref), containing the topological genome, the name if
known and the stability of the net. In case of error, the exception is reported.

It is strongly recommended to toggle warnings off (through [`toggle_warning`](@ref)) and
not to export any file since those actions may critically reduce performance,
especially for numerous files.

If `save` is set, the result is also stored in a julia serialized file located at
"\\\$path/../results_\\\$i" where `i` is the lowest integer such that this path does
not already exist at the start of the computation.
While processing, this path will be used to create a directory storing the
current state of the computation: to continue an interrupted computation, simply
pass this temporary directory as the path. If `autoclean` is set, this directory
is removed at the end if the computation was successful.

If `save` is set and `autoclean` is unset, the directory of temporary files will
be renamed into "\\\$path/../results_\\\$i.OLD\\\$j".
"""
function determine_topology_dataset(path, save, autoclean, options::Options)
    if isdirpath(path)
        path = dirname(path)
    end
    if startswith(basename(path), "results_") && isfile(joinpath(path), "data")
        resultdir = path
        path = only(readlines(joinpath(path, "data")))
        alreadydone = String[]
        for _f in readdir(resultdir; join=true)
            basename(_f) == "data" && continue
            io = open(_f, "r+")
            pos = max(filesize(io) - 2, 0)
            seek(io, pos)
            while pos > 0 && read(io, Char) != '\n'
                pos -= 1
                seek(io, pos)
            end
            truncate(io, pos+(pos>0)) # remove the last line if incomplete
            close(io)
            for l in eachline(_f)
                _splits = split(l, '/')
                isempty(last(_splits)) && continue
                push!(alreadydone, join(_splits[1:end-2], '/'))
            end
        end
        files = recursive_readdir(path)
        setdiff!(files, alreadydone)
    else
        resultdir = tmpexportname(dirname(path), "", "results", "")
        mkdir(resultdir)
        open(joinpath(resultdir, "data"), "w") do f
            println(f, path)
        end
        files = recursive_readdir(path)
    end
    resultdir::String
    files::Vector{String}

    @threads for file in files
        f = joinpath(path, file)
        # threadid() == 1 && @show f # to find infinite loops: the last one printed is probably running

        genomes::Vector{Tuple{Vector{Int},TopologyResult}} = try
            topological_genome(UnderlyingNets(parse_chemfile(f, options)))
        catch e
            (options.throw_error || isinterrupt(e)) && rethrow()
            [(Int[], TopologyResult(string(e)))]
        end
        if isempty(genomes)
            push!(genomes, (Int[], TopologyResult("")))
        end
        for (j, (_, genome)) in enumerate(genomes)
            newname = length(genomes) == 1 ? file * '/' : file * '/' * string(j)
            open(joinpath(resultdir, string(threadid())), "a") do results
                io = IOContext(results, :compact => true)
                println(io, newname, '/', genome)
            end
        end
    end

    ret = Pair{String,TopologyResult}[]
    for _f in readdir(resultdir; join=true)
        basename(_f) == "data" && continue
        for l in eachline(_f)
            isempty(l) && continue
            splits = split(l, '/')
            _genome = pop!(splits)
            isempty(splits[end]) && pop!(splits)
            push!(ret, Pair(join(splits, '/'), parse(TopologyResult, _genome)))
        end
    end
    result::Dict{String,TopologyResult} = Dict(ret)
    if save
        i = 0
        tmpresultdir = resultdir*".OLD"*string(i)
        while ispath(tmpresultdir)
            i += 1
            tmpresultdir = resultdir*".OLD"*string(i)
        end
        mv(resultdir, tmpresultdir)
        Serialization.serialize(resultdir, result)
        println("Topologies of ", path, " saved at ", resultdir)
        if autoclean
            rm(tmpresultdir; recursive=true)
        else
            println("Temporary files kept at ", tmpresultdir)
        end
    elseif autoclean
        rm(resultdir; recursive=true)
    else
        println("Temporary files kept at ", resultdir)
    end
    return result
end
function determine_topology_dataset(path, save=true, autoclean=true; kwargs...)
    opts, restore_warns = db_options(; kwargs...)
    ret = determine_topology_dataset(path, save, autoclean, opts)
    restore_warns && (DOWARN[] = true)
    ret
end



"""
    guess_topology_dataset(path, save, autoclean, options::Options)
    guess_topology_dataset(path, save=true, autoclean=true; kwargs...)

Given a path to a directory containing structure input files, guess the topology of each
structure within the directory using [`guess_topology`](@ref).
Return a dictionary linking each file name to the result.
The result is the corresponding topology name, if known, or the topological
genome preceded by an "UNKNOWN" mention otherwise. In case of error, the result
is the exception preceded by a "FAILED with" mention. Finally, if the input does
not represent a periodic structure, the result is "non-periodic".

It is strongly recommended to toggle warnings off (through [`toggle_warning`](@ref)) and
not to export any file since those actions may critically reduce performance,
especially for numerous files.

The `save` and `autoclean` arguments work identically to their counterpart for
[`determine_topology_dataset`](@ref).
"""
function guess_topology_dataset(path, save, autoclean, options::Options)
    if isdirpath(path)
        path = dirname(path)
    end
    if startswith(basename(path), "guessed_") && isfile(joinpath(path), "data")
        resultdir = path
        path = only(readlines(joinpath(path, "data")))
        alreadydone = String[]
        for _f in readdir(resultdir; join=true)
            basename(_f) == "data" && continue
            io = open(_f, "r+")
            pos = max(filesize(io) - 2, 0)
            seek(io, pos)
            while pos > 0 && read(io, Char) != '\n'
                pos -= 1
                seek(io, pos)
            end
            truncate(io, pos+(pos>0)) # remove the last line if incomplete
            close(io)
            for l in eachline(_f)
                isempty(l) && continue
                _splits = split(l, '/')
                push!(alreadydone, join(_splits[1:end-1], '/'))
            end
        end
        files = recursive_readdir(path)
        setdiff!(files, alreadydone)
    else
        resultdir = tmpexportname(dirname(path), "", "guessed", "")
        mkdir(resultdir)
        open(joinpath(resultdir, "data"), "w") do f
            println(f, path)
        end
        files = recursive_readdir(path)
    end
    resultdir::String
    files::Vector{String}

    @threads for file in files
        f = joinpath(path, file)

        genome::TopologicalGenome = try
            guess_topology(f, options)
        catch e
            (options.throw_error || isinterrupt(e)) && rethrow()
            TopologicalGenome(escape_string(string(e)::String))
        end
        open(joinpath(resultdir, string(threadid())), "a") do results
            println(results, file, '/', genome)
        end
    end

    ret = Pair{String,TopologicalGenome}[]
    for _f in readdir(resultdir; join=true)
        basename(_f) == "data" && continue
        for l in eachline(_f)
            splits = split(l, '/')
            isempty(last(splits)) && @show _f, l
            _genome = pop!(splits)
            push!(ret, Pair(join(splits, '/'), parse(TopologicalGenome, _genome)))
        end
    end
    result::Dict{String,TopologicalGenome} = Dict(ret)
    if save
        i = 0
        tmpresultdir = resultdir*".OLD"*string(i)
        while ispath(tmpresultdir)
            i += 1
            tmpresultdir = resultdir*".OLD"*string(i)
        end
        mv(resultdir, tmpresultdir)
        Serialization.serialize(resultdir, result)
        println("Guessed topologies of ", path, " saved at ", resultdir)
        if autoclean
            rm(tmpresultdir; recursive=true)
        else
            println("Temporary files kept at ", tmpresultdir)
        end
    elseif autoclean
        rm(resultdir; recursive=true)
    else
        println("Temporary files kept at ", resultdir)
    end
    return result
end
function guess_topology_dataset(path, save=true, autoclean=true; kwargs...)
    opts, restore_warns = db_options(; structure=StructureType.Guess, kwargs...)
    ret = guess_topology_dataset(path, save, autoclean, opts)
    restore_warns && (DOWARN[] = true)
    ret
end
