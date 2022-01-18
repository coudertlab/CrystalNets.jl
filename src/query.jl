## CrystalNets main API

"""
    topological_genome(net::CrystalNet{D,T})::String where {D,T}

Compute the topological genome of a net. The topological genome is an invariant
if the net, meaning that it does not depend on its representation. It is also
the string representation of a D-periodic graph such that
`PeriodicGraph{D}(topological_genome(net))` is isomorphic to `net.graph` (except
possibly if the `ignore_types` option is unset).

!!! info
    Options must be passed directly within `net`.
"""
function topological_genome(net::CrystalNet{D,T})::String where {D,T}
    isempty(net.pos) && return "non-periodic"
    if net.options.ignore_types
        net = CrystalNet{D,T}(net.cell, fill(Symbol(""), length(net.types)), net.pos,
                              net.graph, net.options)
    end
    if !net.options.skip_minimize
        flag = true
        try
            net = minimize(net)
            flag = false
        catch e
            if T == BigInt || !(e isa OverflowError || e isa InexactError)
                rethrow()
            end
        end
        if flag # not in the catch to avoid a StackOverflow of errors in case something goes wrong
            newnet = CrystalNet{D,widen(soft_widen(T))}(net; ignore_types=false)
            return topological_genome(newnet)
        end
    end
    if isempty(net.options._pos) # could not be exported before
        export_default(net, "net", net.options.name, net.options.export_net; repeats=2)
    end

    try
        return topological_key(net)
    catch e
        if T == BigInt || !(e isa OverflowError || e isa InexactError)
            rethrow()
        end
    end
    return topological_genome(CrystalNet{D,widen(soft_widen(T))}(net; skip_minimize=true))
end
topological_genome(net::CrystalNet{0,T}) where {T} = "non-periodic"

"""
    topological_genome(g::PeriodicGraph, options::Options=Options())
    topological_genome(g::PeriodicGraph; kwargs...)

Compute the topological genome of a periodic graph.
"""
function topological_genome(g::PeriodicGraph, options::Options)
    net = CrystalNet(g, options)
    return topological_genome(net)
end
topological_genome(g::PeriodicGraph; kwargs...) = topological_genome(g, Options(; kwargs...))

function _loop_group!(ex, id, net, group)
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
                currallowed = first(group)[2].options.dimensions
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
    topological_genome(group::CrystalNetGroup)::Vector{Tuple{Vector{Int},String}}

Compute the topological genome of each subnet stored in `group`.

!!! info
    Options must be passed directly within the subnets.
"""
function topological_genome(group::CrystalNetGroup)
    ret = Tuple{Vector{Int},String}[]
    @loop_group for (id, net) in group
        push!(ret, (id, topological_genome(net)))
    end
    return ret
end

"""
    recognize_topology(genome::AbstractString, arc=CRYSTAL_NETS_ARCHIVE)
    recognize_topology(genomes::Vector{Tuple{Vector{Int},String}}, arc=CRYSTAL_NETS_ARCHIVE)

Attempt to recognize a topological genome from an archive of known genomes.

The second form is used on the output of topological_genome(::CrystalNetGroup).
"""
function recognize_topology(genome::AbstractString, arc=CRYSTAL_NETS_ARCHIVE)
    genome == "non-periodic" && return genome
    get(arc, genome, startswith(genome, "unstable") ? genome : "UNKNOWN "*genome)
end

function recognize_topology(genomes::Vector{Tuple{Vector{Int},String}}, arc=CRYSTAL_NETS_ARCHIVE)
    return [(src, recognize_topology(genome, arc)) for (src, genome) in genomes]
end

"""
    determine_topology(path, options::Options)
    determine_topology(path; kwargs...)

Compute the topology of the structure described in the file located at `path`.
This is essentially equivalent to calling
`recognize_topology(topological_genome(CrystalNetGroup(parse_chemfile(path, options))))`.
In the case where the structure is not made of interpenetrating nets, return the name of
the only net.
"""
function determine_topology(path, options::Options)
    genomes::Vector{Tuple{Vector{Int},String}} =
        topological_genome(CrystalNetGroup(parse_chemfile(path, options)))
    if length(genomes) == 1
        return recognize_topology(genomes[1][2])
    end
    length(genomes) == 0 && return "non-periodic"
    return recognize_topology(genomes)
end
determine_topology(path; kwargs...) = determine_topology(path, Options(; kwargs...))

"""
    determine_topologies(path, options::Options)
    determine_topologies(path; kwargs...)

Compute the topology of the files at `path`.
Return a triplet `(tops, genomes, failed)` where:
- `tops` is a dictionary linking each file name to its topology name, if recognised, or
  topological genome otherwise.
  If an input file "X" contains interpenetrating nets, they will automatically be separated
  and each subnet (called "X/1", "X/2", ...) will be linked to its corresponding topology.
- `genomes` contain the topological genomes that have been encountered but were not part of
  the archive.
- `failed` is a dictionary linking each file name to the tuple `(e, bt)` that stores the
  information on the cause of the failure of CrystalNets. It can be visualised with
  `Base.showerror(stdout, e, bt)`.
"""
function determine_topologies(path, options::Options)
    dircontent = isdir(path) ? collect(enumerate(readdir(path; join=true))) : [(1, path)]
    n = length(dircontent)
    ret = [Pair{String,String}[] for _ in 1:n]
    newgenomes = [String[] for _ in 1:n]
    failed::Vector{Union{Missing, Pair{String, Tuple{Exception, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}}}} =
        fill(missing, n)
    @threads for (i, f) in dircontent
        name = last(splitdir(f))
        # println(name)
        genomes::Vector{Tuple{Vector{Int},String}} = try
            topological_genome(CrystalNetGroup(parse_chemfile(f, options)))
        catch e
            if e isa InterruptException ||
              (e isa TaskFailedException && e.task.result isa InterruptException)
                rethrow()
            end
            failed[i] = (name => (e, catch_backtrace()))
            Tuple{Vector{Int},String}[]
        end
        if isempty(genomes)
            push!(genomes, (Int[], "non-periodic"))
        end
        for (j, (_, genome)) in enumerate(genomes)
            newname = length(genomes) == 1 ? name : name * '/' * string(j)
            id = recognize_topology(genome)
            if startswith(id, "UNKNOWN")
                push!(newgenomes[i], id[9:end])
            end
            push!(ret[i], (newname => id))
        end
    end
    newgens = sort!(collect(Iterators.flatten(newgenomes)))
    return Dict(Iterators.flatten(ret)), unique!(newgens),
           Dict(skipmissing(failed))
end
determine_topologies(path; kwargs...) = determine_topologies(path, db_options(; kwargs...))



macro ifvalidgenomereturn(opts, msg, skipcrystal=false)
    crystaldef = skipcrystal ? nothing : :(crystal = parse_chemfile(path, $opts))
    msg = "(found by $msg)"
    ifprintinfo = isempty(msg) ? nothing : :(CrystalNets.@ifwarn @info $msg)

    return esc(quote
    $crystaldef
    group = CrystalNetGroup(crystal)
    dim, subnets = isempty(group.D3) ? isempty(group.D2) ? (1, group.D1) : (2, group.D2) : (3, group.D3)
    for (_, net) in subnets
        sig = string(net.graph)::String
        sig = get(encountered_graphs, sig, sig)
        if haskey(encountered_genomes, sig)
            encountered_genomes[sig] += 1
        else
            genome = recognize_topology(topological_genome(net)::String)
            if dim == maxdim && !startswith(genome, "UNKNOWN") &&
                    !startswith(genome, "unstable") && !startswith(genome, "non-periodic")
                $ifprintinfo
                return genome
            end
            encountered_graphs[sig] = genome
            encountered_genomes[genome] = 1
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
function guess_topology(path, defopts)
    maxdim = maximum(defopts.dimensions; init=0)
    encountered_graphs = Dict{String,String}()
    encountered_genomes = Dict{String,Int}()
    if defopts.name == "unnamed"
        defopts = Options(defopts; name=splitext(splitdir(path)[2])[1])
    end

    @ifvalidgenomereturn defopts "using clustering_mode=ClusteringMode.Guess" false

    # dryrun = Dict{Symbol,Union{Nothing,Set{Symbol}}}()
    # crystal = parse_chemfile(path, Options(defopts; dryrun))
    atoms = Set{Symbol}(crystal.types)

    if any(x -> (at = atomic_numbers[x]; (ismetal[at] | ismetalloid[at])), atoms)
        @ifvalidgenomereturn Options(defopts; clustering_mode=ClusteringMode.MOF) "using MOF clusters"
    end

    defopts = Options(; clustering_mode=ClusteringMode.Auto)
    # if defopts.cutoff_coeff == 0.75 # if default cutoff was used, try larger
    #     @ifvalidgenomereturn Options(defopts; cutoff_coeff=0.85) "using longer cutoff"
    # end
    # invalidatoms = union(get(dryrun, :invalidatoms, Set{Symbol}()),
    #                      get(dryrun, :suspect_smallbonds, Set{Symbol}()))
    # hashydrogens = :H ∈ atoms
    # if hashydrogens
    #     @ifvalidgenomereturn Options(defopts; ignore_atoms=(:H,)) "ignoring H"
    #     # delete!(invalidatoms, :H)
    # end
    # for a in invalidatoms
    #     @ifvalidgenomereturn Options(defopts; ignore_atoms=tuple(a)) "ignoring $a"
    #     if hashydrogens
    #         @ifvalidgenomereturn Options(defopts; ignore_atoms=(a, :H)) "ignoring H and $a"
    #     end
    # end
    # if any(==(:C), atoms)# && :C ∉ invalidatoms # organic solvent
    #     @ifvalidgenomereturn Options(defopts; ignore_atoms=(:C,)) "ignoring C"
    # end
    # if haskey(dryrun, :try_Input_bonds)
    #     @ifvalidgenomereturn Options(defopts; bonding_mode=BondingMode.Input) "using input bonds"
    # end
    # @ifvalidgenomereturn Options(defopts; ignore_low_occupancy=true) "removing atoms with occupancy < 0.5"
    if :Al ∈ atoms && :P ∈ atoms # ALPO datastructure
        # flag = false
        # for i in vertices(crystal.graph)
        #     crystal.types[i] === :O || continue
        #     currt = :O
        #     for x in neighbors(crystal.graph, i)
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
    #     @ifvalidgenomereturn Options(defopts; bonding_mode=BondingMode.Guess) "discarding input bonds and guessing them"
    # end
    # if haskey(dryrun, :collisions)
    #     @ifvalidgenomereturn Options(defopts; authorize_pruning=false) "retaining all colliding atoms"
    # end

    # Finally, if everything fails, return the one encountered most
    maxencounter = 0
    most_plausible_genome = ""
    for (genome, i) in encountered_genomes
        if i > maxencounter
            maxencounter = i
            most_plausible_genome = genome
        end
    end
    return isempty(most_plausible_genome) ? begin
            issetequal(defopts.dimensions, [1,2,3]) ? "non-periodic" :
                "FAILED (perhaps because no net with suitable dimension in $(defopts.dimensions) was found)"
            end : most_plausible_genome
end
guess_topology(path; kwargs...) = guess_topology(path, Options(clustering_mode=ClusteringMode.Guess; kwargs...))


"""
    guess_topologies(path, options)
    guess_topologies(path; kwargs...)

Attempt to determine the topology of the files stored in the directory at `path`
by using `guess_topology` on each.

Similarly to `recognize_topologies`, return a triplet `(tops, genomes, failed)`
where `tops` is a dictionary linking the name of the files to the name of the
topology (or the topological genome if none), `genomes` contain the genomes
not part of the archive and `failed` is a dictionary linking the name of the
structure to the error and backtrace justifying its failure.
"""
function guess_topologies(path, options)
    dircontent = isdir(path) ? collect(enumerate(readdir(path))) : [(1, path)]
    ret = Vector{Pair{String,String}}(undef, length(dircontent))
    newgenomes::Vector{Union{Missing, String}} = fill(missing, length(dircontent))
    failed::Vector{Union{Missing, Pair{String, Tuple{Exception, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}}}} =
        fill(missing, length(dircontent))
    @threads for (i, f) in dircontent
        name = last(splitdir(f))
        # println(name)
        result::String = try
            guess_topology(joinpath(path, f), options)
        catch e
            if e isa InterruptException ||
              (e isa TaskFailedException && e.task.result isa InterruptException)
                rethrow()
            end
            failed[i] = (name => (e, catch_backtrace()))
            "FAILED"
        end
        ret[i] = (name => result)
        if startswith(result, "UNKNOWN")
            newgenomes[i] = result[9:end]
        end
    end
    newgens = sort!(collect(skipmissing(newgenomes)))
    return Dict(ret), unique!(newgens), Dict(skipmissing(failed))
end
guess_topologies(path; kwargs...) = guess_topologies(path, db_options(; kwargs...))


"""
    topologies_dataset(path, save, autoclean, options::Options)
    topologies_dataset(path, save=true, autoclean=true; kwargs...)

Given a path to a directory containing structure input files, compute the
topology of each structure within the directory.
Return a dictionary linking each file name to the result.
The result is the corresponding topology name, if known, or the topological
genome preceded by an "UNKNOWN" mention otherwise. In case of error, the result
is the exception preceded by a "FAILED with" mention. Finally, if the input does
not represent a periodic structure, the result is "non-periodic".

This function is similar to [`determine_topologies`](@ref)`, but targets larger datasets,
for which performance is critical. In particular, no attempt to recognise the
found topologies is performed: only the topological key is returned.

It is strongly recommended to toggle warnings off (through `toggle_warning`) and
not to export any file since those actions may critically reduce performance,
especially for numerous files.

If `save` is set, the result is also stored in a julia serialized file located at
"\$path/../results_\$i" where `i` is the lowest integer such that this path does
not already exist at the start of the computation.
While processing, this path will be used to create a directory storing the
current state of the computation: to continue an interrupted computation, simply
pass this temporary directory as the path. If `autoclean` is set, this directory
is removed at the end if the computation was successful.

If `save` is set and `autoclean` is unset, the directory of temporary files will
be renamed into "\$path/../results_\$i.OLD\$j".
"""
function topologies_dataset(path, save, autoclean, options::Options)
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
                isempty(l) && continue
                _i = 1; while l[_i] != '/'; _i+=1; end
                push!(alreadydone, joinpath(path, l[1:_i-1]))
            end
        end
        files = readdir(path; join=true, sort=true)
        setdiff!(files, alreadydone)
    else
        resultdir = tmpexportname(dirname(path), "", "results", "")
        mkdir(resultdir)
        open(joinpath(resultdir, "data"), "w") do f
            println(f, path)
        end
        files = readdir(path; join=true, sort=true)
    end

    @threads for f in files
        file = splitdir(f)[2]

        genomes::Vector{Tuple{Vector{Int},String}} = try
            topological_genome(CrystalNetGroup(parse_chemfile(f, options)))
        catch e
            if e isa InterruptException ||
              (e isa TaskFailedException && e.task.result isa InterruptException)
                rethrow()
            end
            [(Int[], "FAILED with "*escape_string(string(e)))]
        end
        if isempty(genomes)
            push!(genomes, (Int[], "non-periodic"))
        end
        for (j, (_, genome)) in enumerate(genomes)
            newname = length(genomes) == 1 ? file * '/' : file * '/' * string(j)
            open(joinpath(resultdir, string(threadid())), "a") do results
                println(results, newname, '/', genome)
            end
        end
    end

    ret = Pair{String,String}[]
    for _f in readdir(resultdir; join=true)
        basename(_f) == "data" && continue
        for l in eachline(_f)
            isempty(l) && continue
            _i = 1; while l[_i] != '/'; _i+=1; end
            _j = _i+1; while l[_j] != '/'; _j+=1; end
            push!(ret, Pair(l[1:_j-1-(_i==_j)], l[_j+1:end]))
        end
    end
    result::Dict{String,String} = Dict(ret)
    if save
        i = 0
        tmpresultdir = resultdir*".OLD"*string(i)
        while ispath(tmpresultdir)
            i += 1
            tmpresultdir = resultdir*".OLD"*string(i)
        end
        mv(resultdir, tmpresultdir)
        Serialization.serialize(resultdir, result)
        println("Topologies of $path saved at $resultdir")
        if autoclean
            rm(tmpresultdir; recursive=true)
        else
            println("Temporary files kept at $tmpresultdir")
        end
    elseif autoclean
        rm(resultdir; recursive=true)
    else
        println("Temporary files kept at $resultdir")
    end
    return result
end
function topologies_dataset(path, save=true, autoclean=true; kwargs...)
    topologies_dataset(path, save, autoclean, db_options(; kwargs...))
end
