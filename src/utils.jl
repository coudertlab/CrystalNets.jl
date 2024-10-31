## Common utility functions

@static if VERSION < v"1.8-"
    # copy-pasted from the implementation of @lazy_str in Base
    macro lazy_str(text)
        parts = Any[]
        lastidx = idx = 1
        while (idx = findnext('$', text, idx)) !== nothing
            lastidx < idx && push!(parts, text[lastidx:idx-1])
            idx += 1
            expr, idx = Meta.parseatom(text, idx; filename=string(__source__.file))
            push!(parts, esc(expr))
            lastidx = idx
        end
        lastidx <= lastindex(text) && push!(parts, text[lastidx:end])
        :(string($(parts...)))
    end
end

"""
    @toggleassert expression

Internal macro used to assert and expression conditionally on a build-time constant.
To toggle on or off these assertions, the constant has to be modified in the source code
and the module rebuilt afterwards.
"""
macro toggleassert(ex...)
    # Change the following boolean and reload CrystalNets to enable/disable assertions
    ENABLE_ASSERT = true
    if ENABLE_ASSERT
        return esc(quote
            @assert $(ex...)
        end)
    end
    return nothing
end

function no_metadata_metafmt(level::Logging.LogLevel, _module, group, id, file, line)
    @nospecialize
    color = (@static VERSION < v"1.12.0-" ? Logging : Base.CoreLogging).default_logcolor(level)
    prefix = string(level == Logging.Warn ? "Warning" : string(level), ':')
    return color, prefix, ""
end
const minimal_logger = Logging.ConsoleLogger(; meta_formatter=no_metadata_metafmt)


macro ifwarn(ex)
    return quote
        if (DOWARN[]::Bool)
            @static if VERSION < v"1.7-"
                $(esc(ex))
            else
                Base.with_logger(minimal_logger) do
                    $(esc(ex))
                end
            end
        end
    end
end

macro iferror(ex)
    return quote
        if (DOERROR[]::Bool)
            @static if VERSION < v"1.7-"
                $(esc(ex))
            else
                Base.with_logger(minimal_logger) do
                    $(esc(ex))
                end
            end
        end
    end
end


function recursive_readdir!(stored, prefix, path)
    for f in readdir(path; sort=false, join=false)
        name = joinpath(prefix, f)
        file = joinpath(path, f)
        if isdir(file)
            recursive_readdir!(stored, name, file)
        elseif !isempty(splitext(f)[2]) # otherwise, not a crystallographic file
            push!(stored, name)
        end
    end
    nothing
end
function recursive_readdir(path)
    ret = String[]
    while isdirpath(path)
        path = path[1:end-1]
    end
    sizehint!(ret, length(readdir(path; sort=false, join=false)))
    recursive_readdir!(ret, "", path)
    return ret
end


function tmpexportname(path, pre, name, ext)
    if name isa Nothing
        return tempname(path; cleanup=false)*ext
    end
    pre0 = pre*(name::String)
    fastpath = joinpath(path, pre0*ext)
    ispath(fastpath) || return fastpath
    paths = Set{String}(readdir(path; sort=false, join=false))
    i = 1
    pre1 = pre0*"_other"
    x = pre1*string(i)*ext
    while x ∈ paths
        i += 1
        x = pre1*string(i)*ext
    end
    return joinpath(path, x)
end

"""
    export_default(c::Union{PeriodicGraph,CrystalNet,Crystal}, obj=nothing, name=nothing, path=tempdir(); repeats=nothing)

Export a VTF representation of an object at the given `path`.

`obj` is a `String` describing the nature of the object, such as "net", "clusters" or
"subnet" for example. Default is `string(typeof(c))`.

`name` is a `String` inserted in the exported file name. Default is a [`tempname`](https://docs.julialang.org/en/v1/base/file/#Base.Filesystem.tempname).

`repeats` is the maximum distance between a represented atom out of the unit cell and one
inside. Default is between 2 and 6, depending on `obj` and the size of the graph.
"""
export_default(g::PeriodicGraph, args...; kwargs...) = export_default(CrystalNet(g), args...; kwargs...)
function export_default(c, obj=nothing, name=nothing, path=tempdir(); repeats=nothing)
    _repeats = repeats isa Integer ? repeats : begin
        if obj isa AbstractString && (obj == "net" || startswith(obj, "clusters") || startswith(obj, "subnet"))
            2
        else
            nv(c.pge.g) == 0 ? 1 : clamp(fld(600, nv(c.pge.g)), 2, 6)
        end
    end
    if !isempty(path)
        _name = tmpexportname(path, (obj isa Nothing ? string(typeof(c)) : obj)*'_', name, ".vtf")
        truepath = replace(joinpath(path, _name), ('\\' => "/"))
        if obj isa Nothing
            println(lazy"Saving file (representing $(typeof(c))) at $truepath")
        else
            println(lazy"Export of $obj is enabled: saving file at $truepath")
        end
        try
            export_vtf(truepath, c, _repeats)
        catch e
            if e isa SystemError
                @error lazy"Failed to export because of the following error: $e"
            else
                rethrow()
            end
        end
    end
    nothing
end

function db_options(; kwargs...)
    restore_warns = false
    if DOWARN[] && !get(kwargs, :force_warn, false)
        restore_warns = true
        @ifwarn begin
            @warn "Warnings may critically alter performance for this operation and were turned off."
            @info "Use CrystalNets.toggle_warning(false) or the --no-warn option to remove this warning"
            @info "Use the force_warn option to force printing warnings"
        end
        DOWARN[] = false
    end
    # if kwargs explicitly set one of the two, it will take precedence
    return Options(; export_input="", export_attributions="", export_clusters="",
                     export_net="", export_subnets="", kwargs...),
           restore_warns
end


const atomic_numbers = IdDict{Symbol, Int}(
    :H => 1,
    :D => 1,
    :He => 2,
    :Li => 3,
    :Be => 4,
    :B => 5,
    :C => 6,
    :N => 7,
    :O => 8,
    :F => 9,
    :Ne => 10,
    :Na => 11,
    :Mg => 12,
    :Al => 13,
    :Si => 14,
    :P => 15,
    :Pc => 15, # fictitious atom representing an organic P
    :S => 16,
    :Ss => 16, # fictitious atom representing an organic S
    :Cl => 17,
    :Ar => 18,
    :K => 19,
    :Ca => 20,
    :Sc => 21,
    :Ti => 22,
    :V => 23,
    :Cr => 24,
    :Mn => 25,
    :Fe => 26,
    :Co => 27,
    :Ni => 28,
    :Cu => 29,
    :Zn => 30,
    :Ga => 31,
    :Ge => 32,
    :As => 33,
    :Se => 34,
    :Br => 35,
    :Kr => 36,
    :Rb => 37,
    :Sr => 38,
    :Y => 39,
    :Zr => 40,
    :Nb => 41,
    :Mo => 42,
    :Tc => 43,
    :Ru => 44,
    :Rh => 45,
    :Pd => 46,
    :Ag => 47,
    :Cd => 48,
    :In => 49,
    :Sn => 50,
    :Sb => 51,
    :Te => 52,
    :I => 53,
    :Xe => 54,
    :Cs => 55,
    :Ba => 56,
    :La => 57,
    :Ce => 58,
    :Pr => 59,
    :Nd => 60,
    :Pm => 61,
    :Sm => 62,
    :Eu => 63,
    :Gd => 64,
    :Tb => 65,
    :Dy => 66,
    :Ho => 67,
    :Er => 68,
    :Tm => 69,
    :Yb => 70,
    :Lu => 71,
    :Hf => 72,
    :Ta => 73,
    :W => 74,
    :Re => 75,
    :Os => 76,
    :Ir => 77,
    :Pt => 78,
    :Au => 79,
    :Hg => 80,
    :Tl => 81,
    :Pb => 82,
    :Bi => 83,
    :Po => 84,
    :At => 85,
    :Rn => 86,
    :Fr => 87,
    :Ra => 88,
    :Ac => 89,
    :Th => 90,
    :Pa => 91,
    :U => 92,
    :Np => 93,
    :Pu => 94,
    :Am => 95,
    :Cm => 96,
    :Bk => 97,
    :Cf => 98,
    :Es => 99,
    :Fm => 100,
    :Md => 101,
    :No => 102,
    :Lr => 103,
    :Rf => 104,
    :Db => 105,
    :Sg => 106,
    :Bh => 107,
    :Hs => 108,
    :Mt => 109,
    :Ds => 110,
    :Rg => 111,
    :Cn => 112,
    :Nh => 113,
    :Fl => 114,
    :Mc => 115,
    :Lv => 116,
    :Ts => 117,
    :Og => 118,
    :Uue => 119
)

const element_categories = Symbol[ # populated using PeriodicTable.jl
    :nonmetal, :noble, :metal, :metal, :metalloid, :nonmetal, :nonmetal,
    :nonmetal, :halogen, :noble, :metal, :metal, :metal, :metalloid, :nonmetal,
    :nonmetal, :halogen, :noble, :metal, :metal, :metal, :metal, :metal, :metal,
    :metal, :metal, :metal, :metal, :metal, :metal, :metal, :metalloid,
    :metalloid, :nonmetal, :halogen, :noble, :metal, :metal, :metal, :metal,
    :metal, :metal, :metal, :metal, :metal, :metal, :metal, :metal, :metal,
    :metal, :metalloid, :metalloid, :halogen, :noble, :metal, :metal,
    :lanthanide, :lanthanide, :lanthanide, :lanthanide, :lanthanide,
    :lanthanide, :lanthanide, :lanthanide, :lanthanide, :lanthanide,
    :lanthanide, :lanthanide, :lanthanide, :lanthanide, :lanthanide, :metal,
    :metal, :metal, :metal, :metal, :metal, :metal, :metal, :metal, :metal,
    :metal, :metal, :metal, :halogen, :noble, :metal, :metal, :actinide,
    :actinide, :actinide, :actinide, :actinide, :actinide, :actinide, :actinide,
    :actinide, :actinide, :actinide, :actinide, :actinide, :actinide, :actinide,
    :metal, :metal, :metal, :metal, :metal, :metal, :metal, :metal, :metal,
    :metal, :metal, :metal, :metal, :halogen, :noble, :metal]


function string_atomtype(t::Symbol)
    @static if VERSION < v"1.7-"
        replace(replace(string(t), "Pc" => 'P'), "Ss" => 'S')
    else
        replace(string(t), "Pc" => 'P', "Ss" => 'S')
    end
end

function representative_atom(t::Symbol, default::Int=0)
    at = get(atomic_numbers, t, 0)
    at == 0 || return t, at
    t === :* && return t, default
    styp = string(t)
    # @toggleassert length(styp) ≥ 2 not true for fake atoms like G
    t = islowercase(styp[2]) ? Symbol(styp[1:2]) : Symbol(styp[1])
    return t, get(atomic_numbers, t, default)
end

"""
    issingular(x::SMatrix{N,N,T}) where {N,T<:Rational}

Test whether a NxN matrix is singular.
"""
function issingular(x::SMatrix{N,N,T}) where {N,T<:Rational}
    try
        return iszero(det(x))
    catch e
        e isa OverflowError || rethrow()
    end
    return iszero(det(SMatrix{N,N,widen(T)}(x)))
end

function issingular(x::SMatrix{3,3,T,9})::Bool where T
@inbounds begin
    (i, j, k) = iszero(x[1,1]) ? (iszero(x[1,2]) ? (3,1,2)  : (2,1,3)) : (1,2,3)
    iszero(x[1,i]) && return true
    x1i = x[1,i]
    factj = (x[1,j] // x1i)
    factk = (x[1,k] // x1i)
    y11 = x[2,j] - factj * x[2,i]
    y12 = x[3,j] - factj * x[3,i]
    y21 = x[2,k] - factk * x[2,i]
    y22 = x[3,k] - factk * x[3,i]
    return widemul(y11, y22) == widemul(y12, y21)
    # This can overflow so the input matrix should already have a wide enough type
end
end

function issingular(x::SMatrix{2,2,T,4})::Bool where T
@inbounds begin
    if iszero(x[1,1])
        (iszero(x[1,2]) || iszero(x[2,1])) && return true
        return false
    end
    x12 = x[1,2] // x[1,1]
    return x[2,2] == widemul(x[2,1], x12)
    # This can overflow so the input matrix should already have a wide enough type
end
end

function issingular(x::SMatrix{1,1})::Bool
    iszero(@inbounds x[1,1])
end

function isrank3(x::AbstractMatrix{T}) where T<:Rational
    _n, m = size(x)
    @toggleassert _n == 3
    m < 3 && return false
    cols = collect(eachcol(x))
    u1 = pop!(cols)
    while iszero(u1) && !isempty(cols)
        u1 = pop!(cols)
    end
    n = length(cols)
    n < 2 && return false
    k = iszero(@inbounds u1[1]) ? iszero(@inbounds u1[2]) ? 3 : 2 : 1
    u1k = @inbounds u1[k]
    u2 = u1
    i = 1
    @inbounds while i < n
        u2 = cols[i]
        (u2[k]//u1k)*u1 == u2 || break
        i+=1
    end
    i == n && return false
    @inbounds for j in i+1:n
        a = SMatrix{3,3,T,9}(hcat(u1, u2, cols[j]))
        issingular(a) || return true
    end
    return false
end

function isrank2(x::AbstractMatrix{T}) where T<:Rational
    _n, n = size(x)
    @toggleassert _n == 2
    n < 2 && return false
    cols = collect(eachcol(x))
    u1 = pop!(cols)
    while iszero(u1) && !isempty(cols)
        u1 = pop!(cols)
    end
    length(cols) == 0 && return false
    # x is at least rank 1
    k = ifelse(iszero(@inbounds u1[1]), 2, 1)
    u1k = @inbounds u1[k]
    @inbounds for u2 in cols
        (u2[k]//u1k)*u1 == u2 || return true # found a non-colinear vector, thus rank 2
    end
    return false # only colinear vectors: rank 1
end

function isrank1(x::AbstractMatrix{T}) where T<:Rational
    @toggleassert size(x, 1) == 1
    @inbounds for y in eachcol(x)
        iszero(y[1]) || return true
    end
    return false
end


"""
    back_to_unit(r::Rational)

Return `x - floor(Int, x)`
"""
function back_to_unit(r::Rational)
    return Base.unsafe_rational(mod(numerator(r), denominator(r)), denominator(r))
end


"""
    nextword(l::String, k::Int)

Return the triplet of indices `(i, j, x)` such that `l[i:j]` is the next word in the
string `l` after position `k`.
Use `k = x` to get the following word, and so forth.

`(0, 0, 0)` is returned if there is no next word.
"""
function nextword(l::String, i::Int)
    n = lastindex(l)
    i == n && return (0, 0, 0)

    i::Int = nextind(l, i)
    start::Int = 0
    while i <= n
        if isspace(l[i])
            i = nextind(l, i)
        else
            if l[i] == '#'
                i = findnext(isequal('\n'), l, i)
                i === nothing && return (0, 0, 0)
            else
                start = i
                break
            end
        end
    end
    start == 0 && return (0, 0, 0)

    inquote = l[start] == '\'' || l[start] == '"'
    quotesymb = l[start]
    inmultiline = l[start] == ';' && l[prevind(l,start)] == '\n'
    instatus = inmultiline | inquote
    instatus && (start = nextind(l, start); i = start)
    while i <= n
        c = l[i]
        if !isspace(c)
            if instatus
                if inmultiline && c == ';' && l[prevind(l,i)] == '\n'
                    return (start, prevind(l, i-1), i)
                elseif inquote && c == quotesymb && i != n
                    nextl = l[nextind(l, i)]
                    if isspace(nextl) || nextl == '#'
                        return (start, prevind(l, i), i)
                    end
                end
            elseif c == '#' && isspace(l[prevind(l,i)])
                return (start, prevind(l, i), prevind(l, i))
            end
            i = nextind(l, i)
        elseif !instatus
            return (start, prevind(l, i), i)
        else
            i = nextind(l, i)
        end
    end
    if instatus
        if inmultiline
            error(lazy"Invalid syntax: opening multiline field at position $start is not closed")
        end
        if inquote
            error(lazy"Invalid syntax: opening quote $quotesymb at position $start is not closed")
        end
    end
    return (start, n, n)
end


function angle(p1, p2)
    return acosd(clamp(dot(p1, p2)/(norm(p1)*norm(p2)), -1.0, 1.0))
end
function dihedral(p1, p2, p3)
    β = angle(cross(p1, p2), cross(p2, p3))
    return isnan(β) ? 0.0 : β
end


function isinterrupt(@nospecialize(e))
    return e isa InterruptException || (e isa TaskFailedException && e.task.result isa InterruptException)
end
function isoverfloworinexact(@nospecialize(e))
    (e isa OverflowError || e isa InexactError) && return true
    if e isa TaskFailedException
        result = e.task.result
        return result isa OverflowError || result isa InexactError
    end
    if e isa CompositeException
        all(isoverfloworinexact, e) && return true
        return false
    end
    return false
end


# PlainChangesIterator for iterating over relevant permutation of vertices in unstable nets

"""
    PlainChangesIterator

Iterator over indices `i` of "plain changes", i.e. such that, starting from a vector of
`n` distinct elements and permutating successively each pair `(i, i+1)` with `i` given by
iterating over `PlainChangesIterator(n)` makes the vector undergo all the ``n!``
permutations.

See also: Knuth's P" algorithm for the method of plain changes, taken from
https://mathoverflow.net/a/279874, and the Steinhaus–Johnson–Trotter algorithm wikipedia
page https://en.wikipedia.org/wiki/Steinhaus%E2%80%93Johnson%E2%80%93Trotter_algorithm
"""
struct PlainChangesIterator
    n::Int
end
Base.length(pci::PlainChangesIterator) = factorial(pci.n) - 1
Base.eltype(::Type{PlainChangesIterator}) = Int

_init_state(pci) = (zeros(Int, pci.n), ones(Int8, pci.n))
function Base.iterate(pci::PlainChangesIterator, state=nothing)
    n = pci.n
    c, d = state===nothing ? _init_state(pci) : state
    j = n
    s = 0
    @label P4
    cj = c[j]
    q = cj + d[j]
    q < 0 && @goto P7
    q == j && @goto P6
    c[j] = q
    a = j - cj + s; b = j - q + s
    return min(a, b), (c, d)
    @goto P4
    @label P6
    j == 1 && return nothing
    s += 1
    @label P7
    d[j] = -d[j]
    j -= 1
    @goto P4
end

struct ContiguousPlainChangesIterator
    ranges::Vector{UnitRange{Int}}
    pcis::Vector{PlainChangesIterator}
end
function ContiguousPlainChangesIterator(ranges::Vector{UnitRange{Int}})
    pci = [PlainChangesIterator(length(r)) for r in ranges]
    ContiguousPlainChangesIterator(ranges, pci)
end

Base.IteratorSize(::Type{ContiguousPlainChangesIterator}) = Base.SizeUnknown()
Base.eltype(::Type{ContiguousPlainChangesIterator}) = Int

function Base.iterate(cpci::ContiguousPlainChangesIterator, _state=nothing)
    T = Tuple{Vector{Int}, Vector{Int8}}
    states = _state===nothing ? Union{Nothing,T}[_init_state(pci) for pci in cpci.pcis] : _state
    for (i, state) in enumerate(states)
        state===nothing && continue
        res = iterate(cpci.pcis[i], states[i])
        if res===nothing
            states[i] = nothing
            continue
        end
        a, newstate = res
        for j in 1:(i-1)
            states[j] = _init_state(cpci.pcis[j])
        end
        states[i] = newstate
        return first(cpci.ranges[i]) + a - 1, states
    end
    return nothing
end
