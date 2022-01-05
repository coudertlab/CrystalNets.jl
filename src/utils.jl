## Common utility functions

macro toggleassert(ex...)
    # Change the following boolean and reload CrystalNets to enable assertions
    ENABLE_ASSERT = false
    if ENABLE_ASSERT
        return esc(quote
            @assert $(ex...)
        end)
    end
    return nothing
end

function no_metadata_metafmt(level::Logging.LogLevel, _module, group, id, file, line)
    @nospecialize
    color = Logging.default_logcolor(level)
    prefix = string(level == Logging.Warn ? "Warning" : string(level), ':')
    return color, prefix, ""
end
const minimal_logger = Logging.ConsoleLogger(; meta_formatter=no_metadata_metafmt)


macro ifwarn(ex)
    return quote
        if (DOWARN[]::Bool)
            with_logger(minimal_logger) do
                $(esc(ex))
            end
        end
    end
end

function tmpexportname(path, pre, name, ext)
    if name isa Nothing
        return tempname(path; cleanup=false)*ext
    end
    i = 0
    pre = pre*name*'_'
    x = pre*string(i)*ext
    paths = Set{String}(readdir(path; sort=false, join=false))
    while x ∈ paths
        i += 1
        x = pre*string(i)*ext
    end
    return joinpath(path, x)
end

export_default(g::PeriodicGraph, args...; kwargs...) = export_default(CrystalNet(g), args...; kwargs...)
function export_default(c, obj=nothing, _name=nothing, path=tempdir(); repeats=6)
    if !isempty(path)
        name = tmpexportname(path, (obj isa Nothing ? string(typeof(c)) : obj)*'_', _name, ".vtf")
        truepath = replace(joinpath(path, name), ('\\' => "/"))
        if obj isa Nothing
            println("Saving file (representing $(typeof(c))) at $truepath")
        else
            println("Export of $obj is enabled: saving file at $truepath")
        end
        try
            export_vtf(truepath, c, repeats)
        catch e
            if e isa SystemError
                @error "Failed to export because of the following error: $e"
            else
                rethrow()
            end
        end
    end
    nothing
end

function db_options(; kwargs...)
    if DOWARN[]
        @error "Warnings may critically alter performance."
        @info "Use CrystalNets.toggle_warning(false) or the --no-warn option to remove them"
    end
    # if kwargs explicitly set one of the two, it will take precedence
    return Options(; export_input="", export_attributions="", export_clusters="", 
                     export_net="", kwargs...)
end


const atomic_numbers = IdDict{Symbol, Int}(
    :H => 1,
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
    :S => 16,
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


"""
    soft_widen(::Type)

Internal function used to selectively widen small integer and rational types.
This is useful to avoid overflow without sacrificing too much efficiency by
always having to resolve to very large types.
"""
soft_widen(::Type{T}) where {T} = T
soft_widen(::Type{Int32}) = Int64
soft_widen(::Type{Int16}) = Int32
soft_widen(::Type{Int8}) = Int16
soft_widen(::Type{Rational{T}}) where {T} = Rational{soft_widen(T)}

"""
    issingular(x::SMatrix{N,N,T}) where {N,T<:Rational}

Test whether a NxN matrix is singular. The input matrix must have a wide enough
type to avoid overflows. If an overflow happens however, it is detected and results
in an error (rather than silently corrupting the result).
"""
function issingular(x::SMatrix{N,N,T}) where {N,T<:Rational}
    try
        return iszero(det(x))
    catch e
        @toggleassert e isa OverflowError
        return iszero(det(SMatrix{N,N,widen(T)}(x)))
    end
end

function issingular(x::SMatrix{3,3,T})::Bool where T<:Rational
@inbounds begin
    (i, j, k) = iszero(x[1,1]) ? (iszero(x[1,2]) ? (3,1,2)  : (2,1,3)) : (1,2,3)
    iszero(x[1,i]) && return true
    U = widen(T)
    x1i = U(x[1,i])
    factj = (x[1,j] // x1i)
    factk = (x[1,k] // x1i)
    y11 = x[2,j] - factj * x[2,i]
    y12 = x[3,j] - factj * x[3,i]
    y21 = x[2,k] - factk * x[2,i]
    y22 = x[3,k] - factk * x[3,i]
    try
        return y11 * y22 == y12 * y21
    catch e
        @toggleassert e isa OverflowError
        return widemul(y11, y22) == widemul(y12, y21)
    end
    # This can overflow so the input matrix should already have a wide enough type
end
end

function issingular(x::SMatrix{2,2,T})::Bool where T<:Rational
@inbounds begin
    if iszero(x[1,1])
        (iszero(x[1,2]) || iszero(x[2,1])) && return true
        return false
    end
    U = widen(T)
    x12 = x[1,2] // U(x[1,1])
    return x[2,2] == try
        x[2,1] * x12
    catch e
        @toggleassert e isa OverflowError
        widemul(x[2,1], x12)
    end
    # This can overflow so the input matrix should already have a wide enough type
end
end

function issingular(x::SMatrix{1,1,T})::Bool where T<:Rational
    iszero(@inbounds x[1,1])
end

function isrank3(x::AbstractMatrix{T}) where T<:Rational
    _n, n = size(x)
    @toggleassert _n == 3
    n < 3 && return false
    cols = collect(eachcol(x))
    u1 = popfirst!(cols)
    u2 = u1
    while iszero(u1) && !isempty(cols)
        u1 = popfirst!(cols)
    end
    n = length(cols)
    k = iszero(@inbounds u1[1]) ? iszero(@inbounds u1[2]) ? 3 : 2 : 1
    i = 1
    @inbounds while i < n
        u2 = cols[i]
        widen(u2[k]//u1[k])*u1 == u2 || break
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
    u1 = popfirst!(cols)
    u2 = u1
    while iszero(u1) && !isempty(cols)
        u1 = popfirst!(cols)
    end
    n = length(cols)
    k = iszero(@inbounds u1[1]) ? 2 : 1
    i = 1
    @inbounds while i < n
        u2 = cols[i]
        widen(u2[k]//u1[k])*u1 == u2 || break
        i+=1
    end
    i == n && return false
    return true
end

function isrank1(x::AbstractMatrix{T}) where T<:Rational
    _n, n = size(x)
    @toggleassert _n == 1
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
    nextword(l, k)

Return the triplet of indices `(i, j, x)` such that `l[i:j]` is the next word in the
string `l` after position `k`.
Use `k = x` to get the following word, and so forth.

`(0, 0, 0)` is returned if there is no next word.
"""
function nextword(l, i)
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
                elseif inquote && c == quotesymb && i != n &&
                        (isspace(l[nextind(l, i)]) || l[nextind(l, i)] == '#')
                    return (start, prevind(l, i), i)
                end
            elseif c == '#'
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
            error("Invalid syntax: opening multiline field at position $start is not closed")
        end
        if inquote
            error("Invalid syntax: opening quote $quotesymb at position $start is not closed")
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

