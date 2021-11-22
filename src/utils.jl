## Common utility functions

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

function export_default(c, obj=string(typeof(c)), _name=nothing, path=tempdir())
    if !isempty(path)
        name = tmpexportname(path, obj*'_', _name, ".vtf")
        truepath = replace(joinpath(path, name), ('\\' => "/"))
        println("Export of $obj is enabled: saving file at $truepath")
        try
            export_vtf(truepath, c, 6)
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
        @assert e isa OverflowError
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
        @assert e isa OverflowError
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
        @assert e isa OverflowError
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
    @assert _n == 3
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
    @assert _n == 2
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
    @assert _n == 1
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
