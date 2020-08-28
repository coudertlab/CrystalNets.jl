## Common utility functions

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
    issingular(x::SMatrix{3,3,T,9}) where T<:Rational

Test whether a 3x3 matrix is singular. The input matrix must have a wide enough
type to avoid avorflows. If an overflow happens however, it is detected and results
in an error (rather than silently corrupting the result).
"""
function issingular(x::SMatrix{3,3,T,9}) where T<:Rational
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
    return y11 * y22 == y12 * y21
    # This can overflow so the input matrix should already have a wide enough type
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

    i = nextind(l, i)
    start = 0
    while i <= n
        if isspace(l[i])
            i = nextind(l, i)
        else
            if l[i] == '#'
                i = findnext(isequal('\n'), l, i)
                i == nothing && return (0, 0, 0)
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
            throw("Invalid syntax: opening multiline field at position $start is not closed")
        end
        if inquote
            throw("Invalid syntax: opening quote $quotesymb at position $start is not closed")
        end
    end
    return (start, n, n)
end
