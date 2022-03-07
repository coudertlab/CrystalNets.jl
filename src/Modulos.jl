# Inspired by the Mods.jl package

module Modulos

import Base: ==, +, -, *, inv, /, ^, hash, show, unsigned
import Base.Checked: mul_with_overflow

export Modulo, is_invertible

"""
`Modulo{p,T}(n)` creates a modular number in mod `p` with value `n%p` represented
by the integer type `T`.
Signed `T` will yield better performance.
"""
struct Modulo{p,T<:Integer} <: Real
    value::T
    function Modulo{p,T}(n::Integer) where {p,T}
        @assert p > 0 # This check is elided at runtime since p is a compile-time constant
        return new{p,T}(mod(n, p))
    end
end

(::Type{T})(x::Modulo) where {T<:Integer} = T(x.value)

function hash(x::Modulo{p}, h::UInt64=UInt64(0)) where p
    hash(Integer(x), hash(p, h))
end

==(x::Modulo{p}, y::Modulo{p}) where {p} = Integer(x) == Integer(y)

function +(x::Modulo{p,T1}, y::Modulo{p,T2}) where {p,T1,T2}
    T = promote_type(T1,T2)
    if T1 != T2
        tmin = typemin(T)
        tmax = typemax(T)
        if !(tmin <= typemin(T1) && tmin <= typemin(T2) &&
             tmax >= typemax(T1) && tmax >= typemax(T2))
            T = widen(T)
        end
    end
    # Next are some fast-paths if p is small enough.
    # Note that the checks are elided at runtime since p and T are compile-time constants
    if p < typemax(T) ÷ 2
        Modulo{p,T}((Integer(x) % T) + (Integer(y) % T))
    elseif p < typemax(unsigned(T)) ÷ 2
        U = unsigned(T)
        Modulo{p,T}((Integer(x) % U) + (Integer(y) % U))
    else
        V = widen(T)
        Modulo{p,T}((Integer(x) % V) + (Integer(y) % V))
    end
end

function -(x::Modulo{p,T}) where {p,T<:Signed}
    Modulo{p,T}(-Integer(x))
end
function -(x::Modulo{p,T}) where {p,T<:Unsigned}
    U = signed(T)
    y = Integer(x)
    y > (typemax(U) % T) ? Modulo{p,T}(-(y % widen(U))) : Modulo{p,T}(-signed(y))
end

-(x::Modulo, y::Modulo) = x + (-y)


function *(x::Modulo{p,T1}, y::Modulo{p,T2}) where {p,T1,T2}
    T = promote_type(T1,T2)
    if T1 != T2
        tmin = typemin(T)
        tmax = typemax(T)
        if !(tmin <= typemin(T1) && tmin <= typemin(T2) &&
             tmax >= typemax(T1) && tmax >= typemax(T2))
            T = widen(T)
        end
    end
    r, flag = mul_with_overflow((Integer(x) % T), (Integer(y) % T))
    flag ? Modulo{p,T}(widemul(Integer(x), Integer(y))) : Modulo{p,T}(r)
end


# """
# `is_invertible(x::Modulo)` determines if `x` is invertible.
# """
# is_invertible(x::Modulo{p}) where {p} = return p>1 && isone(gcd(Integer(x), p))

@noinline __throw_notinvertible(x) = error(x, " is not invertible")
inv(x::Modulo{1}) = __throw_notinvertible(x)
function inv(x::Modulo{p,T}) where {p,T}
    g, v, _ = gcdx(Integer(x), p)
    g == 1 || __throw_notinvertible(x)
    return Modulo{p,T}(v)
end

function /(x::Modulo{p}, y::Modulo{p}) where {p}
    return x * inv(y)
end

function ^(x::Modulo{p,T}, k::Integer) where {p,T}
    if k>0
        return Modulo{p,T}(powermod(Integer(x), k, p))
    end
    if k==0
        return one(Modulo{p})
    end
    Modulo{p,T}(powermod(Integer(inv(x)), -k, p))
end


+(x::Modulo{p,T}, k::Integer) where {p,T} = x + Modulo{p,T}(k)
+(k::Integer, x::Modulo) = x + k

-(x::Modulo, k::Integer) = x + (-k)
-(k::Integer, x::Modulo) = (-x) + k

*(x::Modulo{p,T}, k::Integer) where {p,T} = x * Modulo{p,T}(k)
*(k::Integer, x::Modulo) = x * k

/(x::Modulo{p,T}, k::Integer) where {p,T} = x / Modulo{p,T}(k)
/(k::Integer, x::Modulo{p,T}) where {p,T} = Modulo{p,T}(k) / x


function ==(x::Modulo{p,T1}, k::T2) where {p,T1,T2<:Integer}
    T = promote_type(T1,T2)
    T(x) == T(mod(k, p))
end
==(k::Integer, x::Modulo) = x == k

function show(io::IO, x::Modulo{p,T}) where {p,T}
    verbose = get(io, :typeinfo, Any) != Modulo{p,T}
    if verbose
        print(io, Modulo{p,T}, '(')
    end
    print(io, x.value)
    if verbose
        print(io, ')')
    end
end

end # end of module Modulos
