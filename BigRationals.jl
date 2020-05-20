module BigRationals

export BigRational

import Base: numerator, denominator, string, show, ==, hastypemax, iszero, Rational
import Base.GMP: Limb, ClongMax, CulongMax
import Base.GMP.MPZ: mpz_t


mutable struct BigRational <: Real
    num_alloc::Cint
    num_size::Cint
    num_d::Ptr{Limb}
    den_alloc::Cint
    den_size::Cint
    den_d::Ptr{Limb}

    function BigRational()
        b = init!(new())
        finalizer(cglobal((:__gmpq_clear, :libgmp)), b)
        return b
    end
end


# module MPQ

const mpq_t = Ref{BigRational}

gmpq(op::Symbol) = (Symbol(:__gmpq_, op), :libgmp)

for (op) in (:canonicalize, :init)
    op! = Symbol(op, :!)
    @eval $op!(x::BigRational) = (ccall($(gmpq(op)), Cvoid, (mpq_t,), x); x)
end

for (op, T) in (:set_ui => Culong, :set_si => Clong)
    op! = Symbol(op, :!)
    @eval begin
        $op!(x::BigRational, a, b) = (ccall($(gmpq(op)), Cvoid, (mpq_t, $T, Culong), x, a, b); x)
        $op(a, b) = $op!(BigRational(), a, b)
    end
end

for (op, T) in (:set => mpq_t, :set_z => mpz_t, :set_d => Cdouble, :set_num => mpz_t, :set_den => mpz_t)
    op! = Symbol(op, :!)
    @eval $op!(x::BigRational, a) = (ccall($(gmpq(op)), Cvoid, (mpq_t, $T), x, a); x)
end

for op in (:get_num, :get_den)
    op! = Symbol(op, :!)
    @eval begin
        $op!(x::BigInt, a::BigRational) = (ccall($(gmpq(op)), Cvoid, (mpz_t, mpq_t), x, a); x)::BigInt
        $op(a::BigRational) = $op!(BigInt(), a)
    end
end

for (op) in (:inv, :neg, :abs)
    op! = Symbol(op, :!)
    @eval begin
        $op!(x::BigRational, a::BigRational) = (ccall($(gmpq(op)), Cvoid, (mpq_t, mpq_t), x, a); x)
        $op!(x::BigRational) = $op!(x, x)
    end
end

for op in (:set, :set_z, :set_d, :set_num, :inv, :neg, :abs)
    op! = Symbol(op, :!)
    @eval $op(a) = $op!(BigRational(), a)
end

for (op) in (:add, :sub, :mul, :div)
    op! = Symbol(op, :!)
    @eval begin
        $op!(x::BigRational, a::BigRational, b::BigRational) = (ccall($(gmpq(op)), Cvoid, (mpq_t, mpq_t, mpq_t), x, a, b); x)
        $op(a::BigRational, b::BigRational) = $op!(BigRational(), a, b)
        $op!(x::BigRational, b::BigRational) = $op!(x, x, b)
    end
end

for op in (:mul_2exp, :div_2exp)
    op! = Symbol(op, :!)
    @eval begin
        $op!(x::BigRational, a::BigRational, b) = (ccall($(gmpq(op)), Cvoid, (mpq_t, mpq_t, Culong), x, a, b); x)
        $op(a::BigRational, b) = $op!(BigRational(), a, b)
        $op!(x::BigRational, b) = $op!(x, x, b)
    end
end

cmp(a::BigRational, b::BigRational) = Int(ccall((:__gmpq_cmp, :libgmp), Cint, (mpq_t, mpq_t), a, b))
equal(a::BigRational, b::BigRational) = Bool(ccall((:__gmpq_equal, :libgmp), Cint, (mpq_t, mpq_t), a, b))

cmp_z(a::BigRational, b::BigInt) = Int(ccall((:__gmpq_cmp_z, :libgmp), Cint, (mpq_t, mpz_t), a, b))
cmp_ui(a::BigRational, b, c) = Int(ccall((:__gmpq_cmp_ui, :libgmp), Cint, (mpq_t, Culong, Culong), a, b, c))
cmp_si(a::BigRational, b, c) = Int(ccall((:__gmpq_cmp_si, :libgmp), Cint, (mpq_t, Clong, Culong), a, b, c))

get_d(a::BigRational) = ccall((:__gmpq_get_d, :libgmp), Cdouble, (mpq_t,), a)
# TODO add set_f!(a::BigRational, a::MPF)
get_str!(x, a, b::BigRational) = (ccall((:__gmpq_get_str, :libgmp), Ptr{Cchar}, (Ptr{Cchar}, Cint, mpq_t), x, a, b); x)
set_str!(x::BigRational, a, b) = Int(ccall((:__gmpq_set_str, :libgmp), Cint, (mpq_t, Ptr{UInt8}, Cint), x, a, b))

# end # module MPQ

hastypemax(::Type{BigRational}) = true

BigRational(x::Union{Clong,Int32}) = set_si(x, one(Culong))
BigRational(x::ClongMax) = BigRational(Clong(x))
BigRational(x::BigInt) = set_num(x)
BigRational(x::Integer) = BigRational(BigInt(x))

function BigRational(x::Union{Clong,Int32}, y::Union{Culong,UInt32})
    num, den = divgcd(x, y)
    set_si(num, den)
end
function BigRational(x::ClongMax, y::ClongMax)
    num = Clong(flipsign(x, y))
    den = Culong(unsigned(abs(y)))
    BigRational(num, den)
end
BigRational(x::ClongMax, y::CulongMax) = BigRational(Clong(x), Culong(y))

function BigRational(x::Union{Culong,UInt32}, y::Union{Culong,UInt32})
    num, den = divgcd(x, y)
    set_ui(num, den)
end

function BigRational(x::Base.BitUnsignedSmall, y::ClongMax)
    y < 0 && return BigRational(-Clong(signed(widen(x))), Culong(unsigned(-y)))
    BigRational(Culong(x), unsigned(y))
end
function BigRational(x::UInt, y::ClongMax)
    y < 0 && return BigRational(-BigInt(x), -BigInt(y))
    BigRational(x, unsigned(y))
end
BigRational(x::CulongMax, y::CulongMax) = BigRational(Culong(x), Culong(y))

function BigRational(x::Integer, y::Integer)
    iszero(y) && return BigRational(Rational{BigInt}(flipsign(one(Clong), x), zero(Clong)))
    canonicalize!(set_den!(BigRational(x), y))
end

numerator(x::BigRational) = get_num(x)
denominator(x::BigRational) = get_den(x)
BigRational(x::Rational) = set_den!(BigRational(numerator(x)), denominator(x))
Rational{BigInt}(x::BigRational) = Base.unsafe_rational(BigInt, numerator(x), denominator(x))
Rational(x::BigRational) = Rational{BigInt}(x)

function string(x::BigRational; base::Integer = 10, pad::Integer = 1)
    string("BigRational(" * string(numerator(x); base, pad)   * ',' *
                            string(denominator(x); base, pad) * ')')
end
show(io::IO, x::BigRational) = print(io, string(x))
==(x::BigRational, y::BigRational) = equal(x, y)
iszero(x::BigRational) = x.num_size == 0

end
