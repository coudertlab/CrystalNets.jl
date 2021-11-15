

"""
    find_ratbasis(list::AbstractVector{<:StaticVector{N,<:Rational{T}}}) where {N,T}

Given a list of `N`-dimensional rational vectors, return a `N×N` matrix whose columns
are independent vectors from the list.
If no such `N` vectors exist, the returned matrix will not be invertible.
"""
function find_ratbasis(list::AbstractVector{<:StaticVector{N,<:Rational{T}}}) where {N,T}
    U = widen(T)
    n = length(list)
    ratbasis = zero(SizedMatrix{N,N,Rational{U}})
    i = 1
    @inbounds while i <= n
        l = list[i]
        if !iszero(l)
            ratbasis[:,1] .= l
            i += 1
            break
        end
        i += 1
    end
    j = 2
    while i <= n && j <= N
        l = list[i]
        if iszero(l)
            i += 1
            continue
        end
        ratbasis[:,j] .= l
        if rank(@view ratbasis[:,1:j]) == j
            j += 1
        end
        i += 1
    end
    return ratbasis
end

# Specialized version for N = 3
function find_ratbasis(list::AbstractVector{<:StaticVector{3,<:Rational{T}}}) where T
    U = widen(T)
    n = length(list)
    ratbasis = zero(SizedMatrix{3,3,Rational{U}})
    # We start by computing a basis for the i first vectors of v where i (≤ n) is minimal
    i = 1
    @inbounds while i <= n
        if !iszero(list[i])
            ratbasis[:,3] .= list[i]
            i += 1
            break
        end
        i += 1
    end
    @inbounds while i <= n
        v = list[i]
        valid = false
        # v is valid if it is linearly independant of ratbasis[:,3]
        if iszero(v[1])
            valid = (!iszero(ratbasis[1,3]) & (!iszero(v[2]) | !iszero(v[3]))) ||
                    (v[2] * ratbasis[3,3] != v[3] * ratbasis[2,3])
        else
            x = ratbasis[1,3] // v[1]
            valid = ratbasis[2,3] != x * v[2] || ratbasis[3,3] != x * v[3]
        end
        if valid
            ratbasis[:,2] .= v
            i += 1
            break
        end
        i += 1
    end
    d = zero(Rational{U})
    @inbounds while i <= n
        v = list[i]
        ratbasis[:,1] .= v
        d = det(ratbasis)
        iszero(d) || break
        i += 1
    end
    # At this point ratbasis is a basis in which the i-1 first coordinates are
    # [0, 0, 1] followed by a (possibly empty) series of [0, 0, x], then
    # [0, 1, 0] followed by a (possibly empty) series of [0, x, y], then [1, 0, 0]
    @inbounds if d < 0
        ratbasis[:,1] .= .-ratbasis[:,1] # Ensure a direct transformation
    end
    return SMatrix{3,3,Rational{U}}(ratbasis)
end

# Specialized version for N = 2
function find_ratbasis(list::AbstractVector{<:StaticVector{2,<:Rational{T}}}) where T
    U = widen(T)
    n = length(list)
    ratbasis = zero(SizedMatrix{2,2,Rational{U}})
    i = 1
    @inbounds while i <= n
        if !iszero(list[i])
            ratbasis[:,1] .= list[i]
            i += 1
            break
        end
        i += 1
    end
    d = zero(Rational{U})
    @inbounds while i <= n
        v = list[i]
        ratbasis[:,2] .= v
        d = det(ratbasis)
        iszero(d) || break
        i += 1
    end
    @inbounds if d < 0
        ratbasis[:,1] .= .-ratbasis[:,1]
    end
    return SMatrix{2,2,Rational{U}}(ratbasis)
end

"""
    normal_basis_rational(list::AbstractVector{<:StaticVector{N,<:Rational{T}}}) where {N,T}

Given a list of `N`-dimensional rational vectors, return a basis for the space
spanned by integer combinations of these vectors.

This basis is deterministically computed from the input.

It should depend only on the spanned space, not on the exact input (although this
assertion should be considered experimental for now).
"""
function normal_basis_rational(list::AbstractVector{<:StaticVector{N,<:Rational{T}}}) where {N,T}
    U = widen(T)
    n = length(list)
    ratbasis = find_ratbasis(list)
    # ratbasis should always be invertible in our setting: if not, it means that
    # we missed the fact that the dimensionality of the graph was strictly lower
    # than N, which should have been detected much earlier.
    # In this case, the definition of invratbasis will cause a failure.

    lcms = SizedVector{N,T}(ones(T, N))
    # per-column lcm of the vectors expressed in ratbasis
    expressed_in_ratbasis = Vector{SVector{N,Rational{T}}}(undef, n)
    invratbasis = inv(ratbasis) # See comment above in case of error here.
    @inbounds for j in 1:n
        expressed_in_ratbasis[j] = invratbasis * list[j]
        lcms .= lcm.(lcms, denominator.(expressed_in_ratbasis[j]))
    end

    expressed_in_basis = [MVector{N,Int}(undef) for _ in 1:n]
    @inbounds for j in 1:n
        x = expressed_in_ratbasis[j]
        expressed_in_basis[j] .= div.(numerator.(x).*lcms, denominator.(x))
    end
    #= lcms contains the lcm of all denominators of the coordinates of translations
       expressed in ratbasis.
       We now work in an artificial basis in which the coordinates are all integers.
       The integer coordinates are stored in expressed_in_basis.
       expressed_in_basis is now going to be reduced to its normal form.
    =#

    intbasis, _ = PeriodicGraphs.normal_basis(expressed_in_basis)

    # We finally go back from the integer normal form to a rational basis that expresses
    # all the original vectors with integer coordinates
    return SMatrix{N,N,Rational{U}}(ratbasis * [@inbounds(intbasis[i,j] // lcms[i]) for i in 1:N, j in 1:N])

    # newcoords = [SVector{N,Int}(normalization * x) for x in intcoords]
    # @assert (newbasis,) .* newcoords == list

    # return newbasis, newcoords
end

# Specialized version for N = 1
function normal_basis_rational(list::AbstractVector{<:StaticVector{1,<:Rational{T}}}) where T
    nzl = [x[] for x in list if !iszero(x[])]
    den = lcm(denominator.(nzl))
    S = widen(T)
    num = minimum(abs(S(den*x)) for x in nzl)
    return SMatrix{1,1,Rational{S}}(num // den)
end
