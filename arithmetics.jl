import LinearAlgebra
import LinearAlgebra: UpperTriangular, istril
import Base.GMP.MPZ



# Algorithm 3 from George Havas, Bohdan S. Majewski, and Keith R. Matthews,
# "Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction"

function reduce1!(k, i, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ)
    q = @inbounds begin
        if !iszero(a[i])
            div(a[k], a[i], RoundNearest)
        elseif 2*λ[i,k] > D[i+1]
            div(λ[i,k], D[i+1], RoundNearest)
        else
            zero(BigInt)
        end
    end

    @inbounds if !iszero(q)
        MPZ.mul!(tmpa, q, a[i])
        MPZ.sub!(a[k], tmpa)
        @simd for j in 1:m
            MPZ.mul!(tmpB[j], q, B[j,i])
            MPZ.sub!(B[j,k], tmpB[j])
        end
        MPZ.mul!(tmpD, q, D[i+1])
        MPZ.sub!(λ[i,k], tmpD)
        @simd for j in 1:i-1
            MPZ.mul!(tmpλ[j], q, λ[j,i])
            MPZ.sub!(λ[j,k], tmpλ[j])
        end
    end
    nothing
end

function swap!(k, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ, tmpt)
    a[k], a[k-1] = a[k-1], a[k]
    B[:,k], B[:,k-1] = B[:,k-1], B[:,k]
    @simd for j in 1:k-2
        λ[j,k], λ[j,k-1] = λ[j,k-1], λ[j,k]
    end
    @simd for i in k+1:m
        MPZ.mul!(tmpB[i], λ[k-1,i], D[k+1])
        MPZ.mul!(tmpλ[i], λ[k,i], λ[k-1,k])
        MPZ.sub!(tmpt[i], tmpB[i], tmpλ[i])
        MPZ.mul!(tmpB[i], λ[k-1,i], λ[k-1,k])
        MPZ.mul!(tmpλ[i], λ[k,i], D[k-1])
        MPZ.add!(tmpλ[i], tmpB[i])
        MPZ.cdiv_q!(λ[k-1,i], tmpλ[i], D[k])
        MPZ.cdiv_q!(λ[k,i], tmpt[i], D[k])
    end
    MPZ.mul!(tmpD, D[k-1], D[k+1])
    MPZ.pow_ui!(tmpa, λ[k-1,k], UInt(2))
    MPZ.add!(tmpD, tmpa)
    MPZ.cdiv_q!(D[k], tmpD, D[k])
    nothing
end

function extended_gcd(s)
    m = length(s)
    B = BigInt[i==j for i in 1:m, j in 1:m] # LinearAlgebra.I defaults to fill
    λ = BigInt[0 for _ in 1:m, _ in 1:m] # zeros defaults to fill
    # λ = zeros(Int, m, m)
    D = BigInt[1 for _ in 1:m+1]
    a = BigInt.(abs.(s))

    tmpa = BigInt()
    tmpD = BigInt()
    tmpB = [BigInt() for _ in 1:m]
    tmpλ = [BigInt() for _ in 1:m]
    tmpt = [BigInt() for _ in 1:m]

    k = 2
    @inbounds while k <= m
        reduce1!(k, k-1, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ)
        if !iszero(a[k-1]) || ((iszero(a[k-1]) & iszero(a[k])) &&
                               4*(D[k-1]*D[k+1] + λ[k-1,k]^2) < 3*D[k]^2)
            swap!(k, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ, tmpt)
            k -= (k > 2)
        else
            for i in k-2:-1:1
                reduce1!(k, i, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ)
            end
            k += 1
        end
    end
    if a[end] < 0
        MPZ.neg!(a[end])
        foreach(MPZ.neg!, @view B[:,end])
    end
    i = 1
    for x in s
        if x < 0
            MPZ.neg!(B[i,end])
        end
        i+=1
    end
    return a[end], B[:,end]
end



# Normal form for a list of vectors in ℚ³

function nf3D(list::AbstractVector{<:AbstractVector{<:Rational{T}}}) where T
    U = widen(T)
    n = length(list)
    ratbasis = zero(SizedMatrix{3,3,Rational{U}})
    # We start by computing a basis for the i first vectors of v where i (≤ n) is minimal
    i = 1
    while i <= n
        if !iszero(list[i])
            ratbasis[:,3] .= list[i]
            i += 1
            break
        end
        i += 1
    end
    while i <= n
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
    while i <= n
        v = list[i]
        ratbasis[:,1] .= v
        d = det(ratbasis)
        iszero(d) || break
        i += 1
    end
    @assert !iszero(d) # equivalent to @assert !iszero(det(ratbasis))
    # At this point ratbasis is a basis in which the i-1 first coordinates are
    # [0, 0, 1] followed by a (possibly empty) series of [0, 0, x], then
    # [0, 1, 0] followed by a (possibly empty) series of [0, x, y], then [1, 0, 0]
    if d < 0
        ratbasis[:,1] .= .-ratbasis[:,1] # Ensure a direct transformation
    end

    lcms = SizedVector{3,T}([1,1,1])
    # per-column lcm of the vectors expressed in ratbasis
    expressed_in_ratbasis = Vector{SVector{3,Rational{T}}}(undef, n)
    invratbasis = inv(SMatrix{3,3,Rational{U}}(ratbasis))
    for j in 1:n
        expressed_in_ratbasis[j] = invratbasis * list[j]
        lcms .= lcm.(lcms, denominator.(expressed_in_ratbasis[j]))
    end

    expressed_in_basis = [MVector{3,Int}(undef) for _ in 1:n]
    for j in 1:n
        x = expressed_in_ratbasis[j]
        expressed_in_basis[j] .= div.(numerator.(x).*lcms, denominator.(x))
    end
    #= lcms contains the lcm of all denominators of the coordinates of translations
       expressed in ratbasis.
       We now work in an artificial basis in which the coordinates are all integers.
       The integer coordinates are stored in expressed_in_basis.
       expressed_in_basis is now going to be reduced to its normal form.
    =#

    # intcoords = [MVector{3,Int}(undef) for _ in 1:n]
    intbasis = zero(MMatrix{3, 3, Int}) # The artificial integer basis
    for j in 1:3
        d, coefs = Tuple{Int, Vector{Int}}(extended_gcd(x[j] for x in expressed_in_basis))
        for i in 1:n
            intbasis[:,j] .+= coefs[i] .* expressed_in_basis[i]
        end
        @assert intbasis[j,j] > 0
        for i in 1:n
            coord = expressed_in_basis[i][j] .÷ d
            # intcoords[i][j] = coord
            expressed_in_basis[i] .-= coord .* intbasis[:,j]
        end
    end

    # normalization = MMatrix{3,3,Int}(LinearAlgebra.I)
    for (i,j) in ((2,1), (3,1), (3,2))
        # normalization[i,j] = (intbasis[i,j] ÷ intbasis[i,i]) - signbit(intbasis[i,j])
        # intbasis[:,j] .-= normalization[i,j] * intbasis[:,i]
        x = signbit(intbasis[i,j])
        intbasis[:,j] .-= (((intbasis[i,j] + x) ÷ intbasis[i,i]) - x) * intbasis[:,i]
    end

    # We finally go back from the integer normal form to a rational basis that expresses
    # all the original vectors with integer coordinates
    return SMatrix{3,3,Rational{U}}(ratbasis * [intbasis[i,j] // lcms[i] for i in 1:3, j in 1:3])

    # newcoords = [SVector{3,Int}(normalization * x) for x in intcoords]
    # @assert (newbasis,) .* newcoords == list

    # return newbasis, newcoords
    nothing
end
