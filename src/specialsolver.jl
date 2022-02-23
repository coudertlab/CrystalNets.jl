## Specialized solver for sparse symmetric integer linear systems with exact rational solution

include("Modulos.jl")
using .Modulos
import Base.GMP: MPZ

import LinearAlgebra: BlasInt, checknonsingular, LU, tril!, triu!, ipiv2perm,
                      lu, lu!, ldiv!, Factorization, issuccess
using SparseArrays
import SparseArrays: getcolptr

using BigRationals


function rational_lu!(B::SparseMatrixCSC, col_offset, check=true)
    Tf = eltype(B)
    m, n = size(B)
    minmn = min(m, n)
    info = 0
    @inbounds begin
        for k in 1:minmn
            ipiv = getcolptr(B)[k] + col_offset[k]
            piv = nonzeros(B)[ipiv]
            if iszero(piv)
                check && checknonsingular(k-1, Val(false)) # TODO update with Pivot
                return LU{Tf,SparseMatrixCSC{Tf,Int}}(B, collect(1:minmn), convert(BlasInt, k-1))
            end
            Bkkinv = inv(piv)
            @simd for i in ipiv+1:getcolptr(B)[k+1]-1
                nonzeros(B)[i] *= Bkkinv
            end
            for j in k+1:n
                r1 = getcolptr(B)[j]
                r2 = getcolptr(B)[j+1]-1
                r = searchsortedfirst(rowvals(B), k, r1, r2, Base.Forward)
                ((r > r2) || (rowvals(B)[r] != k)) && continue
                Bkj = nonzeros(B)[r]
                for i in ipiv+1:getcolptr(B)[k+1]-1
                    Bik = nonzeros(B)[i]
                    l = i - ipiv

                    while rowvals(B)[l+r] < rowvals(B)[i]
                        r += 1
                    end
                    nonzeros(B)[l+r] -= Bik * Bkj
                end
            end
        end
    end
    check && checknonsingular(info, Val(false))
    return LU{Tf,SparseMatrixCSC{Tf,Int}}(B, Vector{BlasInt}(1:minmn), convert(BlasInt, info))
end

# function lu!(B::SparseMatrixCSC{<:Rational}, ::Val{Pivot} = Val(false);
#                            col_offset, check::Bool = true) where Pivot
function rational_lu!(B::SparseMatrixCSC{BigRational}, col_offset, check::Bool=true)
    Tf = Rational{BigInt}
    m, n = size(B)
    minmn = min(m, n)
    info = 0
    Bkkinv = BigRational()
    tmp = BigRational()
    @inbounds begin
        for k in 1:minmn
            ipiv = getcolptr(B)[k] + col_offset[k]
            piv = nonzeros(B)[ipiv]
            if iszero(piv)
                check && checknonsingular(k-1, Val(false)) # TODO update with Pivot
                return LU{Tf,SparseMatrixCSC{Tf,Int}}(Tf.(B), Vector{BlasInt}(1:minmn), convert(BlasInt, k-1))
            end
            BigRationals.MPQ.inv!(Bkkinv, piv)
            @simd for i in ipiv+1:getcolptr(B)[k+1]-1
                BigRationals.MPQ.mul!(nonzeros(B)[i], Bkkinv)
            end
            for j in k+1:n
                r1 = getcolptr(B)[j]
                r2 = getcolptr(B)[j+1]-1
                r = searchsortedfirst(rowvals(B), k, r1, r2, Base.Forward)
                ((r > r2) || (rowvals(B)[r] != k)) && continue
                Bkj = nonzeros(B)[r]
                for i in ipiv+1:getcolptr(B)[k+1]-1
                    Bik = nonzeros(B)[i]
                    l = i - ipiv

                    while rowvals(B)[l+r] < rowvals(B)[i]
                        r += 1
                    end
                    # Base.GMP.MPZ.mul!(tmp, Bik, Bkj)
                    # Base.GMP.MPZ.sub!(nonzeros(B)[l+r], tmp)
                    BigRationals.MPQ.mul!(tmp, Bik, Bkj)
                    BigRationals.MPQ.sub!(nonzeros(B)[l+r], tmp)
                end
            end
        end
    end
    check && checknonsingular(info, Val(false))
    return LU{Tf,SparseMatrixCSC{Tf,Int}}(Tf.(B), Vector{BlasInt}(1:minmn), convert(BlasInt, info))
end

# function lu(A::SparseMatrixCSC{<:Rational}, pivot::Union{Val{false}, Val{true}} = Val(false); check::Bool = true)
function rational_lu(A::SparseMatrixCSC, check::Bool=true, ::Type{Ti}=BigRational) where {Ti}
    Tf = Ti == BigRational ? Rational{BigInt} : Ti

    Base.require_one_based_indexing(A)
    _I, _J, _V = findnz(A)
    I, J, V = issorted(_J) ? (_I, _J, _V) : begin
        _indices = sortperm(_J)
        @inbounds (_I[_indices], _J[_indices], _V[_indices])
    end
    # @inbounds if !issorted(_J)
    #     indices = sortperm(J)
    #     I = I[indices]; J = J[indices]; V = V[indices]
    # end
    isempty(J) && return LU{Tf,SparseMatrixCSC{Tf,Int}}(Tf.(A), Int[], convert(BlasInt, 0))
    m, n = size(A)
    minmn = min(m, n)
    if J[1] != 1 || I[1] != 1
        check && checknonsingular(1, Val(false)) # TODO update with Pivot
        # return LU{eltype(A), typeof(A)}(A, collect(1:minmn), convert(BlasInt, 1))
        return LU{Tf,SparseMatrixCSC{Tf,Int}}(Tf.(A), collect(1:minmn), convert(BlasInt, 1))
    end

    col_offset = zeros(Int, minmn) # for each col, index of the pivot element
    idx_cols = [[I[i] for i in getcolptr(A)[col+1]-1:-1:getcolptr(A)[col]] for col in 1:minmn]
    # For each column, indices of the non-zeros elements
    in_idx_colscol = falses(n)
    for col in 2:minmn
        sort!(idx_cols[col-1]; rev=true)
        # All idx_cols[x] are sorted by decreasing order for x < col
        # @show idx_cols[col]
        idx_colscol = idx_cols[col]
        in_idx_colscol[idx_colscol] .= true
        for row_j in idx_colscol
            row_j >= col && continue
            col_offset[col] += 1
            # @show idx_cols[row_j]
            idx_colsj = idx_cols[row_j]
            sizcol = length(idx_colscol)
            for row_i in idx_colsj
                if row_i ≤ row_j
                    break # Because the row_i are sorted in decreasing order
                end
                if !in_idx_colscol[row_i]
                    push!(idx_colscol, row_i)
                    in_idx_colscol[row_i] = true
                end
            end
            countadd = length(idx_colscol) - sizcol
            if countadd > 0
                siz = length(I)
                resize!(I, siz + countadd)
                resize!(J, siz + countadd)
                resize!(V, siz + countadd)
                for i in 1:countadd
                    row_i = idx_colscol[sizcol+i]
                    _idx = siz + i
                    J[_idx] = col
                    I[_idx] = row_i
                    V[_idx] = 0
                end
            end
        end
        in_idx_colscol[idx_colscol] .= false
    end
    B = sparse(I, J, Ti.(V)) # TODO update with Pivot
    # lu!(B, col_offset, check)
    rational_lu!(B, col_offset, check)
end

#=
function lu(A::Hermitian{T, <:SparseMatrixCSC{T}}, pivot::Union{Val{false}, Val{true}} = Val(false); check::Bool = true) where {T<:Rational, Pivot}
    lu(ishermitian(A.data) ? A.data : sparse(A), pivot)
end

function Base.getproperty(F::LU{T,<:SparseMatrixCSC{<:Rational{<:Integer}, <:Integer}}, d::Symbol) where T
    m, n = size(F)
    if d === :L
        L = tril!(getfield(F, :factors)[1:m, 1:min(m,n)])
        for i = 1:min(m,n); L[i,i] = one(T); end
        return L
    elseif d === :U
        return triu!(getfield(F, :factors)[1:min(m,n), 1:n])
    elseif d === :p
        return ipiv2perm(getfield(F, :ipiv), m)
    elseif d === :P
        return Matrix{T}(LinearAlgebra.I, m, m)[:,invperm(F.p)]
    else
        getfield(F, d)
    end
end
=#

function forward_substitution!(L::SparseMatrixCSC, b)
    _, n = size(L)
    _, m = size(b)
    @inbounds for col in 1:n
        k = getcolptr(L)[col]
        rowvals(L)[k] == col || checknonsingular(col, Val(true))
        invnzLk = inv(nonzeros(L)[k])
        x = invnzLk .* b[col,:]
        b[col,:] .= x
        #=@simd=# for i in (k+1):getcolptr(L)[col+1]-1
            nzLi = nonzeros(L)[i]
            rvLi = rowvals(L)[i]
            @simd for j in 1:m
                b[rvLi,j] -= nzLi*x[j]
            end
        end
    end
    nothing
end

function backward_substitution!(U::SparseMatrixCSC, b)
    _, n = size(U)
    _, m = size(b)
    @inbounds for col in n:-1:1
        k = getcolptr(U)[col+1]-1
        rowvals(U)[k] == col || checknonsingular(col, Val(true))
        invnzUk = inv(nonzeros(U)[k])
        x = invnzUk .* b[col,:]
        b[col,:] .= x
        #=@simd=# for i in getcolptr(U)[col]:(k-1)
            nzUi = nonzeros(U)[i]
            rvUi = rowvals(U)[i]
            @simd for j in 1:m
                b[rvUi,j] -= nzUi*x[j]
            end
        end
    end
    nothing
end

function linsolve!(F::LU{<:Any,<:AbstractSparseMatrix}, B::Base.StridedVecOrMat)
    TFB = typeof(oneunit(eltype(B)) / oneunit(eltype(F)))
    BB = similar(B, TFB, size(B))
    copyto!(BB, B)
    m, n = size(F)
    minmn = min(m,n)
    L = tril!(getfield(F, :factors)[1:m, 1:minmn])
    for i = 1:minmn; L[i,i] = 1; end
    forward_substitution!(L, BB)
    x = triu!(getfield(F, :factors)[1:minmn, 1:n])
    backward_substitution!(x, BB)
    return BB
end

#=
function ldiv!(F::LU{<:Any,<:AbstractSparseMatrix}, B::Base.StridedVecOrMat)
    forward_substitution!(F.L, B)
    backward_substitution!(F.U, B)
    return B
end
=#

function rational_solve(::Val{N}, A, Y) where N
    B = rational_lu(A, false)
    if !issuccess(B)
        error("Singular exception while equilibrating. Is the graph connected?")
    end
    Z = linsolve!(B, Rational{BigInt}.(Y))
    return hcat(zeros(Rational{Int128}, N), adjoint(Rational{Int128}.(Z)))
    # Rational{Int64} is not enough for tep for instance.
end



# function _inner_dixon_p!(Z::Matrix{Rational{T}}, h, x̄, sqh, tmp) where T
#     for j in eachindex(Z)
#         ua = MPZ.set(h)
#         ub = @inbounds x̄[j]
#         va = zero(T)
#         vb = one(T)
#         k = 0
#         while ub >= sqh
#             k += 1
#             # cpua = deepcopy(ua)
#             # cpub = deepcopy(ub)
#             MPZ.tdiv_qr!(tmp, ua, ua, ub)
#             ua, ub = ub, ua
#             # @toggleassert tmp == cpua ÷ cpub
#             # @toggleassert ua == cpub
#             # @toggleassert ub == cpua - tmp * cpub
#             # cpuc = deepcopy(va)
#             if typemin(Clong) < vb < typemax(Clong)
#                 MPZ.mul_si!(tmp, vb % Clong)
#             else
#                 tmp *= vb
#             end
#             flag = signbit(va)
#             va = abs(va)
#             if va < typemax(Culong)
#                 if flag
#                     MPZ.sub_ui!(tmp, va)
#                 else
#                     MPZ.add_ui!(tmp, va)
#                 end
#                 va, vb = vb, T(tmp)
#             else
#                 va, vb = vb, va + tmp
#             end
#             # @toggleassert vb == cpuc + tmp * va
#         end
#         uc, vc = if T === BigInt
#             Base.divgcd(ub, vb)
#         else
#             ud, vd = Base.divgcd(ub, vb)
#             if T !== BigInt
#                 m = typemin(T)
#                 M = typemax(T)
#                 (m < ud < M && m < vd < M) || return false
#             end
#             (ud % T, vd % T)
#         end
#
#         @inbounds Z[j] = (-2*isodd(k)+1) * Base.checked_den(uc, vc)
#         # @show Z[j]
#         # @toggleassert mod((-1)^isodd(k) * ub, h) == mod(vb * x̄[j], h)
#     end
#     return true
# end

function copyuntil(j, oldZ, ::Type{T}) where T
    Z = similar(oldZ, T)
    for i in eachindex(Z)
        i == j && return Z
        Z[i] = oldZ[i]
    end
    error("Invalid failure of _inner_dixon_p!, please report this error")
    return Z # Does not matter
end


function _inner_dixon_p!(indices, Z::Matrix{Rational{T}}, h, x̄, sqh, tmp) where T
    while !isempty(indices)
        j = pop!(indices)
        ua = MPZ.set(h)
        ub = deepcopy(@inbounds x̄[j])
        va = Int128(0)
        vb = Int128(1)
        k = 0
        while ub >= sqh
            k += 1
            # cpua = deepcopy(ua)
            # cpub = deepcopy(ub)
            MPZ.tdiv_qr!(tmp, ua, ua, ub)
            ua, ub = ub, ua
            # @toggleassert tmp == cpua ÷ cpub
            # @toggleassert ua == cpub
            # @toggleassert ub == cpua - tmp * cpub
            # cpuc = deepcopy(va)
            if typemin(Clong) < vb < typemax(Clong)
                MPZ.mul_si!(tmp, vb % Clong)
            else
                tmp *= vb
            end
            flag = signbit(va)
            va = abs(va)
            if va < typemax(Culong)
                if flag
                    MPZ.sub_ui!(tmp, va)
                else
                    MPZ.add_ui!(tmp, va)
                end
                va, vb = vb, Int128(tmp)
            else
                va, vb = vb, va + tmp
            end
            # @toggleassert vb == cpuc + tmp * va
        end

        uv::Tuple{T,T} = if T === BigInt
            Base.divgcd(ub, vb)
        else
            ud, vd = Base.divgcd(ub, vb)
            m = typemin(T)
            M = typemax(T)
            if !(m < ud < M && m < vd < M)
                push!(indices, j)
                return false
            end
            (ud % T, vd % T)
        end

        @inbounds Z[j] = (-2*isodd(k)+1) * Base.checked_den(uv[1], uv[2])
        # @show Z[j]
        # @toggleassert mod((-1)^isodd(k) * ub, h) == mod(vb * x̄[j], h)
    end
    return true
end

function dixon_p(::Val{N}, A, C::Factorization{Modulo{p,T}}, Y) where {N,p,T}
    λs = [norm(x) for x in eachcol(A)]
    append!(λs, [norm(x) for x in eachcol(Y)])
    partialsort!(λs, N)
    for _ in 1:N
        popfirst!(λs)
    end
    δ::BigFloat = prod(BigFloat, λs; init=one(BigFloat))
    # @show δ
    # @show p
    m = ceil(Int, 2*log(δ / (MathConstants.φ - 1))/log(p))
    @toggleassert m ≥ 1
    # @show m
    B = copy(Y)
    x̄ = BigInt.(linsolve!(C, B))
    X = copy(x̄)
    @toggleassert A * Modulo{p,T}.(X) == B
    h = one(BigInt) # = p^i
    tmp = BigInt()
    for i in 1:m-1
        MPZ.mul_si!(h, p)
        B .= (B .- A*Integer.(X)) .÷ p
        X .= Integer.(linsolve!(C, B))
        @toggleassert A * Modulo{p,T}.(X) == B
        # x̄ .+= h .* X
        @inbounds for j in eachindex(x̄)
            MPZ.mul!(tmp, X[j], h)
            MPZ.add!(x̄[j], tmp)
        end
    end
    MPZ.mul_si!(h, p) # h = p^m
    @toggleassert mod.(A * x̄, h) == mod.(Y, h)
    sqh = MPZ.sqrt(h) # h = p^{m/2}
    typeofZ = Union{Matrix{Rational{Int64}},Matrix{Rational{Int128}},Matrix{Rational{BigInt}}}
    Z::typeofZ = similar(Y, Rational{Int64})
    indices = collect(reverse(eachindex(Z)))
    success = _inner_dixon_p!(indices, Z, h, x̄, sqh, tmp)
    if !success
        Z = copyuntil(first(indices), Z, Rational{Int128})
        success = _inner_dixon_p!(indices, Z, h, x̄, sqh, tmp)
        if !success
            Z = copyuntil(first(indices), Z, Rational{BigInt})
            success = _inner_dixon_p!(indices, Z, h, x̄, sqh, tmp)
            @toggleassert success
        end
    end

    @toggleassert eltype(Y).(A * big.(Z)) == Y
    return Z::typeofZ
    # return hcat(zeros(Rational{Int128}, N), Rational{Int128}.(x̄)')
    # Rational{Int64} is not enough for tep for instance.
end

"""
    dixon_solve(::Val{N}, A, Y) where N

Specialized solver for the linear system `A*X = Y` where `A` is a sparse symmetric
integer `n×n` matrix and `Y` is a dense integer `n×N` matrix, using Dixon's method.

Return `X` as a matrix of `Rational{Int128}`.
"""
function dixon_solve(::Val{N}, A, Y) where N
    # @show time_ns()
    typeofB = Union{
        LU{Modulo{2147483647,Int32},SparseMatrixCSC{Modulo{2147483647,Int32},Int}},
        LU{Modulo{2147483629,Int32},SparseMatrixCSC{Modulo{2147483629,Int32},Int}},
        LU{Modulo{2147483587,Int32},SparseMatrixCSC{Modulo{2147483587,Int32},Int}}
    }
    B::typeofB = rational_lu(A, false, Modulo{2147483647,Int32})
    typeofZ = Union{Matrix{Rational{Int64}},Matrix{Rational{Int128}},Matrix{Rational{BigInt}}}
    if issuccess(B)
        Z = dixon_p(Val(N), A, B, Y)'
    else
        B = rational_lu(A, false, Modulo{2147483629,Int32})
        if issuccess(B)
            Z = dixon_p(Val(N), A, B, Y)'
        else
            B = rational_lu(A, false, Modulo{2147483587,Int32})
            if issuccess(B)
                Z = dixon_p(Val(N), A, B, Y)'
            else # The probability of this being required is *extremely* low
                return rational_solve(Val(N), A, Y)
            end
        end
    end
    return hcat(zeros(eltype(Z), N), Z)::typeofZ
end
