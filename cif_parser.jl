using LightGraphs, GraphPlot

struct AtomID
    id::Int
end

struct Cell
    latticesystem::Symbol
    spacegroup::String
    tablenumber::Int
    a::Float64
    b::Float64
    c::Float64
    α::Float64
    β::Float64
    γ::Float64
end

struct EquivalentPositions end


struct CIF
    natoms::Int
    geometry::Cell
    atoms::Vector{AtomID}
    bonds::BitMatrix
    pos::Matrix{Float64}
end

struct Crystal
    natoms::Int
    atoms::Vector{AtomID}
    bonds::BitMatrix
    pos::Matrix{Float64}
end

function nextword(l, i)
    n = lastindex(l)
    i == n && return (0, 0)

    i = nextind(l, i)
    start = 0
    if l[i] == '\'' || l[i] == '"'
        if l[nextind(l, i)] == '#'
            i = findnext(isequal('\n'), l, i)
            i == nothing && return (0, 0)
        else
            @assert isspace(l[nextind(l, i)])
            i = nextind(l, nextind(l, i))
        end
    end
    while i <= n
        if isspace(l[i])
            i = nextind(l, i)
        else
            if l[i] == '#'
                i = findnext(isequal('\n'), l, i)
                i == nothing && return (0, 0)
            else
                start = i
                break
            end
        end
    end
    start == 0 && return (0, 0)

    inquote = false
    quotesymb = "'"
    while i <= n
        c = l[i]
        if !isspace(c)
            c == '#' && return (start, prevind(l, i))
            if c == '\'' || c == '"'
                inquote = true
                quotesymb = c
                break
            end
            i = nextind(l, i)
        else
            return (start, prevind(l, i))
        end
    end
    if inquote
        while i < n
            c = l[i]
            if c == quotesymb && (isspace(l[nextind(l, i)]) ||
                                 l[nextind(l, i)] == '#')
                return (start+1, prevind(l, i))
            end
            i = nextind(l, i)
        end
        if !l[n] == quotesymb
            throw("Invalid syntax: opening quote $quotesymb at position $start is not closed")
        end
        return (start+1, prevind(l, i))
    end
    return (start, n)
end

function main(file)
    all_data = Dict{String, Any}()
    inloop = false
    loopisspecified = false
    loopspec = String[]
    loop_n = 0

    l = read(file, String)
    i, j = nextword(l, 0)
    while i != 0

        if inloop
            if !loopisspecified
                if l[i] != '_' # This indicates the start of the values
                    loop_n = length(loopspec)
                    loopisspecified = true
                    @show loopspec
                else # The identifier is part of the loop specification
                    push!(loopspec, l[i+1:j])
                    all_data[l[i+1:j]] = Any[]
                    i, j = nextword(l, j); continue
                end
            end

            # From this point, the loop has been specified
            @assert loopisspecified
            if l[i] != '_' && l[i:j] != "loop_"
                for k in 1:loop_n
                    push!(all_data[loopspec[k]], l[i:j])
                    i, j = nextword(l, j)
                end
                continue
            end
            if l[i:j] == "loop_"
                loopisspecified = false
                loopspec = String[]
                i, j = nextword(l, j); continue
            end

            # This point can only be reached if we just quitted a loop
            inloop = false
        end

        @assert !inloop
        if l[i] == '_' # Simple identifier definition
            next_i, next_j = nextword(l, j)
            @assert next_i != 0
            all_data[l[i+1:j]] = l[next_i:next_j]
            j = next_j
        else
            if l[i:j] == "loop_"
                inloop = true
                loopisspecified = false
                loopspec = String[]
            elseif j-i > 4 && l[i:5] == "data_"
                # simply skip
            else
                k = findprev(isequal('\n'), l, i)
                n = count("\n", l[1:k])
                throw("Unkown word \"$(l[i:j])\" at line $n, position $(i-k):$(j-k)")
            end
        end

        i, j = nextword(l, j)
    end

    return all_data
end
