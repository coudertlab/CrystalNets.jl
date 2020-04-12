
include("./types.jl")
using .CIFTypes
using HTTP


function get_symmetries(i)
    s = ""
    while isempty(s)
        r = try
            HTTP.request("GET", "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen?what=text&gnum=$i")
        catch e
            e isa InterruptException && rethrow()
            continue
        end
        r.status != 200 && continue
        s = String(r.body)
    end
    @assert s[1:51] == "<h3>Full set of general positions</h3><hr><br><big>"
    @assert s[end-25:end] == "</big><br>\n</body>\n</html>"
    s = split(s[52:end-26], "</big><br><big>")
    @assert popfirst!(s) == "1 x,y,z"
    equivalents = EquivalentPosition[]
    for i in 2:length(s)+1
        x = s[i-1]
        @assert startswith(x, string(i)*' ')
        push!(equivalents, parse(EquivalentPosition, x[length(string(i))+1:end]))
    end
    return equivalents
end

function get_all_symmetries()
    x = Vector{EquivalentPosition}[]
    @progress for i in 1:230
        push!(x, get_symmetries(i))
    end
    return x
end
