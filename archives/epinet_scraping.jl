using HTTP, Gumbo, Cascadia, StaticArrays, PeriodicGraphs, CrystalNets
using Base.Threads
CrystalNets.toggle_warning(false)
CrystalNets.toggle_export(false)

## Database scraping

const numnets = let numpage, numtables
    numpage = parsehtml(String(HTTP.get("http://epinet.anu.edu.au/searches/1").body))
    numtable = only(eachmatch(sel"div[id=\"search-area\"]", numpage.root))
    parse(Int, only(match(r"matched all (\d+) nets", numtable.children[2].text).captures))
end

function _get_page(getlock, page)
    try
        lock(getlock)
        return HTTP.get(page)
    catch
        return nothing
    finally
        unlock(getlock)
    end
    nothing
end

function get_page(getlock, i)::String
    page = "http://epinet.anu.edu.au/sqc$i"
    while true
        http = _get_page(getlock, page)
        if http isa Nothing
            sleep(0.1)
        else
            return String((http::HTTP.Messages.Response).body)
        end
    end
    error("!")
    nothing
end

function scrape_epinet(numnets)
    nets = Vector{Pair{String,String}}(undef, numnets)
    getlock = ReentrantLock()
    _tmppath = joinpath(dirname(dirname(pathof(CrystalNets))), "archives", "epinet.TMP")
    todo = collect(1:numnets)
    tokeep = trues(numnets)
    for i in 1:nthreads()
        tmppath = _tmppath*string(i)*".arc"
        if isfile(tmppath)
            prev_arc = CrystalNets.parse_arc(tmppath)[2]
            for (v, k) in prev_arc
                @assert startswith(k, "sqc")
                l = parse(Int, k[4:end])
                nets[l] = v => k
                tokeep[l] = false
            end
        else
            open(tmppath, "w") do f
                println(f, "Made by CrystalNets.jl v$(CrystalNets.CRYSTAL_NETS_VERSION)\n")
            end
        end
    end
    todo = todo[tokeep]
    @threads for i in todo
        name = "sqc$i"
        try
            page = parsehtml(get_page(getlock, i))
            table = only(eachmatch(sel"div[id=\"systre_key\"]", page.root))
            _edges = table.children[1].children[1].children[2].children[1].children[1].children[1]
            _cvector = table.children[1].children[1].children[2].children[3].children[1].children[1]
            edges = Vector{PeriodicEdge3D}(undef, length(_edges.children))
            for (j, (edg, cvec)) in enumerate(zip(_edges.children, _cvector.children))
                s = parse(Int, (only(edg.children[1].children).text)::String)
                d = parse(Int, (only(edg.children[2].children).text)::String)
                ofs = SVector{3,Int}(parse(Int, (only(cvec.children[k].children).text)::String) for k in 1:3)
                edges[j] = PeriodicEdge3D(s, d, ofs)
            end
            g = PeriodicGraph3D(edges)
            genome = topological_genome(g)
            nets[i] = genome => name
            open(_tmppath*string(threadid())*".arc", "a") do f
                println(f, "key\t", genome, "\nid\t", name, '\n')
            end
            nothing
        catch e
            println("FAILED FOR $name with $e")
        end
        nothing
    end

    path = joinpath(dirname(dirname(pathof(CrystalNets))), "archives", "epinet.arc")
    arc = Dict{String,String}([(k,v) for (k,v) in nets])
    CrystalNets.export_arc(path, false, arc)
    for i in 1:nthreads()
        rm(_tmppath*string(i)*".arc")
    end
    arc
end

scrape_epinet(numnets)


## List of "known" nets in EPINET

function known_nets()
    page = parsehtml(String(HTTP.get("http://epinet.anu.edu.au/known_nets/show_all").body))
    ret = Int[]
    for node in eachmatch(sel"a[title=\"Network details\"]", page.root)
        t = only(node.children).text
        @assert startswith(t, "sqc")
        push!(ret, parse(Int, t[4:end]))
    end
    return ret
end
