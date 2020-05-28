include("systre.jl")
using Random
using Test
using .CIFTypes.PeriodicGraphs
using Base.Threads

function stress_test(g::Union{AbstractString,PeriodicGraph3D}, num=nothing,
                     refgraph=systre(PeriodicGraph3D(g)), skiplong=false)
    ref = string(refgraph)
    control = 50
    m = if num isa Nothing
        clamp(round(Int, begin
            _counter = 0
            timebefore = time_ns()
            while time_ns() - timebefore < 1.5e9
                _ref = string(systre(refgraph))
                @assert _ref == ref
                _counter += 1
            end
            timeafter = time_ns()
            elapsed = (timeafter - timebefore)*1.0e-9 / _counter
            if elapsed < 0.5
                20 + 30 / elapsed
            elseif elapsed < 5
                control = 10
                10 + 19 / elapsed
            elseif elapsed < 30
                control = 3
                3 + 59 / elapsed
            else
                skiplong && throw("(skipped: over 30s)")
                control = 1
                120 / elapsed
            end
        end), 2, 200)
    else
        num
    end
    println("(m = ", m, ')')
    n = nv(refgraph)
    # graph = PeriodicGraph3D(refgraph)
    # counter = 0
    graphs = [PeriodicGraph3D(refgraph) for _ in 1:nthreads()]
    counter = Atomic{Int}(0)
    #Juno.progress() do progressid
        @threads for i in 1:m
        # for i in 1:m
            id = threadid()
            # r = collect(1:n)
            # r[end], r[1], r[10] = r[10], r[end], r[1]
            r = randperm(n)
            # offsets = [zero(SVector{3,Int}) for _ in 1:n]
            offsets = [SVector{3,Int}([rand(-3:3) for _ in 1:3]) for _ in 1:n]
            graph = swap_axes!(offset_representatives!(graphs[id][r], offsets), randperm(3))
            graphs[id] = graph
            res = string(systre(graph))
            if res != ref
                #@info "stress_test main loop" progress="done" _id=progressid
                throw(string("TEST FAILED: different representations for res and ref:\n",
                        "ref = \"", ref, "\"\nres = \"", res, '\"', '\n'))
            end
            iterations = atomic_add!(counter, 1)
            # counter += 1
            #if iterations % (1 + m÷50) == 0
            #    @info "stress_test main loop" progress=iterations/m _id=progressid
            #end
            if iterations % control == 0
            # if counter % control == 0
                GC.gc(true)
                # @info "stress_test main loop" progress=counter/m _id=progressid
            end
        end
        #@info "stress_test main loop" progress="done" _id=progressid
    #end
    # foreach(println, graphs)
    nothing
end
function stress_test(net::CrystalNet, m=nothing)
    stress_test(net.graph, m)
end


function test_arc(path, newarc)
    arc = parse_arc(path)
    @show length(arc)
    @show nthreads()
    # id, _ = popfirst!(arc)
    # while id != "ucn"
    #     id, _ = popfirst!(arc)
    # end
    for (id, key) in arc
        key[1] != '3' && continue
        print("Testing ", id, ' ')
        try
            graph = try
                systre(PeriodicGraph3D(key))
            catch e
                if e == "UNTESTED"
                    println("(untested)")
                    continue
                end
                rethrow()
            end
            # open(newarc, "a") do f
            #     println(f, "key\t\t", string(graph))
            #     println(f, "id\t\t", id)
            #     println(f)
            # end
            stress_test(key, nothing, graph, false)
        catch e
            if e isa String
                println("\"\"\"", e, "\"\"\"")
            elseif e isa InterruptException || e isa TaskFailedException && e.task.exception isa InterruptException
                println("INTERRUPTION CAUGHT")
                rethrow()
            elseif e isa TaskFailedException
                println("TaskFailedException: ", e.task.exception, '\n')
            elseif e isa AssertionError
                println("AssertionError: ", e, '\n')
            else
                println("\nUNKNOWN ERROR happened:")
                showerror(stdout, e)
                println('\n')
            end
        end
        GC.gc(true)
        ccall(:malloc_trim, Cvoid, (Cint,), 0)
    end
    open(newarc, "a") do f
        println(f, '\n')
    end
end

test_arc("/home/liozou/Documents/MOFscripts/RCSRnets-2019-06-01.arc", "nets.arc")
