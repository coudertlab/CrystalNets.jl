using CrystalNets, PeriodicGraphs, ArgParse, LinearAlgebra, SparseArrays,
      StaticArrays, Logging, Tokenize, BigRationals
import Chemfiles

const __bodyfunction__ = Dict{Method,Any}()

# Find keyword "body functions" (the function that contains the body
# as written by the developer, called after all missing keyword-arguments
# have been assigned values), in a manner that doesn't depend on
# gensymmed names.
# `mnokw` is the method that gets called when you invoke it without
# supplying any keywords.
function __lookup_kwbody__(mnokw::Method)
    function getsym(arg)
        isa(arg, Symbol) && return arg
        @assert isa(arg, GlobalRef)
        return arg.name
    end

    f = get(__bodyfunction__, mnokw, nothing)
    if f === nothing
        fmod = mnokw.module
        # The lowered code for `mnokw` should look like
        #   %1 = mkw(kwvalues..., #self#, args...)
        #        return %1
        # where `mkw` is the name of the "active" keyword body-function.
        ast = Base.uncompressed_ast(mnokw)
        if isa(ast, Core.CodeInfo) && length(ast.code) >= 2
            callexpr = ast.code[end-1]
            if isa(callexpr, Expr) && callexpr.head == :call
                fsym = callexpr.args[1]
                if isa(fsym, Symbol)
                    f = getfield(fmod, fsym)
                elseif isa(fsym, GlobalRef)
                    if fsym.mod === Core && fsym.name === :_apply
                        f = getfield(mnokw.module, getsym(callexpr.args[2]))
                    elseif fsym.mod === Core && fsym.name === :_apply_iterate
                        f = getfield(mnokw.module, getsym(callexpr.args[3]))
                    else
                        f = getfield(fsym.mod, fsym.name)
                    end
                else
                    f = missing
                end
            else
                f = missing
            end
        else
            f = missing
        end
        __bodyfunction__[mnokw] = f
    end
    return f
end

function precompile_kwarg(@nospecialize(tt), @nospecialize(kwargtypes))
    let fbody = try __lookup_kwbody__(which(tt)) catch missing end
        if !ismissing(fbody)
            ttu = Base.unwrap_unionall(tt)
            newtt = Tuple{Core.Typeof(fbody), kwargtypes..., ttu.parameters...}
            precompile(Base.rewrap_unionall(newtt, tt))
        end
    end
end

function _precompile_dependencies()
    # Tokenize
    precompile(Tuple{typeof(Base._collect),UnitRange{Int},Tokenize.Lexers.Lexer{IOBuffer, Tokenize.Tokens.Token},Base.HasEltype,Base.SizeUnknown})

    # Graphs
    precompile(Tuple{typeof(Graphs.floyd_warshall_shortest_paths),Graphs.SimpleGraph{Int},Graphs.DefaultDistance})
    precompile(Tuple{Type{Graphs.SimpleGraph},Vector{Graphs.SimpleGraphs.SimpleEdge{Int}}})

    # PeriodicGraphs
    precompile(Tuple{Type{Dict{PeriodicEdge3D, Nothing}}})
    precompile(Tuple{Type{Dict{PeriodicVertex3D, Nothing}}})
    precompile(Tuple{typeof(Base._unique!),typeof(identity),Vector{PeriodicEdge3D},Set{PeriodicEdge3D},Int,Int})
    precompile(Tuple{typeof(Base.ht_keyindex),Dict{PeriodicVertex3D, Nothing},PeriodicVertex3D})
    precompile(Tuple{typeof(copyto!),Vector{PeriodicEdge3D},PeriodicGraphs.PeriodicEdgeIter{3}})
    precompile(Tuple{typeof(deleteat!),Vector{PeriodicVertex3D},BitVector})
    precompile(Tuple{typeof(has_edge),PeriodicGraph3D,PeriodicEdge3D})
    precompile(Tuple{typeof(offset_representatives!),PeriodicGraph3D,Vector{SVector{3,Int}}})
    precompile(Tuple{typeof(searchsortedfirst),Vector{PeriodicVertex3D},PeriodicVertex3D,Int,Int,Base.Order.ForwardOrdering})
    precompile(Tuple{typeof(setindex!),Dict{PeriodicEdge3D, Nothing},Nothing,PeriodicEdge3D})
    precompile(Tuple{typeof(sort!),Vector{PeriodicEdge3D},Int,Int,Base.Sort.MergeSortAlg,Base.Order.ForwardOrdering,Vector{PeriodicEdge3D}})
    precompile(Tuple{typeof(sort!),Vector{PeriodicVertex3D},Int,Int,Base.Sort.MergeSortAlg,Base.Order.ForwardOrdering,Vector{PeriodicVertex3D}})
    precompile(Tuple{typeof(union!),Set{PeriodicVertex3D},Vector{PeriodicVertex3D}})
    precompile(Tuple{typeof(Base.setindex!),Dict{PeriodicEdge3D, Nothing},Nothing,PeriodicEdge3D})
    precompile(Tuple{typeof(Base.union!),Set{PeriodicVertex3D},Vector{PeriodicVertex3D}})

    # SparseArrays
    precompile(Tuple{typeof(*),SparseArrays.SparseMatrixCSC{Int, Int},Matrix{Rational{BigInt}}})
    precompile(Tuple{typeof(Base.copyto_unaliased!),IndexCartesian,SubArray{Int, 1, SparseArrays.SparseMatrixCSC{Int, Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, false},IndexLinear,Vector{Int}})
    precompile(Tuple{typeof(Base.mightalias),SubArray{Int, 1, SparseArrays.SparseMatrixCSC{Int, Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, false},Vector{Int}})
    precompile(Tuple{typeof(LinearAlgebra.mul!),Matrix{BigInt},SparseArrays.SparseMatrixCSC{Int, Int},Matrix{BigInt},Bool,Bool})
    precompile(Tuple{typeof(SparseArrays.dimlub),Vector{Int}})
    precompile(Tuple{typeof(SparseArrays.findnz),SparseArrays.SparseMatrixCSC{Int, Int}})
    precompile(Tuple{typeof(SparseArrays.sparse_check_length),String,Vector{Int},Int,Type})
    precompile(Tuple{typeof(SparseArrays.spzeros),Type{Int},Type{Int},Int,Int})
    precompile(Tuple{typeof(getindex),SparseArrays.SparseMatrixCSC{Int, Int},UnitRange{Int},UnitRange{Int}})

    # Logging
    precompile(Tuple{typeof(Base.CoreLogging.handle_message),Logging.ConsoleLogger,Any,Any,Any,Any,Any,Any,Any})
    precompile(Tuple{typeof(Base.CoreLogging.shouldlog),Logging.ConsoleLogger,Base.CoreLogging.LogLevel,Module,Symbol,Symbol})
    precompile(Tuple{typeof(Logging.default_metafmt),Base.CoreLogging.LogLevel,Any,Any,Any,Any,Any})
    precompile(Tuple{typeof(Logging.termlength),SubString{String}})
    let fbody = try __lookup_kwbody__(which(Base.CoreLogging.handle_message, (Logging.ConsoleLogger,Any,Any,Any,Any,Any,Any,Any,))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Any,typeof(Base.CoreLogging.handle_message),Logging.ConsoleLogger,Any,Any,Any,Any,Any,Any,Any,))
        end
    end

    # ArgParse
    precompile(Tuple{Core.kwftype(typeof(ArgParse.Type)),Any,Type{ArgParse.ArgParseSettings}})
    precompile(Tuple{Core.kwftype(typeof(ArgParse.add_arg_field!)),Any,typeof(ArgParse.add_arg_field!),ArgParse.ArgParseSettings,Vector{T} where T<:AbstractString})
    precompile(Tuple{typeof(ArgParse.parse1_optarg!),ArgParse.ParserState,ArgParse.ArgParseSettings,ArgParse.ArgParseField,Any,AbstractString})
    precompile(Tuple{typeof(ArgParse.preparse!),Channel,ArgParse.ParserState,ArgParse.ArgParseSettings})
    precompile(Tuple{typeof(ArgParse.print_group),IO,Vector{T} where T,AbstractString,Int,Int,AbstractString,AbstractString,AbstractString})
    let fbody = try __lookup_kwbody__(which(ArgParse.show_help, (ArgParse.ArgParseSettings,))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Any,typeof(ArgParse.show_help),ArgParse.ArgParseSettings,))
        end
    end
    let fbody = try __lookup_kwbody__(which(any, (Function,Vector{ArgParse.ArgParseGroup},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Function,typeof(any),Function,Vector{ArgParse.ArgParseGroup},))
        end
    end

    # LinearAlgebra
    precompile(Tuple{Type{Matrix{Float64}},LinearAlgebra.UniformScaling{Bool},Tuple{Int, Int}})
    precompile(Tuple{typeof(LinearAlgebra.norm),Vector{Float64},Int})
    precompile(Tuple{typeof(eltype),LinearAlgebra.Adjoint{Rational{Int}, Matrix{Rational{Int}}}})
    precompile(Tuple{typeof(isone),Matrix{Int32}})
    precompile(Tuple{typeof(hcat),Vector{Rational{Int}},LinearAlgebra.Adjoint{Rational{Int}, Matrix{Rational{Int}}}})
    let fbody = try __lookup_kwbody__(which(LinearAlgebra.rank, (Matrix{Int},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Float64,Float64,typeof(LinearAlgebra.rank),Matrix{Int},))
        end
    end

    # Base
    precompile(Tuple{Core.kwftype(typeof(Base.with_output_color)),NamedTuple{(:bold,), Tuple{Bool}},typeof(Base.with_output_color),Function,Symbol,IOContext{Base.TTY},String,Vararg{Any, N} where N})
    precompile(Tuple{Type{Base.IteratorSize},Base.Iterators.ProductIterator{Tuple{UnitRange{Int}, UnitRange{Int}}}})
    precompile(Tuple{Type{Dict{Symbol, Any}},Base.Iterators.Pairs{Symbol, Any, Tuple{Symbol, Symbol, Symbol}, NamedTuple{(:help, :metavar, :required), Tuple{String, String, Bool}}}})
    precompile(Tuple{Type{Dict{Symbol, Any}},Base.Iterators.Pairs{Symbol, Any, Tuple{Symbol, Symbol}, NamedTuple{(:help, :action), Tuple{String, Symbol}}}})
    precompile(Tuple{Type{Dict{Symbol, Any}},Base.Iterators.Pairs{Symbol, String, Tuple{Symbol, Symbol}, NamedTuple{(:help, :metavar), Tuple{String, String}}}})
    precompile(Tuple{Type{Tuple{Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}, Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}, Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}}},Vector{Tuple{Vector{Any}, Vector{Any}}}})
    for T in (Int32, Int64, Int128)
        precompile(Tuple{Type{SubArray},IndexLinear,Matrix{Rational{T}},Tuple{Base.Slice{Base.OneTo{Int}}, Int},Tuple{Bool}})
        precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, Type{Rational{T}}, Tuple{Matrix{Rational{Int}}}}})
    end
    precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(big), Tuple{Matrix{Rational{Int}}}}})
    precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(denominator), Tuple{Matrix{Rational{Int}}}}})
    precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(numerator), Tuple{Matrix{Rational{Int}}}}})
    precompile(Tuple{typeof(Base.deepcopy_internal),NTuple{9, BigFloat},IdDict{Any, Any}})
    precompile(Tuple{typeof(Base.display_error),Base.TTY,ErrorException,Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}})
    precompile(Tuple{typeof(Base.grow_to!),Dict{Symbol, Any},Tuple{Pair{Symbol, String}, Pair{Symbol, Bool}, Pair{Symbol, String}},Int})
    precompile(Tuple{typeof(Base.grow_to!),Dict{Symbol, Any},Tuple{Pair{Symbol, String}, Pair{Symbol, Symbol}},Int})
    precompile(Tuple{typeof(Base.print_to_string),Int128,Vararg{Any, N} where N})
    precompile(Tuple{typeof(Base.reducedim_init),Function,typeof(min),Matrix{Float64},Int})
    precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Int},Expr,Int})
    precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{Any}},Tuple{},Int})
    precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{Int}},Tuple{},Int})
    precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{}},Tuple{Int},Int})
    precompile(Tuple{typeof(Base.typed_hvcat),Type{BigFloat},Tuple{Int, Int, Int},BigFloat,Vararg{Number, N} where N})
    precompile(Tuple{typeof(Base.typed_hvcat),Type{Float64},Tuple{Int, Int, Int},Int,Vararg{Number, N} where N})
    precompile(Tuple{typeof(copy),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Tuple{Base.OneTo{Int}, Base.OneTo{Int}}, Type{Int}, Tuple{Matrix{Rational{BigInt}}}}})
    precompile(Tuple{typeof(maximum),Matrix{Int}})
    precompile(Tuple{typeof(minimum),Matrix{Int}})
    precompile(Tuple{typeof(push!),Vector{String},String,String,String})
    precompile(Tuple{typeof(setindex!),Dict{String, Tuple{Exception, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}},Tuple{ArgumentError, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}},String})
    precompile(Tuple{typeof(string),Int128,String,Vararg{Any, N} where N})
    precompile(Tuple{typeof(vcat),Vector{Expr},Vector{Expr}})
    @static if VERSION >= v"1.6-"
        let fbody = try __lookup_kwbody__(which(Base.print_within_stacktrace, (IOContext{Base.TTY},String,Vararg{Any, N} where N,))) catch missing end
            if !ismissing(fbody)
                precompile(fbody, (Symbol,Bool,typeof(Base.print_within_stacktrace),IOContext{Base.TTY},String,Vararg{Any, N} where N,))
            end
        end
    end
    let fbody = try __lookup_kwbody__(which(all, (Function,Vector{String},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Function,typeof(all),Function,Vector{String},))
        end
    end
    let fbody = try __lookup_kwbody__(which(any, (Function,Vector{AbstractString},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Function,typeof(any),Function,Vector{AbstractString},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sort!, (Vector{Int},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sort!),Vector{Int},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{Symbol},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{Symbol},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{Tuple{Int, Vector{Int}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{Tuple{Int, Vector{Int}}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{Vector{Int}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{Vector{Int}},))
        end
    end

    # StaticArrays
    precompile(Tuple{Type{Vector{SVector{3, Float64}}},Vector{SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true}}})
    precompile(Tuple{Type{Vector{SVector{3, Int}}},Vector{Any}})
    precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{0},StaticArrays.StaticArrayStyle{2}})
    precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},StaticArrays.StaticArrayStyle{1}})
    precompile(Tuple{Type{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(+),Tuple{SVector{3, Int}, SVector{3, Int}}})
    precompile(Tuple{Type{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(-),Tuple{SVector{3, Rational{Int}}, SVector{3, Int}}})
    precompile(Tuple{Type{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(floor),Tuple{Base.RefValue{Type{Int}}, SVector{3, Rational{Int}}}})
    precompile(Tuple{Type{Dict{Int, SVector{3, Int}}}})
    precompile(Tuple{Type{Dict{SVector{3, Int}, Nothing}}})
    precompile(Tuple{Type{SMatrix{3, 3, BigFloat, 9}},NTuple{9, Float64}})
    precompile(Tuple{typeof(==),Vector{Tuple{Int, Int, SVector{3, Rational{Int}}}},Vector{Tuple{Int, Int, SVector{3, Rational{Int}}}}})
    precompile(Tuple{typeof(>),SVector{3, Int},SizedVector{3, Int, 1}})
    precompile(Tuple{typeof(>),SizedVector{3, Int, 1},SVector{3, Int}})
    precompile(Tuple{typeof(Base.Broadcast.broadcasted),Function,SVector{3, Rational{Int}},SVector{3, Int}})
    precompile(Tuple{typeof(StaticArrays.arithmetic_closure),Type{BigFloat}})
    precompile(Tuple{typeof(LinearAlgebra.generic_norm2),SubArray{BigFloat, 1, SMatrix{3, 3, BigFloat, 9}, Tuple{Base.Slice{SOneTo{3}}, Int}, true}})
    precompile(Tuple{typeof(StaticArrays._axes),Size{(3, 3)}})
    precompile(Tuple{typeof(StaticArrays._axes),Size{(3,)}})

    precompile(Tuple{Type{Size},Type{SubArray{Int, 1, Matrix{Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, true}}})
    precompile(Tuple{typeof(<),SVector{3, Int},SizedVector{3, Int, 1}})
    for T in (Int32, Int64, Int128, BigInt)
        precompile(Tuple{Type{Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}},Vector{Any}})
        precompile(Tuple{Type{SVector{3, Int}},Tuple{Rational{T}, Rational{T}, Rational{T}}})
        precompile(Tuple{Type{Dict{Int, Vector{SMatrix{3, 3, Rational{T}, 9}}}}})
        precompile(Tuple{Type{Dict{SVector{9, Rational{T}}, Nothing}}})
        precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(isempty),Tuple{Tuple{Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}}}})
        precompile(Tuple{typeof(<),SVector{3, Rational{T}},SVector{3, Rational{Int}}})
        precompile(Tuple{typeof(<),SVector{3, Rational{Int}},SVector{3, Rational{T}}})
        precompile(Tuple{typeof(==),SVector{3, Rational{Int}},SVector{3, Rational{T}}})
        precompile(Tuple{typeof(==),SMatrix{3, 3, T, 9},LinearAlgebra.UniformScaling{Bool}})
        precompile(Tuple{typeof(LinearAlgebra.generic_matmatmul!),Matrix{Rational{widen(T)}},Char,Char,SizedMatrix{3, 3, Rational{widen(T)}, 2},Matrix{Rational{T}},LinearAlgebra.MulAddMul{true, true, Bool, Bool}})
        precompile(Tuple{typeof(LinearAlgebra.dot),SubArray{BigFloat, 1, SMatrix{3, 3, BigFloat, 9}, Tuple{Base.Slice{SOneTo{3}}, T}, true},SubArray{BigFloat, 1, SMatrix{3, 3, BigFloat, 9}, Tuple{Base.Slice{SOneTo{3}}, T}, true}})
        precompile(Tuple{typeof(Base.Broadcast.broadcasted),Base.Broadcast.Style{Tuple},Function,Tuple{Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}}})
        precompile(Tuple{typeof(Base.Broadcast.materialize!),Base.Broadcast.DefaultArrayStyle{1},SubArray{Rational{T}, 1, Matrix{Rational{T}}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true},Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(-), Tuple{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(+), Tuple{SVector{3, Rational{T}}, SVector{3, Int}}}, SVector{3, Rational{T}}}}})
        precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(-), Tuple{SVector{3, Rational{T}}, SVector{3, Int}}}})
        precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(floor), Tuple{Base.RefValue{Type{Int}}, SVector{3, Rational{T}}}}})
        precompile(Tuple{typeof(append!),Vector{SMatrix{3, 3, Rational{T}, 9}},Vector{SMatrix{3, 3, Rational{T}, 9}}})
        precompile(Tuple{typeof(append!),Vector{SMatrix{3, 3, Rational{T}, 9}},Vector{SMatrix{3, 3, Rational{widen(T)}, 9}}})
        let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Rational{T}}},))) catch missing end
            if !ismissing(fbody)
                precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Rational{T}}},))
            end
        end
        let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Rational{Int128}}},))) catch missing end
            if !ismissing(fbody)
                precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Rational{T}}},))
            end
        end
        precompile(Tuple{typeof(sort!),Vector{Int},Base.Sort.QuickSortAlg,Base.Order.Perm{Base.Order.ForwardOrdering, Vector{SVector{3, Rational{T}}}}})
    end
    precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(//), Tuple{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(round), Tuple{Base.RefValue{Type{Int128}}, Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(*), Tuple{Int128, SVector{3, Float64}}}}}, Int128}}})
    precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(//), Tuple{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(round), Tuple{Base.RefValue{Type{BigInt}}, Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(*), Tuple{BigInt, SVector{3, Float64}}}}}, BigInt}}})
    precompile(Tuple{typeof(Base._foldl_impl),Base.BottomRF{typeof(hcat)},Base.ReshapedArray{Int, 2, SVector{3, Int}, Tuple{}},Set{SVector{3, Int}}})
    precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{StaticArrays.Dynamic}},Tuple{Int},Int})
    precompile(Tuple{typeof(StaticArrays._axes),Size{(9,)}})
    precompile(Tuple{typeof(setindex!),Dict{Int, SVector{3, Int}},SVector{3, Int},Int})
    precompile(Tuple{typeof(setindex!),Dict{SVector{3, Int}, Nothing},Nothing,SVector{3, Int}})
    precompile(Tuple{typeof(which(StaticArrays._getindex,(AbstractArray,Tuple{Vararg{Size, N} where N},Any,)).generator.gen),Any,Any,Any,Any})
    precompile(Tuple{typeof(which(StaticArrays._map,(Any,Vararg{AbstractArray, N} where N,)).generator.gen),Any,Any,Any})
    precompile(Tuple{typeof(which(StaticArrays.combine_sizes,(Tuple{Vararg{Size, N} where N},)).generator.gen),Any,Any})
    let fbody = try __lookup_kwbody__(which(LinearAlgebra.rank, (Base.ReshapedArray{Int, 2, SVector{3, Int}, Tuple{}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Float64,Float64,typeof(LinearAlgebra.rank),Base.ReshapedArray{Int, 2, SVector{3, Int}, Tuple{}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Float64}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Float64}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Rational{BigInt}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Rational{BigInt}}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Rational{Bool}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Rational{Bool}}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sprint, (Function,SVector{3, Int},Vararg{Any, N} where N,))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Nothing,Int,typeof(sprint),Function,SVector{3, Int},Vararg{Any, N} where N,))
        end
    end
    precompile(Tuple{typeof(sort!),Vector{Int},Base.Sort.QuickSortAlg,Base.Order.Perm{Base.Order.ForwardOrdering, Vector{SVector{3, Float64}}}})


    # Chemfiles
    precompile(Tuple{Type{Chemfiles.Atom},Chemfiles.Frame,Int})
    precompile(Tuple{Type{Chemfiles.Atom},String})
    precompile(Tuple{Type{Chemfiles.Frame}})
    precompile(Tuple{Type{Chemfiles.Topology},Chemfiles.Frame})
    precompile(Tuple{Type{Chemfiles.UnitCell},BigFloat,BigFloat,BigFloat,BigFloat,BigFloat,BigFloat})
    precompile(Tuple{Type{Chemfiles.UnitCell},Chemfiles.Frame})
    precompile(Tuple{typeof(Chemfiles.__free),Chemfiles.CxxPointer{Chemfiles.lib.CHFL_ATOM}})
    precompile(Tuple{typeof(Chemfiles.__free),Chemfiles.CxxPointer{Chemfiles.lib.CHFL_CELL}})
    precompile(Tuple{typeof(Chemfiles.__free),Chemfiles.CxxPointer{Chemfiles.lib.CHFL_TOPOLOGY}})
    precompile(Tuple{typeof(Chemfiles.__strip_null),String})
    precompile(Tuple{typeof(Chemfiles.add_atom!),Chemfiles.Frame,Chemfiles.Atom,Vector{Float64},Vector{Float64}})
    precompile(Tuple{typeof(Chemfiles.bonds),Chemfiles.Topology})
    precompile(Tuple{typeof(Chemfiles.matrix),Chemfiles.UnitCell})
    precompile(Tuple{typeof(Chemfiles.positions),Chemfiles.Frame})
    precompile(Tuple{typeof(Chemfiles.set_type!),Chemfiles.Atom,String})

        # CrystalNets
    # precompile(Tuple{CrystalNets.var"#730#threadsfor_fun#121"{Tuple{Symbol}, String, Dict{String, String}, Dict{String, Tuple{Exception, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}}, Dict{String, String}, Vector{String}}})
    for T in (Int32, Int64, Int128, BigInt)
        # precompile(Tuple{CrystalNets.var"#670#threadsfor_fun#114"{CrystalNet{Rational{T}}, Vector{Int}, Base.Threads.SpinLock, Vector{Pair{Int, Tuple{Matrix{Rational{T}}, Vector{Int}}}}, Int, DataType, Vector{Int}}})
        # precompile(Tuple{CrystalNets.var"#685#threadsfor_fun#116"{Rational{T}, Dict{Int, Vector{SMatrix{3, 3, Rational{T}, 9}}}, Base.Threads.SpinLock, DataType, Vector{Pair{Int, Tuple{Matrix{Rational{T}}, Vector{Int}}}}}})

        precompile(Tuple{Type{CrystalNet{Rational{T}}},CrystalNets.Cell,Vector{Symbol},PeriodicGraph3D,Matrix{Rational{T}}})
        precompile(Tuple{typeof(CrystalNets.isrank3),Matrix{Rational{T}}})
        precompile(Tuple{typeof(CrystalNets.minimize),CrystalNet{Rational{T}}})
        precompile(Tuple{typeof(CrystalNets.topological_key),CrystalNet{Rational{T}}})
        precompile(Tuple{typeof(topological_genome),CrystalNet{Rational{T}}})

    end
    precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(CrystalNets.parsestrip), Tuple{Vector{String}}}})
    precompile(Tuple{typeof(CrystalNets.invalid_input_error),String,ErrorException,Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}})
    precompile(Tuple{typeof(CrystalNets.julia_main)})

    precompile(Tuple{typeof(CrystalNets.parse_chemfile),String})
    precompile(Tuple{typeof(CrystalNets.parse_chemfile),String,Bool})
    precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(CrystalNets.parsestrip),Tuple{Vector{String}}})
    precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(CrystalNets.parsestrip), Tuple{Vector{String}}}})
    precompile(Tuple{typeof(Base.deepcopy_internal),Vector{CrystalNets.EquivalentPosition},IdDict{Any, Any}})
    precompile(Tuple{typeof(CrystalNets.isrank3),Matrix{Rational{Int32}}})
    precompile(Tuple{typeof(parse_chemfile),String})
    precompile(Tuple{typeof(recognize_topology),PeriodicGraph3D})
    for T in (Nothing, CrystalNets.Clusters)
        precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}})
        precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, Nothing})
        precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, String})
        precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, Nothing, String})
        precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, String, String})
    end

    precompile(Tuple{typeof(convert), Type{Any}, Vector{Union{String,Vector{String}}}})
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    # _precompile_dependencies()

    inttypes = (Int32, Int64, Int128, BigInt)
    primes = (2147483647, 2147483629, 2147483587)

    graph = PeriodicGraph3D
    types = Vector{Symbol}
    opts = CrystalNets.Options
    cryst = CrystalNets.Crystal{Nothing}
    crystclust = CrystalNets.Crystal{CrystalNets.Clusters}
    clust = CrystalNets.Clusters
    _pos = SVector{3,Float64}
    p_pos = Vector{SVector{3,Float64}}
    pos{D,T} = SVector{D,Rational{T}}
    ppos{D,T} = Vector{SVector{D,Rational{T}}}
    ofs{D} = SVector{D,Int}
    pofs{D} = Vector{SVector{D,Int}}
    cnet{D,T} = CrystalNet{D,Rational{T}}
    mat = SMatrix{3,3,BigFloat,9}
    fmat = SMatrix{3,3,Float64,9}
    sbus = Vector{Vector{PeriodicVertex3D}}
    kinds = CrystalNets.ClusterKinds
    structure = StructureType._StructureType
    clustering = Clustering._Clustering
    neighs = Vector{PeriodicVertex3D}
    edgs = Vector{PeriodicEdge3D}
    cell = CrystalNets.Cell
    cif = CrystalNets.CIF
    bondlist = Vector{Vector{Tuple{Int,Float32}}}
    unets = CrystalNets.UnderlyingNets
    genome = CrystalNets.TopologicalGenome
    result = CrystalNets.TopologyResult
    collisions = Vector{CrystalNets.CollisionNode}
    modulo{P} = CrystalNets.Modulo{P, Int32}
    collision = CrystalNets.CollisionNode
    smat{D,T,L} = SMatrix{D,D,Rational{T},L}

    # kwargs = Base.Pairs{Symbol,V,Tuple{Vararg{Symbol,N}}, NamedTuple{names,T}} where {V,N,names,T<:Tuple{Vararg{Any,N}}}


    # Modulos.jl
    for P in primes
        precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Int,Function})
        precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Int,Function})
        precompile(Tuple{typeof(-),CrystalNets.Modulos.Modulo{P, Int32},CrystalNets.Modulos.Modulo{P, Int32}})
        precompile(Tuple{typeof(/),Int,CrystalNets.Modulos.Modulo{P, Int32}})
        precompile(Tuple{typeof(==),Matrix{CrystalNets.Modulos.Modulo{P, Int32}},Matrix{Int}})
        precompile(Tuple{typeof(Base._unsafe_copyto!),Matrix{CrystalNets.Modulos.Modulo{P, Int32}},Int,Matrix{Int},Int,Int})
        precompile(Tuple{typeof(Base._unsafe_getindex),IndexLinear,Matrix{CrystalNets.Modulos.Modulo{P, Int32}},Int,Base.Slice{Base.OneTo{Int}}})
        precompile(Tuple{typeof(Base.copyto_unaliased!),IndexLinear,SubArray{CrystalNets.Modulos.Modulo{P, Int32}, 1, Matrix{CrystalNets.Modulos.Modulo{P, Int32}}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, true},IndexLinear,Vector{CrystalNets.Modulos.Modulo{P, Int32}}})
        precompile(Tuple{typeof(LinearAlgebra.mul!),Matrix{CrystalNets.Modulos.Modulo{P, Int32}},SparseArrays.SparseMatrixCSC{Int, Int},Matrix{CrystalNets.Modulos.Modulo{P, Int32}},Bool,Bool})
        precompile(Tuple{typeof(SparseArrays._setindex_scalar!),SparseArrays.SparseMatrixCSC{CrystalNets.Modulos.Modulo{P, Int32}, Int},Int,Int,Int})
        precompile(Tuple{typeof(SparseArrays.sparse!),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Int,typeof(+),Vector{Int},Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}}})
        precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Int,Function})
        precompile(Tuple{typeof(SparseArrays.sparse_check_length),String,Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Type})
        precompile(Tuple{typeof(getindex),SparseArrays.SparseMatrixCSC{CrystalNets.Modulos.Modulo{P, Int32}, Int},UnitRange{Int},UnitRange{Int}})
    end


    # output.jl
    precompile(Tuple{typeof(CrystalNets.export_dataline), Base.IOStream, String})
    precompile(Tuple{typeof(CrystalNets.export_vtf), String, cryst, Int, Bool})
    precompile(Tuple{typeof(CrystalNets.export_vtf), String, crystclust, Int, Bool})
    precompile(Tuple{typeof(CrystalNets.export_vtf), String, cryst, Int})
    precompile(Tuple{typeof(CrystalNets.export_vtf), String, crystclust, Int})
    precompile(Tuple{typeof(CrystalNets.export_vtf), String, cryst})
    precompile(Tuple{typeof(CrystalNets.export_vtf), String, crystclust})
    for D in 1:3
        for T in inttypes
            precompile(Tuple{typeof(CrystalNets.export_vtf), String, cnet{D,T}, Int, Bool})
            precompile(Tuple{typeof(CrystalNets.export_vtf), String, cnet{D,T}, Int})
            precompile(Tuple{typeof(CrystalNets.export_vtf), String, cnet{D,T}})
        end
    end
    precompile(Tuple{typeof(CrystalNets.export_cif), String, cryst})
    precompile(Tuple{typeof(CrystalNets.export_cif), String, crystclust})
    precompile(Tuple{typeof(CrystalNets.export_cif), String, cif})
    precompile(Tuple{typeof(CrystalNets.export_cgd), String, cryst})
    precompile(Tuple{typeof(CrystalNets.export_cgd), String, crystclust})
    for D in 1:3
        precompile(Tuple{typeof(CrystalNets.export_cgd), String, PeriodicGraph{D}})
        for T in inttypes
            precompile(Tuple{typeof(CrystalNets.export_cgd), String, cnet{D,T}})
        end
    end
    precompile(Tuple{typeof(CrystalNets.export_attributions), crystclust, String})
    precompile(Tuple{typeof(CrystalNets.export_attributions), crystclust})
    precompile(Tuple{typeof(CrystalNets.export_arc), String, Bool, Dict{String,String}})
    precompile(Tuple{typeof(CrystalNets.export_arc), String, Bool})
    precompile(Tuple{typeof(CrystalNets.export_arc), String})

    # utils.jl
    for T in (Int8, Int16, Int32, Int64, Int128, BigInt)
        precompile(Tuple{typeof(CrystalNets.double_widen), Type{T}})
        precompile(Tuple{typeof(CrystalNets.double_widen), Type{Rational{T}}})
    end
    precompile(Tuple{typeof(CrystalNets.recursive_readdir!), Vector{String}, String, String})
    precompile(Tuple{typeof(CrystalNets.recursive_readdir), String})
    precompile(Tuple{typeof(CrystalNets.tmpexportname), String, String, String, String})
    precompile(Tuple{typeof(CrystalNets.tmpexportname), String, String, Nothing, String})
    for S in (Int, Nothing)
        for D in 1:3
            precompile_kwarg(Tuple{typeof(export_default), PeriodicGraph{D}, String, String, String}, (S,))
            @static if VERSION > v"1.8-"
                precompile_kwarg(Tuple{typeof(export_default), PeriodicGraph{D}, Base.LazyString, String, String}, (S,))
            end
            for T in inttypes
                precompile_kwarg(Tuple{typeof(export_default), cnet{D, T}, String, String, String}, (S,))
                @static if VERSION > v"1.8-"
                    precompile_kwarg(Tuple{typeof(export_default), cnet{D, T}, Base.LazyString, String, String}, (S,))
                end
            end
        end
        precompile_kwarg(Tuple{typeof(export_default), cryst, String, String, String}, (S,))
        precompile_kwarg(Tuple{typeof(export_default), crystclust, String, String, String}, (S,))
        @static if VERSION > v"1.8-"
            precompile_kwarg(Tuple{typeof(export_default), cryst, Base.LazyString, String, String}, (S,))
            precompile_kwarg(Tuple{typeof(export_default), crystclust, Base.LazyString, String, String}, (S,))
        end
    end
    # precompile_kwarg(Tuple{typeof(export_default), PeriodicGraph2D, String, String, String}, (kwargs,))
    # precompile_kwarg(Tuple{typeof(export_default), PeriodicGraph3D, String, String, String}, (kwargs,))
    precompile(Tuple{typeof(CrystalNets.string_atomtype), Symbol})
    precompile(Tuple{typeof(CrystalNets.representative_atom), Symbol, Int})
    for T in inttypes
        for D in 1:3
            precompile(Tuple{typeof(CrystalNets.issingular), smat{D,T,D*D}})
            precompile(Tuple{typeof(CrystalNets.angle), pos{D,T}, pos{D, T}})
            precompile(Tuple{typeof(CrystalNets.dihedral), pos{D,T},pos{D,T},pos{D,T}})
        end
        precompile(Tuple{typeof(CrystalNets.isrank3), Matrix{Rational{T}}})
        precompile(Tuple{typeof(CrystalNets.isrank2), Matrix{Rational{T}}})
        precompile(Tuple{typeof(CrystalNets.back_to_unit), Rational{T}})
    end
    precompile(Tuple{typeof(CrystalNets.nextword), String, Int})
    precompile(Tuple{typeof(CrystalNets.isinterrupt), Any})
    precompile(Tuple{typeof(CrystalNets.isoverfloworinexact), Any})

    # options.jl
    precompile(Tuple{typeof(CrystalNets.clustering_from_num), Int})
    precompile(Tuple{typeof(CrystalNets.clustering_from_num), Int8})
    precompile(Tuple{typeof(CrystalNets.clustering_from_symb), Symbol})
    precompile(Tuple{typeof(parse), Type{clustering}, String})
    precompile(Tuple{typeof(parse), Type{clustering}, SubString{String}})
    precompile(Tuple{typeof(getindex), kinds, Int})
    precompile(Tuple{typeof(getindex), kinds, Symbol})
    precompile(Tuple{typeof(CrystalNets.getmetal), kinds})
    precompile(Tuple{typeof(length), kinds})
    precompile(Tuple{typeof(CrystalNets.ifbooltempdirorempty), String})
    precompile(Tuple{typeof(CrystalNets.ifbooltempdirorempty), Bool})
    precompile(Tuple{Type{CrystalNets.Options}})
    precompile(Tuple{Type{CrystalNets.Options}, CrystalNets.Options})

    # specialsolver.jl
    for Ti in (BigRational, (modulo{P} for P in primes)...)
        precompile(Tuple{typeof(CrystalNets.rational_lu!), SparseMatrixCSC{Ti,Int}, Vector{Int}, Bool})
        precompile(Tuple{typeof(CrystalNets.rational_lu!), SparseMatrixCSC{Ti,Int}, Vector{Int}})
    end
    for P in primes
        precompile(Tuple{typeof(CrystalNets.rational_lu), SparseMatrixCSC{modulo{P},Int}, Bool, Type{modulo{P}}})
        precompile(Tuple{typeof(CrystalNets.rational_lu), SparseMatrixCSC{modulo{P},Int}, Bool, Type{BigRational}})
        precompile(Tuple{typeof(CrystalNets.rational_lu), SparseMatrixCSC{Int,Int}, Bool})
    end
    for Ti in (Rational{BigInt}, (modulo{P} for P in primes)...)
        precompile(Tuple{typeof(CrystalNets.forward_substitution!), SparseMatrixCSC{Ti,Int}, Matrix{Ti}})
        precompile(Tuple{typeof(CrystalNets.backward_substitution!), SparseMatrixCSC{Ti,Int}, Matrix{Ti}})
    end
    @static if VERSION < v"1.8-"
        precompile(Tuple{typeof(CrystalNets.linsolve!), LU{Rational{BigInt},SparseMatrixCSC{Rational{BigInt},Int}}, Matrix{Rational{BigInt}}})
        for P in primes
            precompile(Tuple{typeof(CrystalNets.linsolve!), LU{modulo{P},SparseMatrixCSC{modulo{P},Int}}, Matrix{Int}})
        end
    else
        precompile(Tuple{typeof(CrystalNets.linsolve!), LU{Rational{BigInt},SparseMatrixCSC{Rational{BigInt},Int},Base.OneTo{Int}}, Matrix{Rational{BigInt}}})
        for P in primes
            precompile(Tuple{typeof(CrystalNets.linsolve!), LU{modulo{P},SparseMatrixCSC{modulo{P},Int},Base.OneTo{Int}}, Matrix{Int}})
        end
    end
    for T in (Int64, Int128, BigInt)
        precompile(Tuple{typeof(CrystalNets.copyuntil), Int, Matrix{Rational{T}}, Type{Rational{T}}})
        precompile(Tuple{typeof(CrystalNets._inner_dixon_p!), Vector{Int}, Matrix{Rational{T}}, BigInt, Matrix{BigInt}, BigInt, BigInt})
    end
    for N in 1:3
        precompile(Tuple{typeof(CrystalNets.rational_solve), Val{N}, SparseMatrixCSC{Int,Int}, Matrix{Int}})
        for P in primes
            @static if VERSION < v"1.8-"
                precompile(Tuple{typeof(CrystalNets.dixon_p), Val{N}, SparseMatrixCSC{Int,Int}, LU{modulo{P},SparseMatrixCSC{modulo{P},Int}}, Matrix{Int}})
            else
                precompile(Tuple{typeof(CrystalNets.dixon_p), Val{N}, SparseMatrixCSC{Int,Int}, LU{modulo{P},SparseMatrixCSC{modulo{P},Int},Base.OneTo{Int}}, Matrix{Int}})
            end
        end
        precompile(Tuple{typeof(CrystalNets.dixon_solve), Val{N}, SparseMatrixCSC{Int,Int}, Matrix{Int}})
    end

    # types.jl 1/2
    precompile(Tuple{Type{cif}, Dict{String,Union{String,Vector{String}}}, cell, Vector{Int}, types, Matrix{Float64}, bondlist})
    precompile(Tuple{typeof(CrystalNets.keepinbonds), bondlist, Vector{Int}})
    precompile(Tuple{typeof(CrystalNets.add_to_bondlist!), Vector{Tuple{Int,Float32}}, Int, Float32})
    precompile(Tuple{typeof(CrystalNets.get_bondlist), Vector{Tuple{Int,Float32}}, Int})
    precompile(Tuple{typeof(CrystalNets.sortprune_bondlist!), Vector{Tuple{Int,Float32}}})
    precompile(Tuple{typeof(CrystalNets.remove_partial_occupancy), cif})
    precompile(Tuple{typeof(CrystalNets.prune_collisions), cif})
    precompile(Tuple{typeof(CrystalNets.expand_symmetry), cif})
    precompile(Tuple{typeof(CrystalNets.edges_from_bonds), bondlist, fmat, p_pos})
    precompile(Tuple{Type{clust}, sbus, Vector{Int}, Vector{Int}, pofs{3}, BitVector})
    precompile(Tuple{Type{clust}, Int})
    precompile(Tuple{typeof(isempty), clust})
    precompile(Tuple{typeof(getindex), clust, Vector{Int}})
    precompile(Tuple{typeof(getindex), clust, Base.OneTo{Int}})
    precompile(Tuple{typeof(getindex), clust, Base.UnitRange{Int}})
    precompile(Tuple{Type{cryst}, cell, types, p_pos, graph, opts})
    precompile(Tuple{Type{crystclust}, cell, types, clust, p_pos, graph, opts})
    precompile(Tuple{Type{CrystalNets.Crystal}, cell, types, Nothing, p_pos, graph, opts})
    precompile(Tuple{Type{CrystalNets.Crystal}, cell, types, clust, p_pos, graph, opts})
    precompile(Tuple{Type{CrystalNets.Crystal}})
    precompile(Tuple{typeof(==), cryst, cryst})
    precompile(Tuple{typeof(==), crystclust, crystclust})
    precompile(Tuple{Type{cryst}, cryst})
    precompile(Tuple{Type{cryst}, crystclust})
    precompile(Tuple{Type{crystclust}, cryst, clust})
    precompile(Tuple{Type{crystclust}, crystclust, clust})
    precompile(Tuple{Type{crystclust}, cryst})
    precompile(Tuple{Type{crystclust}, crystclust})
    precompile(Tuple{typeof(CrystalNets.remove_metal_cluster_bonds!), graph, types, opts})
    precompile(Tuple{typeof(CrystalNets.trimmed_crystal), cryst})
    precompile(Tuple{typeof(getindex), cryst, Vector{Int}})
    precompile(Tuple{typeof(getindex), cryst, Base.OneTo{Int}})
    precompile(Tuple{typeof(getindex), cryst, Base.UnitRange{Int}})
    precompile(Tuple{typeof(getindex), crystclust, Vector{Int}})
    precompile(Tuple{typeof(getindex), crystclust, Base.OneTo{Int}})
    precompile(Tuple{typeof(getindex), crystclust, Base.UnitRange{Int}})
    for D in 1:3
        precompile(Tuple{typeof(CrystalNets.equilibrium), PeriodicGraph{D}})
        precompile(Tuple{typeof(CrystalNets.trim_topology), PeriodicGraph{D}})
        for T in inttypes
            precompile(Tuple{Type{cnet{D,T}}, cell, types, ppos{D,T}, PeriodicGraph{D}, opts})
        end
    end
    for D in 1:3
        for T in inttypes
            for D2 in 1:3
                for T2 in inttypes
                    precompile(Tuple{Type{cnet{D,T}}, cnet{D2,T2}})
                end
                precompile(Tuple{Type{cnet{D}}, cnet{D2,T}})
            end
            precompile(Tuple{typeof(ndims), cnet{D,T}})
            precompile(Tuple{Type{cnet{D,T}}, cell, types, PeriodicGraph{D}, ppos{D,T}, opts})
            precompile(Tuple{Type{cnet{D,T}}, cell, opts})
            precompile(Tuple{typeof(show), Base.TTY, cnet{D,T}})
            precompile(Tuple{typeof(show), Base.IOStream, cnet{D,T}})
        end
        precompile(Tuple{Type{cnet{D}}, cell, opts})
    end

    # clustering.jl
    precompile(Tuple{typeof(CrystalNets.regroup_sbus), graph, Vector{Int}, Vector{Int}})
    precompile(Tuple{typeof(CrystalNets._trim_monovalent!), graph})
    precompile(Tuple{typeof(CrystalNets.trim_monovalent), cryst})
    precompile(Tuple{typeof(CrystalNets.trim_monovalent), crystclust})
    precompile(Tuple{typeof(CrystalNets.delete_target_from_list!), sbus, Int})
    precompile(Tuple{typeof(CrystalNets.is_paddlewheel_candidate!), Vector{Union{Missing, Tuple{Symbol, PeriodicVertex3D}}}, sbus, Int, types, Set{Int}})
    # precompile(Tuple{typeof(CrystalNets.bond_carboxylic_acid!), graph, types})
    precompile(Tuple{typeof(CrystalNets.regroup_paddlewheel!), graph, clust, types, Set{Int}})
    precompile(Tuple{typeof(CrystalNets.split_sbu!), clust, graph, Int, Vector{Int}})
    precompile(Tuple{typeof(CrystalNets.reclassify!), clust, Vector{Int}, Int, graph, types, Dict{Symbol,Int}, Int})
    precompile(Tuple{typeof(CrystalNets.add_to_newclass!), Vector{Int}, graph, clust, Int, Int, types, Set{Symbol}})
    precompile(Tuple{typeof(CrystalNets.add_to_merge_or_newclass!), Vector{Int}, Vector{Tuple{ofs{3},Int}}, graph, clust, Set{Int}, Int, Int})
    precompile(Tuple{typeof(CrystalNets.small_cycles_around), graph, p_pos, mat, Int, PeriodicVertex3D, Vector{Int}, Set{Int}})
    precompile(Tuple{typeof(CrystalNets.in_small_cycles_around), graph, p_pos, mat, Int, Vector{Int}, Set{Int}})
    precompile(Tuple{typeof(CrystalNets.detect_organiccycles), Vector{Int}, graph, p_pos, mat, Int, Set{Int}})
    precompile(Tuple{typeof(CrystalNets.group_cycle), Vector{Set{Int}}, types, graph})
    precompile(Tuple{typeof(CrystalNets.identify_metallic_type), Symbol, kinds, Int})
    precompile(Tuple{typeof(CrystalNets.find_sbus), cryst, kinds})
    precompile(Tuple{typeof(CrystalNets.find_sbus), crystclust, kinds})
    precompile(Tuple{typeof(CrystalNets._split_this_sbu!), Vector{Int}, graph, Int, types, Symbol, sbus})
    precompile(Tuple{typeof(CrystalNets.split_special_sbu!), graph, sbus, graph, types, Bool})
    precompile(Tuple{typeof(CrystalNets.split_O_vertices), cryst})
    precompile(Tuple{typeof(CrystalNets.split_O_vertices), crystclust})
    precompile(Tuple{typeof(CrystalNets.identify_clustering), cryst, structure, clustering})
    precompile(Tuple{typeof(CrystalNets.identify_clustering), crystclust, structure, clustering})
    precompile(Tuple{typeof(CrystalNets.order_atomtype), Symbol})
    precompile(Tuple{typeof(CrystalNets._collapse_clusters), cryst, clust})
    precompile(Tuple{typeof(CrystalNets._find_clusters), cryst, Bool, Bool})
    precompile(Tuple{typeof(CrystalNets._find_clusters), crystclust, Bool, Bool})
    precompile(Tuple{typeof(CrystalNets.find_clusters), cryst})
    precompile(Tuple{typeof(CrystalNets.find_clusters), crystclust})
    precompile(Tuple{typeof(CrystalNets.collapse_clusters), cryst})
    precompile(Tuple{typeof(CrystalNets.collapse_clusters), crystclust})
    precompile(Tuple{typeof(CrystalNets.update_new_edgs!), Dict{PeriodicEdge3D,Bool}, neighs, BitMatrix, BitVector, Vector{Int}})
    precompile(Tuple{typeof(sort), SVector{4,Int}})
    precompile(Tuple{typeof(CrystalNets.edges_of_convex_hull), neighs, Int, fmat, p_pos, BitVector, Set{SVector{4,Int}}, BitMatrix})
    precompile(Tuple{typeof(CrystalNets.edges_of_convex_hull), neighs, Int, fmat, p_pos, BitVector, Set{SVector{4,Int}}})
    precompile(Tuple{typeof(CrystalNets.regroup_toremove), cryst, Vector{Int}, Vector{Vector{Int}}, String})
    precompile(Tuple{typeof(CrystalNets.pem_to_pe), cryst})
    precompile(Tuple{typeof(CrystalNets.regroup_vmap), cryst, Vector{Int}, Vector{Int}, String})
    precompile(Tuple{typeof(CrystalNets.allnodes_to_singlenodes), cryst})

    # guessbonds.jl
    precompile(Tuple{typeof(CrystalNets.guess_bonds), p_pos, types, fmat, opts})

    # types.jl, 2/2
    precompile(Tuple{typeof(CrystalNets.separate_components), cryst})
    precompile(Tuple{typeof(CrystalNets.separate_components), crystclust})
    for D in 1:3
        precompile(Tuple{typeof(CrystalNets._collect_net!), Vector{CrystalNet{D}}, Dict{graph,Int}, Int, cryst, clustering})
        precompile(Tuple{typeof(CrystalNets.collect_nets), Vector{cryst}, Val{D}})
    end
    precompile(Tuple{Type{unets}, Vector{Tuple{Vector{Int},Vector{CrystalNet1D}}},
                                  Vector{Tuple{Vector{Int},Vector{CrystalNet2D}}},
                                  Vector{Tuple{Vector{Int},Vector{CrystalNet3D}}}})
    precompile(Tuple{Type{unets}})
    precompile(Tuple{typeof(CrystalNets._repeatgroups!), Expr, Int})
    precompile(Tuple{Type{unets}, cryst})
    precompile(Tuple{Type{unets}, crystclust})
    precompile(Tuple{typeof(CrystalNets.__warn_nonunique), Int})
    precompile(Tuple{typeof(CrystalNets.__throw_interpenetrating), Int})
    precompile(Tuple{typeof(CrystalNets.__throw_multiplenets), Int})
    precompile(Tuple{Type{CrystalNet}, cryst})
    precompile(Tuple{Type{CrystalNet}, crystclust})
    for D in 1:3
        for D2 in 1:3
            precompile(Tuple{Type{CrystalNet{D}}, cell, types, PeriodicGraph{D2}, opts})
        end
        precompile(Tuple{Type{CrystalNet{D}}, unets})
        for T in (Int64, Int128, BigInt)
            precompile(Tuple{typeof(CrystalNets._CrystalNet), cell, types, PeriodicGraph{D}, Matrix{Rational{T}}, opts})
        end
        precompile(Tuple{Type{CrystalNet{D}}, cell, types, PeriodicGraph{D}, opts})
        precompile(Tuple{Type{unets}, PeriodicGraph{D}, opts})
        precompile(Tuple{Type{unets}, Vector{PeriodicEdge{D}}, opts})
        precompile(Tuple{Type{unets}, PeriodicGraph{D}})
        precompile(Tuple{Type{unets}, Vector{PeriodicEdge{D}}})
        precompile(Tuple{Type{CrystalNet}, PeriodicGraph{D}, opts})
        precompile(Tuple{Type{CrystalNet}, Vector{PeriodicEdge{D}}, opts})
        precompile(Tuple{Type{CrystalNet}, PeriodicGraph{D}})
        precompile(Tuple{Type{CrystalNet}, Vector{PeriodicEdge{D}}})
    end
    precompile(Tuple{Type{CrystalNet}, unets})
    precompile(Tuple{Type{unets}, String, opts})
    for D in 1:3
        precompile(Tuple{Type{CrystalNet{D}}, PeriodicGraph{D}, opts})
        precompile(Tuple{Type{CrystalNet{D}}, Vector{PeriodicEdge{D}}, opts})
        precompile(Tuple{Type{CrystalNet{D}}, PeriodicGraph{D}})
        precompile(Tuple{Type{CrystalNet{D}}, Vector{PeriodicEdge{D}}})
        precompile(Tuple{Type{genome}, PeriodicGraph{D}, Nothing, Bool, String})
        precompile(Tuple{Type{genome}, PeriodicGraph{D}, String, Bool, String})
        precompile(Tuple{Type{genome}, PeriodicGraph{D}, Nothing, Bool})
        precompile(Tuple{Type{genome}, PeriodicGraph{D}, String, Bool})
        precompile(Tuple{Type{genome}, PeriodicGraph{D}, Nothing})
        precompile(Tuple{Type{genome}, PeriodicGraph{D}, String})
    end
    precompile(Tuple{Type{genome}, String})
    precompile(Tuple{Type{genome}})
    precompile(Tuple{typeof(==), genome, genome})
    precompile(Tuple{typeof(hash), genome, UInt})
    precompile(Tuple{typeof(show), Base.TTY, genome})
    precompile(Tuple{typeof(show), Base.IOStream, genome})
    precompile(Tuple{typeof(parse), Type{genome}, String})
    precompile(Tuple{typeof(parse), Type{genome}, SubString{String}})
    precompile(Tuple{Type{result}, SizedVector{8,TopologicalGenome}, MVector{8,Int8}, Vector{Int8}})
    precompile(Tuple{Type{result}})
    precompile(Tuple{Type{result}, Vector{Tuple{clustering,Union{clustering,genome}}}})
    precompile(Tuple{typeof(==), result, result})
    precompile(Tuple{typeof(hash), result, UInt})
    precompile(Tuple{typeof(get), CrystalNets.Returns, result, clustering})
    precompile(Tuple{typeof(get), result, clustering, Nothing})
    precompile(Tuple{typeof(getindex), result, clustering})
    precompile(Tuple{typeof(getindex), result, Symbol})
    precompile(Tuple{typeof(getindex), result, Int})
    precompile(Tuple{typeof(setindex!), result, genome, clustering})
    precompile(Tuple{typeof(setindex!), result, clustering, clustering})
    precompile(Tuple{Type{result}, String})
    precompile(Tuple{typeof(setindex!), result, Nothing, clustering})
    precompile(Tuple{typeof(show), Base.TTY, result})
    precompile(Tuple{typeof(show), Base.IOStream, result})
    precompile(Tuple{typeof(iterate), result, Int})
    precompile(Tuple{typeof(iterate), result})
    precompile(Tuple{typeof(eltype), Type{result}})
    precompile(Tuple{typeof(length), result})
    precompile(Tuple{typeof(parse), Type{result}, String})
    precompile(Tuple{typeof(parse), Type{result}, SubString{String}})
    precompile(Tuple{typeof(setindex!), result, Nothing, clustering})

    # input.jl
    precompile(Tuple{typeof(CrystalNets.parse_cif),String})
    precompile(Tuple{typeof(CrystalNets.parsestrip), Float64, String})
    precompile(Tuple{typeof(CrystalNets.parsestrip), Float32, String})
    precompile(Tuple{typeof(CrystalNets.parsestrip), BigFloat, String})
    precompile(Tuple{typeof(CrystalNets.parsestrip), String})
    precompile(Tuple{typeof(CrystalNets.popstring!), Dict{String, Union{Vector{String},String}}, String})
    precompile(Tuple{typeof(CrystalNets.popvecstring!), Dict{String, Union{Vector{String},String}}, String})
    precompile(Tuple{typeof(CrystalNets.CIF), Dict{String, Union{Vector{String},String}}})
    precompile(Tuple{typeof(CrystalNets.CIF), String})
    precompile(Tuple{typeof(CrystalNets.parse_arc), String})
    precompile(Tuple{typeof(CrystalNets.parse_arcs), String, Vector{String}})
    precompile(Tuple{typeof(CrystalNets.parse_arcs), String})
    precompile(Tuple{typeof(CrystalNets.parse_atom_name), String})
    precompile(Tuple{typeof(CrystalNets.parse_atom), String})
    precompile(Tuple{typeof(CrystalNets.chem_atoms), Chemfiles.Residue})
    precompile(Tuple{typeof(CrystalNets.attribute_residues), Vector{Chemfiles.Residue}, Int, Bool})
    precompile(Tuple{typeof(CrystalNets.check_collision), p_pos, fmat})
    precompile(Tuple{typeof(CrystalNets.fix_atom_on_a_bond!), graph, p_pos, fmat})
    precompile(Tuple{typeof(CrystalNets.least_plausible_neighbours), Vector{Float64}, Int})
    precompile(Tuple{typeof(CrystalNets.fix_valence!), PeriodicGraph3D, p_pos, types, Vector{Int}, Vector{Int}, fmat, Val{true}, opts})
    precompile(Tuple{typeof(CrystalNets.fix_valence!), PeriodicGraph3D, p_pos, types, Vector{Int}, Vector{Int}, fmat, Val{false}, opts})
    precompile(Tuple{typeof(CrystalNets.sanitize_removeatoms!), graph, p_pos, types, fmat, opts})
    precompile(Tuple{typeof(CrystalNets.remove_triangles!), graph, p_pos, types, fmat, edgs})
    precompile(Tuple{typeof(CrystalNets.remove_triangles!), graph, p_pos, types, fmat})
    precompile(Tuple{typeof(CrystalNets._detect_bent_bond), graph, p_pos, Int, Int, fmat})
    precompile(Tuple{typeof(CrystalNets.detect_bent_bond), graph, p_pos, Int, Int, fmat})
    precompile(Tuple{typeof(CrystalNets.sanity_checks!), graph, p_pos, types, fmat, opts})
    precompile(Tuple{typeof(CrystalNets._remove_homoatomic_bonds!), graph, types, Set{Symbol}, Bool})
    precompile(Tuple{typeof(CrystalNets._remove_homometallic_bonds!), graph, types, Vector{Int}})
    precompile(Tuple{typeof(CrystalNets.remove_homoatomic_bonds!), graph, types, Set{Symbol}, Bool})
    precompile(Tuple{typeof(CrystalNets.finalize_checks), cell, p_pos, types, Vector{Int}, bondlist, Bool, opts, String})
    precompile(Tuple{typeof(CrystalNets.parse_as_cif), opts, String})
    precompile(Tuple{typeof(CrystalNets.parse_as_chemfile), Chemfiles.Frame, opts, String})
    precompile(Tuple{typeof(CrystalNets.parse_chemfile), String, opts})
    # precompile_kwarg(Tuple{typeof(parse_chemfile), String}, (kwargs,))

    # precompile(Tuple{typeof(CrystalNets.), })

    # arithmetics.jl
    for D in 1:3
        for T in inttypes
            precompile(Tuple{typeof(CrystalNets.find_ratbasis), ppos{D,T}})
            precompile(Tuple{typeof(CrystalNets.normal_basis_rational), ppos{D,T}})
        end
    end

    # stability.jl
    precompile(Tuple{Type{collision}, PeriodicGraph{0}, Int, Vector{Int}})
    precompile(Tuple{typeof(==), collision, collision})
    precompile(Tuple{typeof(length), collision})
    precompile(Tuple{typeof(hash), collision, UInt})
    for D in 1:3
        for T in inttypes
            precompile(Tuple{typeof(CrystalNets.shrink_collisions), cnet{D,T}, Vector{UnitRange{Int}}})
        end
    end
    precompile(Tuple{typeof(CrystalNets._order_collision), SimpleGraph{Int}, Vector{Vector{Int}}})
    for D in 1:3
        precompile(Tuple{typeof(CrystalNets.order_collision), PeriodicGraph{D}, UnitRange{Int}})
    end
    precompile(Tuple{typeof(CrystalNets.collision_utils), Vector{collision}, Int, Nothing})
    precompile(Tuple{typeof(CrystalNets.collision_utils), Vector{collision}, Int, Vector{Int}})
    precompile(Tuple{typeof(CrystalNets.collision_utils), Vector{collision}, Int})
    precompile(Tuple{typeof(CrystalNets.collision_utils), Vector{collision}, Vector{Int}})
    for D in 1:3
        precompile(Tuple{typeof(CrystalNets.expand_collisions), Vector{collision}, PeriodicGraph{D}, Vector{Int}})
        precompile(Tuple{typeof(CrystalNets.unsorted_node), PeriodicGraph{D}, UnitRange{Int}})
        precompile(Tuple{Type{collision}, PeriodicGraph{D}, UnitRange{Int}, Nothing})
        precompile(Tuple{Type{collision}, PeriodicGraph{D}, UnitRange{Int}})
    end
    precompile(Tuple{Type{collision}, collision, Vector{Int}})
    for D in 1:3
        for T in inttypes
            precompile(Tuple{typeof(CrystalNets.collect_collisions), cnet{D,T}})
            precompile(Tuple{typeof(CrystalNets.collision_nodes), cnet{D,T}})
        end
    end

    # topology.jl
    for D in 1:3
        for T in reverse(inttypes)
            precompile(Tuple{typeof(CrystalNets.check_dimensionality), cnet{D,T}})
            precompile(Tuple{typeof(CrystalNets.check_valid_symmetry), cnet{D,T}, pos{D,T}, Vector{collision}, Nothing})
            precompile(Tuple{typeof(CrystalNets.check_valid_symmetry), cnet{D,T}, pos{D,T}, Vector{collision}})
            precompile(Tuple{typeof(CrystalNets.possible_translations), cnet{D,T}})
            precompile(Tuple{typeof(CrystalNets.find_all_valid_translations), cnet{D,T}, Vector{collision}})
            precompile(Tuple{typeof(CrystalNets.minimal_volume_matrix), NTuple{D, Vector{Tuple{Int,Int,pos{D,T}}}}})
            precompile(Tuple{typeof(CrystalNets.reduce_with_matrix), cnet{D,T}, smat{D,T,D*D}, Vector{collision}})
            precompile(Tuple{typeof(CrystalNets.minimize), cnet{D,T}, Vector{collision}})
            precompile(Tuple{typeof(CrystalNets.findfirstbasis), ppos{D,T}})
            precompile(Tuple{typeof(CrystalNets.findbasis), Vector{Tuple{Int,Int,pos{D,T}}}})
            precompile(Tuple{typeof(CrystalNets.candidate_key), cnet{D,T}, Int, smat{D,T,D*D}, Vector{Tuple{Int,Int,pos{D,T}}}})
        end
        precompile(Tuple{typeof(CrystalNets.partition_by_coordination_sequence), PeriodicGraph{D}, Nothing})
        precompile(Tuple{typeof(CrystalNets.partition_by_coordination_sequence), PeriodicGraph{D}, Vector{Vector{Int}}})
        precompile(Tuple{typeof(CrystalNets.partition_by_coordination_sequence), PeriodicGraph{D}})
        for T in reverse(inttypes)
            precompile(Tuple{typeof(CrystalNets.find_initial_candidates), cnet{D,T}, Vector{Vector{Int}}, Vector{Int}})
            precompile(Tuple{typeof(CrystalNets.find_candidates_onlyneighbors), cnet{D,T}, Vector{Vector{Int}}, Vector{Int}})
        end
    end
    for T in inttypes
        precompile(Tuple{typeof(CrystalNets.find_candidates_fallback), cnet{3,T}, Vector{Int}, Vector{Vector{Int}}, Vector{Int}})
        precompile(Tuple{typeof(CrystalNets.extract_through_symmetry), Dict{Int,Vector{smat{3,T,9}}}, Vector{Vector{Int}}, Vector{smat{3,Int,9}}})
        for D in 1:3
            precompile(Tuple{typeof(CrystalNets.find_candidates), cnet{D,T}, Vector{collision}})
            precompile(Tuple{typeof(CrystalNets.topological_key), cnet{D,T}, Vector{collision}})
            precompile(Tuple{typeof(CrystalNets.topological_key), cnet{D,T}})
        end
    end

    # query.jl
    precompile(Tuple{typeof(CrystalNets.recognize_topology), String, Dict{String,String}})
    precompile(Tuple{typeof(CrystalNets.recognize_topology), String})
    for T in inttypes
        precompile(Tuple{typeof(CrystalNets.topological_genome), cnet{0,T}, Vector{collision}})
    end
    for D in 1:3
        precompile(Tuple{typeof(CrystalNets.recognize_topology), PeriodicGraph{D}, Dict{String,String}})
        precompile(Tuple{typeof(CrystalNets.recognize_topology), PeriodicGraph{D}})
        for T in reverse(inttypes)
            precompile(Tuple{typeof(CrystalNets.topological_genome), cnet{D,T}, Vector{collision}})
            precompile(Tuple{typeof(CrystalNets.topological_genome), cnet{D,T}})
        end
    end
    precompile(Tuple{typeof(CrystalNets.topological_genome), unets})
    for D in 1:3
        precompile(Tuple{typeof(CrystalNets.topological_genome), PeriodicGraph{D}, opts})
        precompile(Tuple{typeof(CrystalNets.topological_genome), PeriodicGraph{D}})
    end
    precompile(Tuple{typeof(CrystalNets.topological_genome), String, opts})
    precompile(Tuple{typeof(CrystalNets.topological_genome), String})
    precompile(Tuple{typeof(CrystalNets._loop_group!), Expr, Symbol, Symbol, Expr})
    precompile(Tuple{typeof(CrystalNets.determine_topology), String, opts})
    precompile(Tuple{typeof(CrystalNets.determine_topology), String})
    precompile(Tuple{typeof(CrystalNets.guess_topology), String, opts})
    precompile(Tuple{typeof(CrystalNets.guess_topology), String})

    # The following are not precompiled because execution time always dominates compilation time
    # precompile(Tuple{typeof(CrystalNets.determine_topology_dataset), String, Bool, Bool, opts})
    # precompile(Tuple{typeof(CrystalNets.determine_topology_dataset), String, Bool, Bool})
    # precompile(Tuple{typeof(CrystalNets.determine_topology_dataset), String, Bool})
    # precompile(Tuple{typeof(CrystalNets.determine_topology_dataset), String})
    # precompile(Tuple{typeof(CrystalNets.guess_topology_dataset), String, Bool, Bool, opts})
    # precompile(Tuple{typeof(CrystalNets.guess_topology_dataset), String, Bool, Bool})
    # precompile(Tuple{typeof(CrystalNets.guess_topology_dataset), String, Bool})
    # precompile(Tuple{typeof(CrystalNets.guess_topology_dataset), String})

    # executable.jl
    precompile(Tuple{typeof(CrystalNets.parse_commandline), Vector{String}})
    precompile(Tuple{typeof(CrystalNets.main), Vector{String}})
    precompile(Tuple{typeof(CrystalNets.main), Vector{SubString{String}}})
    precompile(Tuple{typeof(CrystalNets.main), String})
    precompile(Tuple{typeof(CrystalNets.julia_main)})
end

# _precompile_()

#=
if ccall(:jl_generating_output, Cint, ()) == 1 # precompilation
    tw = DOWARN[]; toggle_warning(false)
    te = DOEXPORT[]; toggle_export(false)

    topological_genome(CrystalNet(PeriodicGraph("1 1 1 1")))

    topological_genome(parse(TopologicalGenome, "hcb").genome)

    cifs = joinpath(dirname(dirname(pathof(CrystalNets))), "test", "cif")
    mof5 = joinpath(cifs, "MOF-5.cif")
    determine_topology(mof5, CrystalNets.Options(structure=StructureType.MOF, throw_error=true))

    im19 = joinpath(cifs, "IM-19.cif")
    determine_topology(im19; structure=StructureType.MOF,
                             clusterings=[Clustering.PEM,
                                          Clustering.PE,
                                          Clustering.AllNodes,
                                          Clustering.SingleNodes,
                                          Clustering.Standard],
                             throw_error=true)

    toggle_warning(tw)
    toggle_export(te)
end
=#
