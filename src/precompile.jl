using CrystalNets, PeriodicGraphs, ArgParse, LinearAlgebra, SparseArrays,
      StaticArrays, Logging, Tokenize
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

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    # Tokenize
    Base.precompile(Tuple{typeof(Base._collect),UnitRange{Int},Tokenize.Lexers.Lexer{IOBuffer, Tokenize.Tokens.Token},Base.HasEltype,Base.SizeUnknown})

    # Graphs
    Base.precompile(Tuple{typeof(Graphs.floyd_warshall_shortest_paths),Graphs.SimpleGraph{Int},Graphs.DefaultDistance})
    Base.precompile(Tuple{Type{Graphs.SimpleGraph},Vector{Graphs.SimpleGraphs.SimpleEdge{Int}}})

    # PeriodicGraphs
    Base.precompile(Tuple{Type{Dict{PeriodicEdge3D, Nothing}}})
    Base.precompile(Tuple{Type{Dict{PeriodicVertex3D, Nothing}}})
    Base.precompile(Tuple{typeof(Base._unique!),typeof(identity),Vector{PeriodicEdge3D},Set{PeriodicEdge3D},Int,Int})
    Base.precompile(Tuple{typeof(Base.ht_keyindex),Dict{PeriodicVertex3D, Nothing},PeriodicVertex3D})
    Base.precompile(Tuple{typeof(copyto!),Vector{PeriodicEdge3D},PeriodicGraphs.PeriodicEdgeIter{3}})
    Base.precompile(Tuple{typeof(deleteat!),Vector{PeriodicVertex3D},BitVector})
    Base.precompile(Tuple{typeof(has_edge),PeriodicGraph3D,PeriodicEdge3D})
    Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph3D,Vector{SVector{3,Int}}})
    Base.precompile(Tuple{typeof(searchsortedfirst),Vector{PeriodicVertex3D},PeriodicVertex3D,Int,Int,Base.Order.ForwardOrdering})
    Base.precompile(Tuple{typeof(setindex!),Dict{PeriodicEdge3D, Nothing},Nothing,PeriodicEdge3D})
    Base.precompile(Tuple{typeof(sort!),Vector{PeriodicEdge3D},Int,Int,Base.Sort.MergeSortAlg,Base.Order.ForwardOrdering,Vector{PeriodicEdge3D}})
    Base.precompile(Tuple{typeof(sort!),Vector{PeriodicVertex3D},Int,Int,Base.Sort.MergeSortAlg,Base.Order.ForwardOrdering,Vector{PeriodicVertex3D}})
    Base.precompile(Tuple{typeof(union!),Set{PeriodicVertex3D},Vector{PeriodicVertex3D}})
    Base.precompile(Tuple{typeof(Base.setindex!),Dict{PeriodicEdge3D, Nothing},Nothing,PeriodicEdge3D})
    Base.precompile(Tuple{typeof(Base.union!),Set{PeriodicVertex3D},Vector{PeriodicVertex3D}})

    # SparseArrays
    Base.precompile(Tuple{typeof(*),SparseArrays.SparseMatrixCSC{Int, Int},Matrix{Rational{BigInt}}})
    Base.precompile(Tuple{typeof(Base.copyto_unaliased!),IndexCartesian,SubArray{Int, 1, SparseArrays.SparseMatrixCSC{Int, Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, false},IndexLinear,Vector{Int}})
    Base.precompile(Tuple{typeof(Base.mightalias),SubArray{Int, 1, SparseArrays.SparseMatrixCSC{Int, Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, false},Vector{Int}})
    Base.precompile(Tuple{typeof(LinearAlgebra.mul!),Matrix{BigInt},SparseArrays.SparseMatrixCSC{Int, Int},Matrix{BigInt},Bool,Bool})
    Base.precompile(Tuple{typeof(SparseArrays.dimlub),Vector{Int}})
    Base.precompile(Tuple{typeof(SparseArrays.findnz),SparseArrays.SparseMatrixCSC{Int, Int}})
    Base.precompile(Tuple{typeof(SparseArrays.sparse_check_length),String,Vector{Int},Int,Type})
    Base.precompile(Tuple{typeof(SparseArrays.spzeros),Type{Int},Type{Int},Int,Int})
    Base.precompile(Tuple{typeof(getindex),SparseArrays.SparseMatrixCSC{Int, Int},UnitRange{Int},UnitRange{Int}})

    # Logging
    Base.precompile(Tuple{typeof(Base.CoreLogging.handle_message),Logging.ConsoleLogger,Any,Any,Any,Any,Any,Any,Any})
    Base.precompile(Tuple{typeof(Base.CoreLogging.shouldlog),Logging.ConsoleLogger,Base.CoreLogging.LogLevel,Module,Symbol,Symbol})
    Base.precompile(Tuple{typeof(Logging.default_metafmt),Base.CoreLogging.LogLevel,Any,Any,Any,Any,Any})
    Base.precompile(Tuple{typeof(Logging.termlength),SubString{String}})
    let fbody = try __lookup_kwbody__(which(Base.CoreLogging.handle_message, (Logging.ConsoleLogger,Any,Any,Any,Any,Any,Any,Any,))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Any,typeof(Base.CoreLogging.handle_message),Logging.ConsoleLogger,Any,Any,Any,Any,Any,Any,Any,))
        end
    end

    # ArgParse
    Base.precompile(Tuple{Core.kwftype(typeof(ArgParse.Type)),Any,Type{ArgParse.ArgParseSettings}})
    Base.precompile(Tuple{Core.kwftype(typeof(ArgParse.add_arg_field!)),Any,typeof(ArgParse.add_arg_field!),ArgParse.ArgParseSettings,Vector{T} where T<:AbstractString})
    Base.precompile(Tuple{typeof(ArgParse.parse1_optarg!),ArgParse.ParserState,ArgParse.ArgParseSettings,ArgParse.ArgParseField,Any,AbstractString})
    Base.precompile(Tuple{typeof(ArgParse.preparse!),Channel,ArgParse.ParserState,ArgParse.ArgParseSettings})
    Base.precompile(Tuple{typeof(ArgParse.print_group),IO,Vector{T} where T,AbstractString,Int,Int,AbstractString,AbstractString,AbstractString})
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
    Base.precompile(Tuple{Type{Matrix{Float64}},LinearAlgebra.UniformScaling{Bool},Tuple{Int, Int}})
    Base.precompile(Tuple{typeof(LinearAlgebra.norm),Vector{Float64},Int})
    Base.precompile(Tuple{typeof(eltype),LinearAlgebra.Adjoint{Rational{Int}, Matrix{Rational{Int}}}})
    Base.precompile(Tuple{typeof(isone),Matrix{Int32}})
    Base.precompile(Tuple{typeof(hcat),Vector{Rational{Int}},LinearAlgebra.Adjoint{Rational{Int}, Matrix{Rational{Int}}}})
    let fbody = try __lookup_kwbody__(which(LinearAlgebra.rank, (Matrix{Int},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Float64,Float64,typeof(LinearAlgebra.rank),Matrix{Int},))
        end
    end

    # Base
    Base.precompile(Tuple{Core.kwftype(typeof(Base.with_output_color)),NamedTuple{(:bold,), Tuple{Bool}},typeof(Base.with_output_color),Function,Symbol,IOContext{Base.TTY},String,Vararg{Any, N} where N})
    Base.precompile(Tuple{Type{Base.IteratorSize},Base.Iterators.ProductIterator{Tuple{UnitRange{Int}, UnitRange{Int}}}})
    Base.precompile(Tuple{Type{Dict{Symbol, Any}},Base.Iterators.Pairs{Symbol, Any, Tuple{Symbol, Symbol, Symbol}, NamedTuple{(:help, :metavar, :required), Tuple{String, String, Bool}}}})
    Base.precompile(Tuple{Type{Dict{Symbol, Any}},Base.Iterators.Pairs{Symbol, Any, Tuple{Symbol, Symbol}, NamedTuple{(:help, :action), Tuple{String, Symbol}}}})
    Base.precompile(Tuple{Type{Dict{Symbol, Any}},Base.Iterators.Pairs{Symbol, String, Tuple{Symbol, Symbol}, NamedTuple{(:help, :metavar), Tuple{String, String}}}})
    Base.precompile(Tuple{Type{Tuple{Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}, Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}, Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}}},Vector{Tuple{Vector{Any}, Vector{Any}}}})
    for T in (Int8, Int16, Int32, Int64, Int128)
       Base.precompile(Tuple{Type{SubArray},IndexLinear,Matrix{Rational{T}},Tuple{Base.Slice{Base.OneTo{Int}}, Int},Tuple{Bool}})
       Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, Type{Rational{T}}, Tuple{Matrix{Rational{Int}}}}})
    end
    Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(big), Tuple{Matrix{Rational{Int}}}}})
    Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(denominator), Tuple{Matrix{Rational{Int}}}}})
    Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(numerator), Tuple{Matrix{Rational{Int}}}}})
    Base.precompile(Tuple{typeof(Base.deepcopy_internal),NTuple{9, BigFloat},IdDict{Any, Any}})
    Base.precompile(Tuple{typeof(Base.display_error),Base.TTY,ErrorException,Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}})
    Base.precompile(Tuple{typeof(Base.grow_to!),Dict{Symbol, Any},Tuple{Pair{Symbol, String}, Pair{Symbol, Bool}, Pair{Symbol, String}},Int})
    Base.precompile(Tuple{typeof(Base.grow_to!),Dict{Symbol, Any},Tuple{Pair{Symbol, String}, Pair{Symbol, Symbol}},Int})
    Base.precompile(Tuple{typeof(Base.print_to_string),Int128,Vararg{Any, N} where N})
    Base.precompile(Tuple{typeof(Base.reducedim_init),Function,typeof(min),Matrix{Float64},Int})
    Base.precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Int},Expr,Int})
    Base.precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{Any}},Tuple{},Int})
    Base.precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{Int}},Tuple{},Int})
    Base.precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{}},Tuple{Int},Int})
    Base.precompile(Tuple{typeof(Base.typed_hvcat),Type{BigFloat},Tuple{Int, Int, Int},BigFloat,Vararg{Number, N} where N})
    Base.precompile(Tuple{typeof(Base.typed_hvcat),Type{Float64},Tuple{Int, Int, Int},Int,Vararg{Number, N} where N})
    Base.precompile(Tuple{typeof(copy),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Tuple{Base.OneTo{Int}, Base.OneTo{Int}}, Type{Int}, Tuple{Matrix{Rational{BigInt}}}}})
    Base.precompile(Tuple{typeof(maximum),Matrix{Int}})
    Base.precompile(Tuple{typeof(minimum),Matrix{Int}})
    Base.precompile(Tuple{typeof(push!),Vector{String},String,String,String})
    Base.precompile(Tuple{typeof(setindex!),Dict{String, Tuple{Exception, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}},Tuple{ArgumentError, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}},String})
    Base.precompile(Tuple{typeof(string),Int128,String,Vararg{Any, N} where N})
    Base.precompile(Tuple{typeof(vcat),Vector{Expr},Vector{Expr}})
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
    Base.precompile(Tuple{Type{Vector{SVector{3, Float64}}},Vector{SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true}}})
    Base.precompile(Tuple{Type{Vector{SVector{3, Int}}},Vector{Any}})
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{0},StaticArrays.StaticArrayStyle{2}})
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},StaticArrays.StaticArrayStyle{1}})
    Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(+),Tuple{SVector{3, Int}, SVector{3, Int}}})
    Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(-),Tuple{SVector{3, Rational{Int}}, SVector{3, Int}}})
    Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(floor),Tuple{Base.RefValue{Type{Int}}, SVector{3, Rational{Int}}}})
    Base.precompile(Tuple{Type{Dict{Int, SVector{3, Int}}}})
    Base.precompile(Tuple{Type{Dict{SVector{3, Int}, Nothing}}})
    Base.precompile(Tuple{Type{SMatrix{3, 3, BigFloat, 9}},NTuple{9, Float64}})
    Base.precompile(Tuple{typeof(==),Vector{Tuple{Int, Int, SVector{3, Rational{Int}}}},Vector{Tuple{Int, Int, SVector{3, Rational{Int}}}}})
    Base.precompile(Tuple{typeof(>),SVector{3, Int},SizedVector{3, Int, 1}})
    Base.precompile(Tuple{typeof(>),SizedVector{3, Int, 1},SVector{3, Int}})
    Base.precompile(Tuple{typeof(Base.Broadcast.broadcasted),Function,SVector{3, Rational{Int}},SVector{3, Int}})
    Base.precompile(Tuple{typeof(StaticArrays.arithmetic_closure),Type{BigFloat}})
    Base.precompile(Tuple{typeof(LinearAlgebra.generic_norm2),SubArray{BigFloat, 1, SMatrix{3, 3, BigFloat, 9}, Tuple{Base.Slice{SOneTo{3}}, Int}, true}})
    Base.precompile(Tuple{typeof(StaticArrays._axes),Size{(3, 3)}})
    Base.precompile(Tuple{typeof(StaticArrays._axes),Size{(3,)}})

    Base.precompile(Tuple{Type{Size},Type{SubArray{Int, 1, Matrix{Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, true}}})
    Base.precompile(Tuple{typeof(<),SVector{3, Int},SizedVector{3, Int, 1}})
    for T in (Int8, Int16, Int32, Int64, Int128)
        Base.precompile(Tuple{Type{Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}},Vector{Any}})
        Base.precompile(Tuple{Type{SVector{3, Int}},Tuple{Rational{T}, Rational{T}, Rational{T}}})
        Base.precompile(Tuple{Type{Dict{Int, Vector{SMatrix{3, 3, Rational{T}, 9}}}}})
        Base.precompile(Tuple{Type{Dict{SVector{9, Rational{T}}, Nothing}}})
        Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(isempty),Tuple{Tuple{Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}}}})
        Base.precompile(Tuple{typeof(<),SVector{3, Rational{T}},SVector{3, Rational{Int}}})
        Base.precompile(Tuple{typeof(<),SVector{3, Rational{Int}},SVector{3, Rational{T}}})
        Base.precompile(Tuple{typeof(==),SVector{3, Rational{Int}},SVector{3, Rational{T}}})
        Base.precompile(Tuple{typeof(==),SMatrix{3, 3, T, 9},LinearAlgebra.UniformScaling{Bool}})
        Base.precompile(Tuple{typeof(LinearAlgebra.generic_matmatmul!),Matrix{Rational{widen(T)}},Char,Char,SizedMatrix{3, 3, Rational{widen(T)}, 2},Matrix{Rational{T}},LinearAlgebra.MulAddMul{true, true, Bool, Bool}})
        Base.precompile(Tuple{typeof(LinearAlgebra.dot),SubArray{BigFloat, 1, SMatrix{3, 3, BigFloat, 9}, Tuple{Base.Slice{SOneTo{3}}, T}, true},SubArray{BigFloat, 1, SMatrix{3, 3, BigFloat, 9}, Tuple{Base.Slice{SOneTo{3}}, T}, true}})
        Base.precompile(Tuple{typeof(Base.Broadcast.broadcasted),Base.Broadcast.Style{Tuple},Function,Tuple{Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}}})
        Base.precompile(Tuple{typeof(Base.Broadcast.materialize!),Base.Broadcast.DefaultArrayStyle{1},SubArray{Rational{CrystalNets.soft_widen(T)}, 1, Matrix{Rational{CrystalNets.soft_widen(T)}}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true},Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(-), Tuple{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(+), Tuple{SVector{3, Rational{T}}, SVector{3, Int}}}, SVector{3, Rational{T}}}}})
        Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(-), Tuple{SVector{3, Rational{T}}, SVector{3, Int}}}})
        Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(floor), Tuple{Base.RefValue{Type{Int}}, SVector{3, Rational{T}}}}})
        Base.precompile(Tuple{typeof(append!),Vector{SMatrix{3, 3, Rational{T}, 9}},Vector{SMatrix{3, 3, Rational{T}, 9}}})
        Base.precompile(Tuple{typeof(append!),Vector{SMatrix{3, 3, Rational{T}, 9}},Vector{SMatrix{3, 3, Rational{widen(T)}, 9}}})
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
        Base.precompile(Tuple{typeof(sort!),Vector{Int},Base.Sort.QuickSortAlg,Base.Order.Perm{Base.Order.ForwardOrdering, Vector{SVector{3, Rational{T}}}}})
     end
     Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(//), Tuple{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(round), Tuple{Base.RefValue{Type{Int128}}, Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(*), Tuple{Int128, SVector{3, Float64}}}}}, Int128}}})
     Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(//), Tuple{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(round), Tuple{Base.RefValue{Type{BigInt}}, Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(*), Tuple{BigInt, SVector{3, Float64}}}}}, BigInt}}})
     Base.precompile(Tuple{typeof(Base._foldl_impl),Base.BottomRF{typeof(hcat)},Base.ReshapedArray{Int, 2, SVector{3, Int}, Tuple{}},Set{SVector{3, Int}}})
     Base.precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{StaticArrays.Dynamic}},Tuple{Int},Int})
     Base.precompile(Tuple{typeof(StaticArrays._axes),Size{(9,)}})
     Base.precompile(Tuple{typeof(setindex!),Dict{Int, SVector{3, Int}},SVector{3, Int},Int})
     Base.precompile(Tuple{typeof(setindex!),Dict{SVector{3, Int}, Nothing},Nothing,SVector{3, Int}})
     Base.precompile(Tuple{typeof(which(StaticArrays._getindex,(AbstractArray,Tuple{Vararg{Size, N} where N},Any,)).generator.gen),Any,Any,Any,Any})
     Base.precompile(Tuple{typeof(which(StaticArrays._map,(Any,Vararg{AbstractArray, N} where N,)).generator.gen),Any,Any,Any})
     Base.precompile(Tuple{typeof(which(StaticArrays.combine_sizes,(Tuple{Vararg{Size, N} where N},)).generator.gen),Any,Any})
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
     Base.precompile(Tuple{typeof(sort!),Vector{Int},Base.Sort.QuickSortAlg,Base.Order.Perm{Base.Order.ForwardOrdering, Vector{SVector{3, Float64}}}})


    # Chemfiles
    Base.precompile(Tuple{Type{Chemfiles.Atom},Chemfiles.Frame,Int})
    Base.precompile(Tuple{Type{Chemfiles.Atom},String})
    Base.precompile(Tuple{Type{Chemfiles.Frame}})
    Base.precompile(Tuple{Type{Chemfiles.Topology},Chemfiles.Frame})
    Base.precompile(Tuple{Type{Chemfiles.UnitCell},BigFloat,BigFloat,BigFloat,BigFloat,BigFloat,BigFloat})
    Base.precompile(Tuple{Type{Chemfiles.UnitCell},Chemfiles.Frame})
    Base.precompile(Tuple{typeof(Chemfiles.__free),Chemfiles.CxxPointer{Chemfiles.lib.CHFL_ATOM}})
    Base.precompile(Tuple{typeof(Chemfiles.__free),Chemfiles.CxxPointer{Chemfiles.lib.CHFL_CELL}})
    Base.precompile(Tuple{typeof(Chemfiles.__free),Chemfiles.CxxPointer{Chemfiles.lib.CHFL_TOPOLOGY}})
    Base.precompile(Tuple{typeof(Chemfiles.__strip_null),String})
    Base.precompile(Tuple{typeof(Chemfiles.add_atom!),Chemfiles.Frame,Chemfiles.Atom,Vector{Float64},Vector{Float64}})
    Base.precompile(Tuple{typeof(Chemfiles.bonds),Chemfiles.Topology})
    Base.precompile(Tuple{typeof(Chemfiles.matrix),Chemfiles.UnitCell})
    Base.precompile(Tuple{typeof(Chemfiles.positions),Chemfiles.Frame})
    Base.precompile(Tuple{typeof(Chemfiles.set_type!),Chemfiles.Atom,String})


    # CrystalNets
    Base.precompile(Tuple{Core.kwftype(typeof(CrystalNets.recognize_topologies)),NamedTuple{(:ignore_atoms,), Tuple{Tuple{Symbol}}},typeof(CrystalNets.recognize_topologies),String})
    # Base.precompile(Tuple{CrystalNets.var"#730#threadsfor_fun#121"{Tuple{Symbol}, String, Dict{String, String}, Dict{String, Tuple{Exception, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}}, Dict{String, String}, Vector{String}}})
    for T in (Bool, Int8, Int16, Int32, Int64, Int128, BigInt)
       # Base.precompile(Tuple{CrystalNets.var"#670#threadsfor_fun#114"{CrystalNet{Rational{T}}, Vector{Int}, Base.Threads.SpinLock, Vector{Pair{Int, Tuple{Matrix{Rational{CrystalNets.soft_widen(T)}}, Vector{Int}}}}, Int, DataType, Vector{Int}}})
       # Base.precompile(Tuple{CrystalNets.var"#685#threadsfor_fun#116"{Rational{T}, Dict{Int, Vector{SMatrix{3, 3, Rational{CrystalNets.soft_widen(T)}, 9}}}, Base.Threads.SpinLock, DataType, Vector{Pair{Int, Tuple{Matrix{Rational{CrystalNets.soft_widen(T)}}, Vector{Int}}}}}})

       Base.precompile(Tuple{Type{CrystalNet{Rational{T}}},CrystalNets.Cell,Vector{Symbol},PeriodicGraph3D,Matrix{Rational{T}}})
       Base.precompile(Tuple{typeof(CrystalNets.isrank3),Matrix{Rational{T}}})
       Base.precompile(Tuple{typeof(CrystalNets.minimize),CrystalNet{Rational{T}}})
       Base.precompile(Tuple{typeof(CrystalNets.topological_key),CrystalNet{Rational{T}}})
       Base.precompile(Tuple{typeof(topological_genome),CrystalNet{Rational{T}}})

    end
    Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(CrystalNets.parsestrip), Tuple{Vector{String}}}})
    Base.precompile(Tuple{typeof(CrystalNets.invalid_input_error),String,ErrorException,Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}})
    Base.precompile(Tuple{typeof(CrystalNets.julia_main)})
    Base.precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{2147483647, Int32}},Int,Int,Function})
    Base.precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{2147483629, Int32}},Int,Int,Function})
    Base.precompile(Tuple{typeof(CrystalNets.parse_chemfile),String})
    Base.precompile(Tuple{typeof(CrystalNets.parse_chemfile),String,Bool})
    Base.precompile(Tuple{Core.kwftype(typeof(CrystalNets.parse_chemfile)),NamedTuple{(:ignore_atoms,), Tuple{Tuple{Symbol}}},typeof(CrystalNets.recognize_topologies),String})
    Base.precompile(Tuple{Core.kwftype(typeof(CrystalNets.parse_chemfile)),NamedTuple{(:ignore_atoms,), Tuple{Tuple{Symbol}}},typeof(CrystalNets.recognize_topologies),String,Bool})
    Base.precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(CrystalNets.parsestrip),Tuple{Vector{String}}})
    Base.precompile(Tuple{typeof(-),CrystalNets.Modulos.Modulo{2147483647, Int32},CrystalNets.Modulos.Modulo{2147483647, Int32}})
    Base.precompile(Tuple{typeof(/),Int,CrystalNets.Modulos.Modulo{2147483647, Int32}})
    Base.precompile(Tuple{typeof(==),Matrix{CrystalNets.Modulos.Modulo{2147483647, Int32}},Matrix{Int}})
    Base.precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(CrystalNets.parsestrip), Tuple{Vector{String}}}})
    Base.precompile(Tuple{typeof(Base._unsafe_copyto!),Matrix{CrystalNets.Modulos.Modulo{2147483647, Int32}},Int,Matrix{Int},Int,Int})
    Base.precompile(Tuple{typeof(Base._unsafe_getindex),IndexLinear,Matrix{CrystalNets.Modulos.Modulo{2147483647, Int32}},Int,Base.Slice{Base.OneTo{Int}}})
    Base.precompile(Tuple{typeof(Base.copyto_unaliased!),IndexLinear,SubArray{CrystalNets.Modulos.Modulo{2147483647, Int32}, 1, Matrix{CrystalNets.Modulos.Modulo{2147483647, Int32}}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, true},IndexLinear,Vector{CrystalNets.Modulos.Modulo{2147483647, Int32}}})
    Base.precompile(Tuple{typeof(Base.deepcopy_internal),Vector{CrystalNets.EquivalentPosition},IdDict{Any, Any}})
    Base.precompile(Tuple{typeof(CrystalNets.isrank3),Matrix{Rational{Int32}}})
    Base.precompile(Tuple{typeof(LinearAlgebra.mul!),Matrix{CrystalNets.Modulos.Modulo{2147483647, Int32}},SparseArrays.SparseMatrixCSC{Int, Int},Matrix{CrystalNets.Modulos.Modulo{2147483647, Int32}},Bool,Bool})
    Base.precompile(Tuple{typeof(SparseArrays._setindex_scalar!),SparseArrays.SparseMatrixCSC{CrystalNets.Modulos.Modulo{2147483647, Int32}, Int},Int,Int,Int})
    Base.precompile(Tuple{typeof(SparseArrays.sparse!),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{2147483647, Int32}},Int,Int,typeof(+),Vector{Int},Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{2147483647, Int32}},Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{2147483647, Int32}}})
    Base.precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{2147483647, Int32}},Int,Int,Function})
    Base.precompile(Tuple{typeof(SparseArrays.sparse_check_length),String,Vector{CrystalNets.Modulos.Modulo{2147483647, Int32}},Int,Type})
    Base.precompile(Tuple{typeof(getindex),SparseArrays.SparseMatrixCSC{CrystalNets.Modulos.Modulo{2147483647, Int32}, Int},UnitRange{Int},UnitRange{Int}})
    Base.precompile(Tuple{typeof(parse_chemfile),String})
    Base.precompile(Tuple{typeof(recognize_topology),PeriodicGraph3D})
    for T in (Nothing, CrystalNets.Clusters)
        Base.precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}})
        Base.precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, Nothing})
        Base.precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, String})
        Base.precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, Nothing, String})
        Base.precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, String, String})
    end
end

#_precompile_()
