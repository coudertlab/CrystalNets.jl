using CrystalNets, PeriodicGraphs, ArgParse, LinearAlgebra, SparseArrays,
      StaticArrays, Logging, Tokenize, BigRationals
import Chemfiles

macro enforce(expr) # strong @assert
    msg = string(expr)
    return :($(esc(expr)) ? $(nothing) : throw(AssertionError($msg)))
end


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

# TODO: check whether it works
function precompile_kwarg(@nospecialize(tt), @nospecialize(kwargtypes))
    let fbody = try __lookup_kwbody__(which(tt)) catch; missing end
        @enforce !ismissing(fbody)
        ttu = Base.unwrap_unionall(tt)
        newtt = Tuple{Core.Typeof(fbody), kwargtypes..., ttu.parameters...}
        #=@enforce=# precompile(Base.rewrap_unionall(newtt, tt))
    end
end

function _precompile_dependencies()
    # Tokenize
    @enforce precompile(Tuple{typeof(Base._collect),UnitRange{Int},Tokenize.Lexers.Lexer{IOBuffer, Tokenize.Tokens.Token},Base.HasEltype,Base.SizeUnknown})

    # Graphs
    @enforce precompile(Tuple{typeof(Graphs.floyd_warshall_shortest_paths),Graphs.SimpleGraph{Int},Graphs.DefaultDistance})
    @enforce precompile(Tuple{Type{Graphs.SimpleGraph},Vector{Graphs.SimpleGraphs.SimpleEdge{Int}}})

    # PeriodicGraphs
    @enforce precompile(Tuple{Type{Dict{PeriodicEdge3D, Nothing}}})
    @enforce precompile(Tuple{Type{Dict{PeriodicVertex3D, Nothing}}})
    @enforce precompile(Tuple{typeof(Base._unique!),typeof(identity),Vector{PeriodicEdge3D},Set{PeriodicEdge3D},Int,Int})
    @enforce precompile(Tuple{typeof(Base.ht_keyindex),Dict{PeriodicVertex3D, Nothing},PeriodicVertex3D})
    @enforce precompile(Tuple{typeof(copyto!),Vector{PeriodicEdge3D},PeriodicGraphs.PeriodicEdgeIter{3}})
    @enforce precompile(Tuple{typeof(deleteat!),Vector{PeriodicVertex3D},BitVector})
    @enforce precompile(Tuple{typeof(has_edge),PeriodicGraph3D,PeriodicEdge3D})
    @enforce precompile(Tuple{typeof(offset_representatives!),PeriodicGraph3D,Vector{SVector{3,Int}}})
    @enforce precompile(Tuple{typeof(searchsortedfirst),Vector{PeriodicVertex3D},PeriodicVertex3D,Int,Int,Base.Order.ForwardOrdering})
    @enforce precompile(Tuple{typeof(setindex!),Dict{PeriodicEdge3D, Nothing},Nothing,PeriodicEdge3D})
    @enforce precompile(Tuple{typeof(sort!),Vector{PeriodicEdge3D},Int,Int,Base.Sort.MergeSortAlg,Base.Order.ForwardOrdering,Vector{PeriodicEdge3D}})
    @enforce precompile(Tuple{typeof(sort!),Vector{PeriodicVertex3D},Int,Int,Base.Sort.MergeSortAlg,Base.Order.ForwardOrdering,Vector{PeriodicVertex3D}})
    @enforce precompile(Tuple{typeof(union!),Set{PeriodicVertex3D},Vector{PeriodicVertex3D}})
    @enforce precompile(Tuple{typeof(Base.setindex!),Dict{PeriodicEdge3D, Nothing},Nothing,PeriodicEdge3D})
    @enforce precompile(Tuple{typeof(Base.union!),Set{PeriodicVertex3D},Vector{PeriodicVertex3D}})

    # SparseArrays
    @enforce precompile(Tuple{typeof(*),SparseArrays.SparseMatrixCSC{Int, Int},Matrix{Rational{BigInt}}})
    @enforce precompile(Tuple{typeof(Base.copyto_unaliased!),IndexCartesian,SubArray{Int, 1, SparseArrays.SparseMatrixCSC{Int, Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, false},IndexLinear,Vector{Int}})
    @enforce precompile(Tuple{typeof(Base.mightalias),SubArray{Int, 1, SparseArrays.SparseMatrixCSC{Int, Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, false},Vector{Int}})
    @enforce precompile(Tuple{typeof(LinearAlgebra.mul!),Matrix{BigInt},SparseArrays.SparseMatrixCSC{Int, Int},Matrix{BigInt},Bool,Bool})
    @enforce precompile(Tuple{typeof(SparseArrays.dimlub),Vector{Int}})
    @enforce precompile(Tuple{typeof(SparseArrays.findnz),SparseArrays.SparseMatrixCSC{Int, Int}})
    @enforce precompile(Tuple{typeof(SparseArrays.sparse_check_length),String,Vector{Int},Int,Type})
    @enforce precompile(Tuple{typeof(SparseArrays.spzeros),Type{Int},Type{Int},Int,Int})
    @enforce precompile(Tuple{typeof(getindex),SparseArrays.SparseMatrixCSC{Int, Int},UnitRange{Int},UnitRange{Int}})

    # Logging
    @enforce precompile(Tuple{typeof(Base.CoreLogging.handle_message),Logging.ConsoleLogger,Any,Any,Any,Any,Any,Any,Any})
    @enforce precompile(Tuple{typeof(Base.CoreLogging.shouldlog),Logging.ConsoleLogger,Base.CoreLogging.LogLevel,Module,Symbol,Symbol})
    @enforce precompile(Tuple{typeof(Logging.default_metafmt),Base.CoreLogging.LogLevel,Any,Any,Any,Any,Any})
    @enforce precompile(Tuple{typeof(Logging.termlength),SubString{String}})
    let fbody = try __lookup_kwbody__(which(Base.CoreLogging.handle_message, (Logging.ConsoleLogger,Any,Any,Any,Any,Any,Any,Any,))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Any,typeof(Base.CoreLogging.handle_message),Logging.ConsoleLogger,Any,Any,Any,Any,Any,Any,Any,))
        end
    end

    # ArgParse
    @enforce precompile(Tuple{Core.kwftype(typeof(ArgParse.Type)),Any,Type{ArgParse.ArgParseSettings}})
    @enforce precompile(Tuple{Core.kwftype(typeof(ArgParse.add_arg_field!)),Any,typeof(ArgParse.add_arg_field!),ArgParse.ArgParseSettings,Vector{T} where T<:AbstractString})
    @enforce precompile(Tuple{typeof(ArgParse.parse1_optarg!),ArgParse.ParserState,ArgParse.ArgParseSettings,ArgParse.ArgParseField,Any,AbstractString})
    @enforce precompile(Tuple{typeof(ArgParse.preparse!),Channel,ArgParse.ParserState,ArgParse.ArgParseSettings})
    @enforce precompile(Tuple{typeof(ArgParse.print_group),IO,Vector{T} where T,AbstractString,Int,Int,AbstractString,AbstractString,AbstractString})
    let fbody = try __lookup_kwbody__(which(ArgParse.show_help, (ArgParse.ArgParseSettings,))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Any,typeof(ArgParse.show_help),ArgParse.ArgParseSettings,))
        end
    end
    let fbody = try __lookup_kwbody__(which(any, (Function,Vector{ArgParse.ArgParseGroup},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Function,typeof(any),Function,Vector{ArgParse.ArgParseGroup},))
        end
    end

    # LinearAlgebra
    @enforce precompile(Tuple{Type{Matrix{Float64}},LinearAlgebra.UniformScaling{Bool},Tuple{Int, Int}})
    @enforce precompile(Tuple{typeof(LinearAlgebra.norm),Vector{Float64},Int})
    @enforce precompile(Tuple{typeof(eltype),LinearAlgebra.Adjoint{Rational{Int}, Matrix{Rational{Int}}}})
    @enforce precompile(Tuple{typeof(isone),Matrix{Int32}})
    @enforce precompile(Tuple{typeof(hcat),Vector{Rational{Int}},LinearAlgebra.Adjoint{Rational{Int}, Matrix{Rational{Int}}}})
    let fbody = try __lookup_kwbody__(which(LinearAlgebra.rank, (Matrix{Int},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Float64,Float64,typeof(LinearAlgebra.rank),Matrix{Int},))
        end
    end

    # Base
    @enforce precompile(Tuple{Core.kwftype(typeof(Base.with_output_color)),NamedTuple{(:bold,), Tuple{Bool}},typeof(Base.with_output_color),Function,Symbol,IOContext{Base.TTY},String,Vararg{Any, N} where N})
    @enforce precompile(Tuple{Type{Base.IteratorSize},Base.Iterators.ProductIterator{Tuple{UnitRange{Int}, UnitRange{Int}}}})
    @enforce precompile(Tuple{Type{Dict{Symbol, Any}},Base.Iterators.Pairs{Symbol, Any, Tuple{Symbol, Symbol, Symbol}, NamedTuple{(:help, :metavar, :required), Tuple{String, String, Bool}}}})
    @enforce precompile(Tuple{Type{Dict{Symbol, Any}},Base.Iterators.Pairs{Symbol, Any, Tuple{Symbol, Symbol}, NamedTuple{(:help, :action), Tuple{String, Symbol}}}})
    @enforce precompile(Tuple{Type{Dict{Symbol, Any}},Base.Iterators.Pairs{Symbol, String, Tuple{Symbol, Symbol}, NamedTuple{(:help, :metavar), Tuple{String, String}}}})
    @enforce precompile(Tuple{Type{Tuple{Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}, Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}, Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}}},Vector{Tuple{Vector{Any}, Vector{Any}}}})
    for T in (Int32, Int64, Int128)
        @enforce precompile(Tuple{Type{SubArray},IndexLinear,Matrix{Rational{T}},Tuple{Base.Slice{Base.OneTo{Int}}, Int},Tuple{Bool}})
        @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, Type{Rational{T}}, Tuple{Matrix{Rational{Int}}}}})
    end
    @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(big), Tuple{Matrix{Rational{Int}}}}})
    @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(denominator), Tuple{Matrix{Rational{Int}}}}})
    @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Nothing, typeof(numerator), Tuple{Matrix{Rational{Int}}}}})
    @enforce precompile(Tuple{typeof(Base.deepcopy_internal),NTuple{9, BigFloat},IdDict{Any, Any}})
    @enforce precompile(Tuple{typeof(Base.display_error),Base.TTY,ErrorException,Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}})
    @enforce precompile(Tuple{typeof(Base.grow_to!),Dict{Symbol, Any},Tuple{Pair{Symbol, String}, Pair{Symbol, Bool}, Pair{Symbol, String}},Int})
    @enforce precompile(Tuple{typeof(Base.grow_to!),Dict{Symbol, Any},Tuple{Pair{Symbol, String}, Pair{Symbol, Symbol}},Int})
    @enforce precompile(Tuple{typeof(Base.print_to_string),Int128,Vararg{Any, N} where N})
    @enforce precompile(Tuple{typeof(Base.reducedim_init),Function,typeof(min),Matrix{Float64},Int})
    @enforce precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Int},Expr,Int})
    @enforce precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{Any}},Tuple{},Int})
    @enforce precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{Int}},Tuple{},Int})
    @enforce precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{}},Tuple{Int},Int})
    @enforce precompile(Tuple{typeof(Base.typed_hvcat),Type{BigFloat},Tuple{Int, Int, Int},BigFloat,Vararg{Number, N} where N})
    @enforce precompile(Tuple{typeof(Base.typed_hvcat),Type{Float64},Tuple{Int, Int, Int},Int,Vararg{Number, N} where N})
    @enforce precompile(Tuple{typeof(copy),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{2}, Tuple{Base.OneTo{Int}, Base.OneTo{Int}}, Type{Int}, Tuple{Matrix{Rational{BigInt}}}}})
    @enforce precompile(Tuple{typeof(maximum),Matrix{Int}})
    @enforce precompile(Tuple{typeof(minimum),Matrix{Int}})
    @enforce precompile(Tuple{typeof(push!),Vector{String},String,String,String})
    @enforce precompile(Tuple{typeof(setindex!),Dict{String, Tuple{Exception, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}},Tuple{ArgumentError, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}},String})
    @enforce precompile(Tuple{typeof(string),Int128,String,Vararg{Any, N} where N})
    @enforce precompile(Tuple{typeof(vcat),Vector{Expr},Vector{Expr}})
    @static if VERSION >= v"1.6-"
        let fbody = try __lookup_kwbody__(which(Base.print_within_stacktrace, (IOContext{Base.TTY},String,Vararg{Any, N} where N,))) catch missing end
            if !ismissing(fbody)
                @enforce precompile(fbody, (Symbol,Bool,typeof(Base.print_within_stacktrace),IOContext{Base.TTY},String,Vararg{Any, N} where N,))
            end
        end
    end
    let fbody = try __lookup_kwbody__(which(all, (Function,Vector{String},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Function,typeof(all),Function,Vector{String},))
        end
    end
    let fbody = try __lookup_kwbody__(which(any, (Function,Vector{AbstractString},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Function,typeof(any),Function,Vector{AbstractString},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sort!, (Vector{Int},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sort!),Vector{Int},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{Symbol},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{Symbol},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{Tuple{Int, Vector{Int}}},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{Tuple{Int, Vector{Int}}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{Vector{Int}},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{Vector{Int}},))
        end
    end

    # StaticArrays
    @enforce precompile(Tuple{Type{Vector{SVector{3, Float64}}},Vector{SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true}}})
    @enforce precompile(Tuple{Type{Vector{SVector{3, Int}}},Vector{Any}})
    @enforce precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{0},StaticArrays.StaticArrayStyle{2}})
    @enforce precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},StaticArrays.StaticArrayStyle{1}})
    @enforce precompile(Tuple{Type{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(+),Tuple{SVector{3, Int}, SVector{3, Int}}})
    @enforce precompile(Tuple{Type{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(-),Tuple{SVector{3, Rational{Int}}, SVector{3, Int}}})
    @enforce precompile(Tuple{Type{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(floor),Tuple{Base.RefValue{Type{Int}}, SVector{3, Rational{Int}}}})
    @enforce precompile(Tuple{Type{Dict{Int, SVector{3, Int}}}})
    @enforce precompile(Tuple{Type{Dict{SVector{3, Int}, Nothing}}})
    @enforce precompile(Tuple{Type{SMatrix{3, 3, BigFloat, 9}},NTuple{9, Float64}})
    @enforce precompile(Tuple{typeof(==),Vector{Tuple{Int, Int, SVector{3, Rational{Int}}}},Vector{Tuple{Int, Int, SVector{3, Rational{Int}}}}})
    @enforce precompile(Tuple{typeof(>),SVector{3, Int},SizedVector{3, Int, 1}})
    @enforce precompile(Tuple{typeof(>),SizedVector{3, Int, 1},SVector{3, Int}})
    @enforce precompile(Tuple{typeof(Base.Broadcast.broadcasted),Function,SVector{3, Rational{Int}},SVector{3, Int}})
    @enforce precompile(Tuple{typeof(StaticArrays.arithmetic_closure),Type{BigFloat}})
    @enforce precompile(Tuple{typeof(LinearAlgebra.generic_norm2),SubArray{BigFloat, 1, SMatrix{3, 3, BigFloat, 9}, Tuple{Base.Slice{SOneTo{3}}, Int}, true}})
    @enforce precompile(Tuple{typeof(StaticArrays._axes),Size{(3, 3)}})
    @enforce precompile(Tuple{typeof(StaticArrays._axes),Size{(3,)}})

    @enforce precompile(Tuple{Type{Size},Type{SubArray{Int, 1, Matrix{Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, true}}})
    @enforce precompile(Tuple{typeof(<),SVector{3, Int},SizedVector{3, Int, 1}})
    for T in (Int32, Int64, Int128, BigInt)
        @enforce precompile(Tuple{Type{Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}},Vector{Any}})
        @enforce precompile(Tuple{Type{SVector{3, Int}},Tuple{Rational{T}, Rational{T}, Rational{T}}})
        @enforce precompile(Tuple{Type{Dict{Int, Vector{SMatrix{3, 3, Rational{T}, 9}}}}})
        @enforce precompile(Tuple{Type{Dict{SVector{9, Rational{T}}, Nothing}}})
        @enforce precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.Style{Tuple}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(isempty),Tuple{Tuple{Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}}}})
        @enforce precompile(Tuple{typeof(<),SVector{3, Rational{T}},SVector{3, Rational{Int}}})
        @enforce precompile(Tuple{typeof(<),SVector{3, Rational{Int}},SVector{3, Rational{T}}})
        @enforce precompile(Tuple{typeof(==),SVector{3, Rational{Int}},SVector{3, Rational{T}}})
        @enforce precompile(Tuple{typeof(==),SMatrix{3, 3, T, 9},LinearAlgebra.UniformScaling{Bool}})
        @enforce precompile(Tuple{typeof(LinearAlgebra.generic_matmatmul!),Matrix{Rational{widen(T)}},Char,Char,SizedMatrix{3, 3, Rational{widen(T)}, 2},Matrix{Rational{T}},LinearAlgebra.MulAddMul{true, true, Bool, Bool}})
        @enforce precompile(Tuple{typeof(LinearAlgebra.dot),SubArray{BigFloat, 1, SMatrix{3, 3, BigFloat, 9}, Tuple{Base.Slice{SOneTo{3}}, T}, true},SubArray{BigFloat, 1, SMatrix{3, 3, BigFloat, 9}, Tuple{Base.Slice{SOneTo{3}}, T}, true}})
        @enforce precompile(Tuple{typeof(Base.Broadcast.broadcasted),Base.Broadcast.Style{Tuple},Function,Tuple{Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}, Vector{Tuple{Int, Int, SVector{3, Rational{T}}}}}})
        @enforce precompile(Tuple{typeof(Base.Broadcast.materialize!),Base.Broadcast.DefaultArrayStyle{1},SubArray{Rational{T}, 1, Matrix{Rational{T}}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true},Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(-), Tuple{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(+), Tuple{SVector{3, Rational{T}}, SVector{3, Int}}}, SVector{3, Rational{T}}}}})
        @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(-), Tuple{SVector{3, Rational{T}}, SVector{3, Int}}}})
        @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(floor), Tuple{Base.RefValue{Type{Int}}, SVector{3, Rational{T}}}}})
        @enforce precompile(Tuple{typeof(append!),Vector{SMatrix{3, 3, Rational{T}, 9}},Vector{SMatrix{3, 3, Rational{T}, 9}}})
        @enforce precompile(Tuple{typeof(append!),Vector{SMatrix{3, 3, Rational{T}, 9}},Vector{SMatrix{3, 3, Rational{widen(T)}, 9}}})
        let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Rational{T}}},))) catch missing end
            if !ismissing(fbody)
                @enforce precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Rational{T}}},))
            end
        end
        let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Rational{Int128}}},))) catch missing end
            if !ismissing(fbody)
                @enforce precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Rational{T}}},))
            end
        end
        @enforce precompile(Tuple{typeof(sort!),Vector{Int},Base.Sort.QuickSortAlg,Base.Order.Perm{Base.Order.ForwardOrdering, Vector{SVector{3, Rational{T}}}}})
    end
    @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(//), Tuple{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(round), Tuple{Base.RefValue{Type{Int128}}, Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(*), Tuple{Int128, SVector{3, Float64}}}}}, Int128}}})
    @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(//), Tuple{Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(round), Tuple{Base.RefValue{Type{BigInt}}, Base.Broadcast.Broadcasted{StaticArrays.StaticArrayStyle{1}, Nothing, typeof(*), Tuple{BigInt, SVector{3, Float64}}}}}, BigInt}}})
    @enforce precompile(Tuple{typeof(Base._foldl_impl),Base.BottomRF{typeof(hcat)},Base.ReshapedArray{Int, 2, SVector{3, Int}, Tuple{}},Set{SVector{3, Int}}})
    @enforce precompile(Tuple{typeof(Base.setindex_widen_up_to),Vector{Tuple{StaticArrays.Dynamic}},Tuple{Int},Int})
    @enforce precompile(Tuple{typeof(StaticArrays._axes),Size{(9,)}})
    @enforce precompile(Tuple{typeof(setindex!),Dict{Int, SVector{3, Int}},SVector{3, Int},Int})
    @enforce precompile(Tuple{typeof(setindex!),Dict{SVector{3, Int}, Nothing},Nothing,SVector{3, Int}})
    @enforce precompile(Tuple{typeof(which(StaticArrays._getindex,(AbstractArray,Tuple{Vararg{Size, N} where N},Any,)).generator.gen),Any,Any,Any,Any})
    @enforce precompile(Tuple{typeof(which(StaticArrays._map,(Any,Vararg{AbstractArray, N} where N,)).generator.gen),Any,Any,Any})
    @enforce precompile(Tuple{typeof(which(StaticArrays.combine_sizes,(Tuple{Vararg{Size, N} where N},)).generator.gen),Any,Any})
    let fbody = try __lookup_kwbody__(which(LinearAlgebra.rank, (Base.ReshapedArray{Int, 2, SVector{3, Int}, Tuple{}},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Float64,Float64,typeof(LinearAlgebra.rank),Base.ReshapedArray{Int, 2, SVector{3, Int}, Tuple{}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Float64}},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Float64}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Rational{BigInt}}},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Rational{BigInt}}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SVector{3, Rational{Bool}}},))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Base.Sort.QuickSortAlg,Function,Function,Nothing,Base.Order.ForwardOrdering,typeof(sortperm),Vector{SVector{3, Rational{Bool}}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(sprint, (Function,SVector{3, Int},Vararg{Any, N} where N,))) catch missing end
        if !ismissing(fbody)
            @enforce precompile(fbody, (Nothing,Int,typeof(sprint),Function,SVector{3, Int},Vararg{Any, N} where N,))
        end
    end
    @enforce precompile(Tuple{typeof(sort!),Vector{Int},Base.Sort.QuickSortAlg,Base.Order.Perm{Base.Order.ForwardOrdering, Vector{SVector{3, Float64}}}})


    # Chemfiles
    @enforce precompile(Tuple{Type{Chemfiles.Atom},Chemfiles.Frame,Int})
    @enforce precompile(Tuple{Type{Chemfiles.Atom},String})
    @enforce precompile(Tuple{Type{Chemfiles.Frame}})
    @enforce precompile(Tuple{Type{Chemfiles.Topology},Chemfiles.Frame})
    @enforce precompile(Tuple{Type{Chemfiles.UnitCell},BigFloat,BigFloat,BigFloat,BigFloat,BigFloat,BigFloat})
    @enforce precompile(Tuple{Type{Chemfiles.UnitCell},Chemfiles.Frame})
    @enforce precompile(Tuple{typeof(Chemfiles.__free),Chemfiles.CxxPointer{Chemfiles.lib.CHFL_ATOM}})
    @enforce precompile(Tuple{typeof(Chemfiles.__free),Chemfiles.CxxPointer{Chemfiles.lib.CHFL_CELL}})
    @enforce precompile(Tuple{typeof(Chemfiles.__free),Chemfiles.CxxPointer{Chemfiles.lib.CHFL_TOPOLOGY}})
    @enforce precompile(Tuple{typeof(Chemfiles.__strip_null),String})
    @enforce precompile(Tuple{typeof(Chemfiles.add_atom!),Chemfiles.Frame,Chemfiles.Atom,Vector{Float64},Vector{Float64}})
    @enforce precompile(Tuple{typeof(Chemfiles.bonds),Chemfiles.Topology})
    @enforce precompile(Tuple{typeof(Chemfiles.matrix),Chemfiles.UnitCell})
    @enforce precompile(Tuple{typeof(Chemfiles.positions),Chemfiles.Frame})
    @enforce precompile(Tuple{typeof(Chemfiles.set_type!),Chemfiles.Atom,String})

        # CrystalNets
    # @enforce precompile(Tuple{CrystalNets.var"#730#threadsfor_fun#121"{Tuple{Symbol}, String, Dict{String, String}, Dict{String, Tuple{Exception, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}}, Dict{String, String}, Vector{String}}})
    for T in (Int32, Int64, Int128, BigInt)
        # @enforce precompile(Tuple{CrystalNets.var"#670#threadsfor_fun#114"{CrystalNet{Rational{T}}, Vector{Int}, Base.Threads.SpinLock, Vector{Pair{Int, Tuple{Matrix{Rational{T}}, Vector{Int}}}}, Int, DataType, Vector{Int}}})
        # @enforce precompile(Tuple{CrystalNets.var"#685#threadsfor_fun#116"{Rational{T}, Dict{Int, Vector{SMatrix{3, 3, Rational{T}, 9}}}, Base.Threads.SpinLock, DataType, Vector{Pair{Int, Tuple{Matrix{Rational{T}}, Vector{Int}}}}}})

        @enforce precompile(Tuple{Type{CrystalNet{Rational{T}}},CrystalNets.Cell{Rational{Int}},Vector{Symbol},PeriodicGraph3D,Matrix{Rational{T}}})
        @enforce precompile(Tuple{typeof(CrystalNets.isrank3),Matrix{Rational{T}}})
        @enforce precompile(Tuple{typeof(CrystalNets.minimize),CrystalNet{Rational{T}}})
        @enforce precompile(Tuple{typeof(CrystalNets.topological_key),CrystalNet{Rational{T}}})
        @enforce precompile(Tuple{typeof(topological_genome),CrystalNet{Rational{T}}})

    end
    @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(CrystalNets.parsestrip), Tuple{Vector{String}}}})
    @enforce precompile(Tuple{typeof(CrystalNets.invalid_input_error),String,ErrorException,Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}})
    @enforce precompile(Tuple{typeof(CrystalNets.julia_main)})

    @enforce precompile(Tuple{typeof(CrystalNets.parse_chemfile),String})
    @enforce precompile(Tuple{typeof(CrystalNets.parse_chemfile),String,Bool})
    @enforce precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes},typeof(CrystalNets.parsestrip),Tuple{Vector{String}}})
    @enforce precompile(Tuple{typeof(Base.Broadcast.materialize),Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(CrystalNets.parsestrip), Tuple{Vector{String}}}})
    @enforce precompile(Tuple{typeof(CrystalNets.isrank3),Matrix{Rational{Int32}}})
    @enforce precompile(Tuple{typeof(parse_chemfile),String})
    @enforce precompile(Tuple{typeof(recognize_topology),PeriodicGraph3D})
    for T in (Nothing, CrystalNets.Clusters)
        @enforce precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}})
        @enforce precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, Nothing})
        @enforce precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, String})
        @enforce precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, Nothing, String})
        @enforce precompile(Tuple{typeof(CrystalNets.ifexport), CrystalNets.Crystal{T}, String, String})
    end

    @enforce precompile(Tuple{typeof(convert), Type{Any}, Vector{Union{String,Vector{String}}}})
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
    p_pos = Vector{SVector{3,Float64}}
    pos{D,T} = SVector{D,Rational{T}}
    ppos{D,T} = Vector{SVector{D,Rational{T}}}
    ofs{D} = SVector{D,Int}
    pofs{D} = Vector{SVector{D,Int}}
    pge{D,T} = PeriodicGraphEmbedding{D,T}
    symmgroup{T} = SymmetryGroup3D{Rational{T}}
    cnet{D,T} = CrystalNet{D,Rational{T}}
    mat = SMatrix{3,3,BigFloat,9}
    fmat = SMatrix{3,3,Float64,9}
    sbus = Vector{Vector{PeriodicVertex3D}}
    kinds = CrystalNets.ClusterKinds
    structure = StructureType._StructureType
    clustering = Clustering._Clustering
    neighs = Vector{PeriodicVertex3D}
    edgs = Vector{PeriodicEdge3D}
    cell = PeriodicGraphEmbeddings.Cell{Rational{Int}}
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
        @enforce precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Int,Function})
        @enforce precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Int,Function})
        @enforce precompile(Tuple{typeof(-),CrystalNets.Modulos.Modulo{P, Int32},CrystalNets.Modulos.Modulo{P, Int32}})
        @enforce precompile(Tuple{typeof(/),Int,CrystalNets.Modulos.Modulo{P, Int32}})
        @enforce precompile(Tuple{typeof(==),Matrix{CrystalNets.Modulos.Modulo{P, Int32}},Matrix{Int}})
        @enforce precompile(Tuple{typeof(Base._unsafe_copyto!),Matrix{CrystalNets.Modulos.Modulo{P, Int32}},Int,Matrix{Int},Int,Int})
        @enforce precompile(Tuple{typeof(Base._unsafe_getindex),IndexLinear,Matrix{CrystalNets.Modulos.Modulo{P, Int32}},Int,Base.Slice{Base.OneTo{Int}}})
        @enforce precompile(Tuple{typeof(Base.copyto_unaliased!),IndexLinear,SubArray{CrystalNets.Modulos.Modulo{P, Int32}, 1, Matrix{CrystalNets.Modulos.Modulo{P, Int32}}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, true},IndexLinear,Vector{CrystalNets.Modulos.Modulo{P, Int32}}})
        @enforce precompile(Tuple{typeof(LinearAlgebra.mul!),Matrix{CrystalNets.Modulos.Modulo{P, Int32}},SparseArrays.SparseMatrixCSC{Int, Int},Matrix{CrystalNets.Modulos.Modulo{P, Int32}},Bool,Bool})
        @enforce precompile(Tuple{typeof(SparseArrays._setindex_scalar!),SparseArrays.SparseMatrixCSC{CrystalNets.Modulos.Modulo{P, Int32}, Int},Int,Int,Int})
        @enforce precompile(Tuple{typeof(SparseArrays.sparse!),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Int,typeof(+),Vector{Int},Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}}})
        @enforce precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Int,Function})
        @enforce precompile(Tuple{typeof(SparseArrays.sparse_check_length),String,Vector{CrystalNets.Modulos.Modulo{P, Int32}},Int,Type})
        @enforce precompile(Tuple{typeof(getindex),SparseArrays.SparseMatrixCSC{CrystalNets.Modulos.Modulo{P, Int32}, Int},UnitRange{Int},UnitRange{Int}})
    end


    # output.jl
    @enforce precompile(Tuple{typeof(CrystalNets.export_dataline), Base.IOStream, String})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.export_vtf), String, cryst, Int, Bool})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.export_vtf), String, crystclust, Int, Bool})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.export_vtf), String, cryst, Int})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.export_vtf), String, crystclust, Int})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.export_vtf), String, cryst})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.export_vtf), String, crystclust})
    for D in 1:3
        for T in inttypes
            @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.export_vtf), String, cnet{D,T}, Int, Bool})
            @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.export_vtf), String, cnet{D,T}, Int})
            @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.export_vtf), String, cnet{D,T}})
        end
    end
    @enforce precompile(Tuple{typeof(CrystalNets.export_cif), String, cryst})
    @enforce precompile(Tuple{typeof(CrystalNets.export_cif), String, crystclust})
    @enforce precompile(Tuple{typeof(CrystalNets.export_cif), String, cif})
    @enforce precompile(Tuple{typeof(CrystalNets.export_cgd), String, cryst})
    @enforce precompile(Tuple{typeof(CrystalNets.export_cgd), String, crystclust})
    for D in 1:3
        @enforce precompile(Tuple{typeof(CrystalNets.export_cgd), String, PeriodicGraph{D}})
        for T in inttypes
            @enforce precompile(Tuple{typeof(CrystalNets.export_cgd), String, cnet{D,T}})
        end
    end
    @enforce precompile(Tuple{typeof(CrystalNets.export_attributions), crystclust, String})
    @enforce precompile(Tuple{typeof(CrystalNets.export_attributions), crystclust})
    @enforce precompile(Tuple{typeof(CrystalNets.export_arc), String, Dict{String,String}})
    @enforce precompile(Tuple{typeof(CrystalNets.export_arc), String, Nothing})
    @enforce precompile(Tuple{typeof(CrystalNets.export_arc), String})

    # utils.jl
    @enforce precompile(Tuple{typeof(CrystalNets.recursive_readdir!), Vector{String}, String, String})
    @enforce precompile(Tuple{typeof(CrystalNets.recursive_readdir), String})
    @enforce precompile(Tuple{typeof(CrystalNets.tmpexportname), String, String, String, String})
    @enforce precompile(Tuple{typeof(CrystalNets.tmpexportname), String, String, Nothing, String})
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
    @enforce precompile(Tuple{typeof(CrystalNets.string_atomtype), Symbol})
    @enforce precompile(Tuple{typeof(CrystalNets.representative_atom), Symbol, Int})
    for T in inttypes
        for D in 1:3
            @enforce precompile(Tuple{typeof(CrystalNets.issingular), smat{D,T,D*D}})
            @enforce precompile(Tuple{typeof(CrystalNets.angle), pos{D,T}, pos{D, T}})
            @enforce precompile(Tuple{typeof(CrystalNets.dihedral), pos{D,T},pos{D,T},pos{D,T}})
        end
        @enforce precompile(Tuple{typeof(CrystalNets.isrank3), Matrix{Rational{T}}})
        @enforce precompile(Tuple{typeof(CrystalNets.isrank2), Matrix{Rational{T}}})
        @enforce precompile(Tuple{typeof(CrystalNets.back_to_unit), Rational{T}})
    end
    @enforce precompile(Tuple{typeof(CrystalNets.nextword), String, Int})
    @enforce precompile(Tuple{typeof(CrystalNets.isinterrupt), Any})
    @enforce precompile(Tuple{typeof(CrystalNets.isoverfloworinexact), Any})

    # options.jl
    @enforce precompile(Tuple{typeof(CrystalNets.clustering_from_num), Int})
    @enforce precompile(Tuple{typeof(CrystalNets.clustering_from_num), Int8})
    @enforce precompile(Tuple{typeof(CrystalNets.clustering_from_symb), Symbol})
    @enforce precompile(Tuple{typeof(parse), Type{clustering}, String})
    @enforce precompile(Tuple{typeof(parse), Type{clustering}, SubString{String}})
    @enforce precompile(Tuple{typeof(getindex), kinds, Int})
    @enforce precompile(Tuple{typeof(getindex), kinds, Symbol})
    @enforce precompile(Tuple{typeof(CrystalNets.getmetal), kinds})
    @enforce precompile(Tuple{typeof(length), kinds})
    @enforce precompile(Tuple{typeof(CrystalNets.ifbooltempdirorempty), String})
    @enforce precompile(Tuple{typeof(CrystalNets.ifbooltempdirorempty), Bool})
    @enforce precompile(Tuple{Type{CrystalNets.Options}})
    @enforce precompile(Tuple{Type{CrystalNets.Options}, CrystalNets.Options})

    # specialsolver.jl
    for Ti in (BigRational, (modulo{P} for P in primes)...)
        @enforce precompile(Tuple{typeof(CrystalNets.rational_lu!), SparseMatrixCSC{Ti,Int}, Vector{Int}, Bool})
        @enforce precompile(Tuple{typeof(CrystalNets.rational_lu!), SparseMatrixCSC{Ti,Int}, Vector{Int}})
    end
    for P in primes
        @enforce precompile(Tuple{typeof(CrystalNets.rational_lu), SparseMatrixCSC{modulo{P},Int}, Bool, Type{modulo{P}}})
        @enforce precompile(Tuple{typeof(CrystalNets.rational_lu), SparseMatrixCSC{modulo{P},Int}, Bool, Type{BigRational}})
        @enforce precompile(Tuple{typeof(CrystalNets.rational_lu), SparseMatrixCSC{Int,Int}, Bool})
    end
    for Ti in (Rational{BigInt}, (modulo{P} for P in primes)...)
        @enforce precompile(Tuple{typeof(CrystalNets.forward_substitution!), SparseMatrixCSC{Ti,Int}, Matrix{Ti}})
        @enforce precompile(Tuple{typeof(CrystalNets.backward_substitution!), SparseMatrixCSC{Ti,Int}, Matrix{Ti}})
    end
    @static if VERSION < v"1.8-"
        @enforce precompile(Tuple{typeof(CrystalNets.linsolve!), LU{Rational{BigInt},SparseMatrixCSC{Rational{BigInt},Int}}, Matrix{Rational{BigInt}}})
        for P in primes
            @enforce precompile(Tuple{typeof(CrystalNets.linsolve!), LU{modulo{P},SparseMatrixCSC{modulo{P},Int}}, Matrix{Int}})
        end
    else
        @enforce precompile(Tuple{typeof(CrystalNets.linsolve!), LU{Rational{BigInt},SparseMatrixCSC{Rational{BigInt},Int},Base.OneTo{Int}}, Matrix{Rational{BigInt}}})
        for P in primes
            @enforce precompile(Tuple{typeof(CrystalNets.linsolve!), LU{modulo{P},SparseMatrixCSC{modulo{P},Int},Base.OneTo{Int}}, Matrix{Int}})
        end
    end
    for T in (Int64, Int128, BigInt)
        @enforce precompile(Tuple{typeof(CrystalNets.copyuntil), Int, Matrix{Rational{T}}, Type{Rational{T}}})
        @enforce precompile(Tuple{typeof(CrystalNets._inner_dixon_p!), Vector{Int}, Matrix{Rational{T}}, BigInt, Matrix{BigInt}, BigInt, BigInt})
    end
    for N in 1:3
        @enforce precompile(Tuple{typeof(CrystalNets.rational_solve), Val{N}, SparseMatrixCSC{Int,Int}, Matrix{Int}})
        for P in primes
            @static if VERSION < v"1.8-"
                @enforce precompile(Tuple{typeof(CrystalNets.dixon_p), Val{N}, SparseMatrixCSC{Int,Int}, LU{modulo{P},SparseMatrixCSC{modulo{P},Int}}, Matrix{Int}})
            else
                @enforce precompile(Tuple{typeof(CrystalNets.dixon_p), Val{N}, SparseMatrixCSC{Int,Int}, LU{modulo{P},SparseMatrixCSC{modulo{P},Int},Base.OneTo{Int}}, Matrix{Int}})
            end
        end
        @enforce precompile(Tuple{typeof(CrystalNets.dixon_solve), Val{N}, SparseMatrixCSC{Int,Int}, Matrix{Int}})
    end

    # types.jl 1/2
    @enforce precompile(Tuple{Type{cif}, Dict{String,Union{String,Vector{String}}}, cell, Vector{Int}, types, Matrix{Float64}, bondlist})
    @enforce precompile(Tuple{typeof(CrystalNets.keepinbonds), bondlist, Vector{Int}})
    @enforce precompile(Tuple{typeof(CrystalNets.add_to_bondlist!), Vector{Tuple{Int,Float32}}, Int, Float32})
    @enforce precompile(Tuple{typeof(CrystalNets.get_bondlist), Vector{Tuple{Int,Float32}}, Int})
    @enforce precompile(Tuple{typeof(CrystalNets.sortprune_bondlist!), Vector{Tuple{Int,Float32}}})
    @enforce precompile(Tuple{typeof(CrystalNets.remove_partial_occupancy), cif})
    @enforce precompile(Tuple{typeof(CrystalNets.prune_collisions), cif})
    @enforce precompile(Tuple{typeof(CrystalNets.expand_symmetry), cif})
    @enforce precompile(Tuple{typeof(CrystalNets.edges_from_bonds), bondlist, fmat, p_pos})
    @enforce precompile(Tuple{Type{clust}, sbus, Vector{Int}, Vector{Int}, pofs{3}, BitVector})
    @enforce precompile(Tuple{Type{clust}, Int})
    @enforce precompile(Tuple{typeof(isempty), clust})
    @enforce precompile(Tuple{typeof(getindex), clust, Vector{Int}})
    @enforce precompile(Tuple{typeof(getindex), clust, Base.OneTo{Int}})
    @enforce precompile(Tuple{typeof(getindex), clust, Base.UnitRange{Int}})
    @enforce precompile(Tuple{Type{cryst}, cell, types, p_pos, graph, opts})
    @enforce precompile(Tuple{Type{crystclust}, cell, types, clust, p_pos, graph, opts})
    @enforce precompile(Tuple{Type{CrystalNets.Crystal}, pge{3,Float64}, types, Nothing, opts})
    @enforce precompile(Tuple{Type{CrystalNets.Crystal}, pge{3,Float64}, types, clust, opts})
    @enforce precompile(Tuple{Type{CrystalNets.Crystal}, cryst})
    @enforce precompile(Tuple{Type{CrystalNets.Crystal}, crystclust})
    @enforce precompile(Tuple{typeof(==), cryst, cryst})
    @enforce precompile(Tuple{typeof(==), crystclust, crystclust})
    @enforce precompile(Tuple{Type{cryst}, cryst})
    @enforce precompile(Tuple{Type{cryst}, crystclust})
    @enforce precompile(Tuple{Type{crystclust}, cryst, clust})
    @enforce precompile(Tuple{Type{crystclust}, crystclust, clust})
    @enforce precompile(Tuple{Type{crystclust}, cryst})
    @enforce precompile(Tuple{Type{crystclust}, crystclust})
    @enforce precompile(Tuple{typeof(CrystalNets.remove_metal_cluster_bonds!), graph, types, opts})
    @enforce precompile(Tuple{typeof(CrystalNets.trimmed_crystal), cryst})
    @enforce precompile(Tuple{typeof(getindex), cryst, Vector{Int}})
    @enforce precompile(Tuple{typeof(getindex), cryst, Base.OneTo{Int}})
    @enforce precompile(Tuple{typeof(getindex), cryst, Base.UnitRange{Int}})
    @enforce precompile(Tuple{typeof(getindex), crystclust, Vector{Int}})
    @enforce precompile(Tuple{typeof(getindex), crystclust, Base.OneTo{Int}})
    @enforce precompile(Tuple{typeof(getindex), crystclust, Base.UnitRange{Int}})
    for D in 1:3
        @enforce precompile(Tuple{typeof(CrystalNets.equilibrium), PeriodicGraph{D}})
        @enforce precompile(Tuple{typeof(CrystalNets.trim_topology), PeriodicGraph{D}})
        for T in inttypes
            @enforce precompile(Tuple{Type{cnet{D,T}}, cell, types, ppos{D,T}, PeriodicGraph{D}, opts})
        end
    end
    for D in 1:3
        for T in inttypes
            for D2 in 1:3
                for T2 in inttypes
                    @enforce precompile(Tuple{Type{cnet{D,T}}, cnet{D2,T2}})
                end
                @enforce precompile(Tuple{Type{CrystalNet{D}}, cnet{D2,T}})
            end
            @enforce precompile(Tuple{typeof(ndims), cnet{D,T}})
            @enforce precompile(Tuple{Type{cnet{D,T}}, PeriodicGraphEmbedding{D,Rational{T}}, types, opts})
            @enforce precompile(Tuple{Type{cnet{D,T}}, cell, opts})
            @enforce precompile(Tuple{typeof(show), Base.TTY, cnet{D,T}})
            @enforce precompile(Tuple{typeof(show), Base.IOStream, cnet{D,T}})
        end
        @enforce precompile(Tuple{Type{CrystalNet{D}}, cell, opts})
    end

    # clustering.jl
    @enforce precompile(Tuple{typeof(CrystalNets.regroup_sbus), graph, Vector{Int}, Vector{Int}})
    @enforce precompile(Tuple{typeof(CrystalNets._trim_monovalent!), graph})
    @enforce precompile(Tuple{typeof(CrystalNets.trim_monovalent), cryst})
    @enforce precompile(Tuple{typeof(CrystalNets.trim_monovalent), crystclust})
    @enforce precompile(Tuple{typeof(CrystalNets.delete_target_from_list!), sbus, Int})
    @enforce precompile(Tuple{typeof(CrystalNets.is_paddlewheel_candidate!), Vector{Union{Missing, Tuple{Symbol, PeriodicVertex3D}}}, sbus, Int, types, Set{Int}})
    # @enforce precompile(Tuple{typeof(CrystalNets.bond_carboxylic_acid!), graph, types})
    @enforce precompile(Tuple{typeof(CrystalNets.regroup_paddlewheel!), graph, clust, types, Set{Int}})
    @enforce precompile(Tuple{typeof(CrystalNets.split_sbu!), clust, graph, Int, Vector{Int}})
    @enforce precompile(Tuple{typeof(CrystalNets.reclassify!), clust, Vector{Int}, Int, graph, types, Dict{Symbol,Int}, Int})
    @enforce precompile(Tuple{typeof(CrystalNets.add_to_newclass!), Vector{Int}, graph, clust, Int, Int, types, Set{Symbol}})
    @enforce precompile(Tuple{typeof(CrystalNets.add_to_merge_or_newclass!), Vector{Int}, Vector{Tuple{ofs{3},Int}}, graph, clust, Set{Int}, Int, Int})
    @enforce precompile(Tuple{typeof(CrystalNets.small_cycles_around!), Dict{PeriodicEdge3D,Vector{Int}}, Int, graph, p_pos, mat, Int, PeriodicVertex3D, Vector{Int}, Set{Int}, Int})
    @enforce precompile(Tuple{typeof(CrystalNets.detect_organiccycles), Vector{Int}, graph, p_pos, mat, Int, Set{Int}})
    @enforce precompile(Tuple{typeof(CrystalNets.group_cycle), Vector{Set{Int}}, types, graph})
    @enforce precompile(Tuple{typeof(CrystalNets.identify_metallic_type), Symbol, kinds, Int})
    @enforce precompile(Tuple{typeof(CrystalNets.find_sbus), cryst, kinds})
    @enforce precompile(Tuple{typeof(CrystalNets.find_sbus), crystclust, kinds})
    @enforce precompile(Tuple{typeof(CrystalNets._split_this_sbu!), Vector{Int}, graph, Int, types, Symbol, sbus})
    @enforce precompile(Tuple{typeof(CrystalNets.split_special_sbu!), graph, sbus, graph, types, Bool})
    @enforce precompile(Tuple{typeof(CrystalNets.split_O_vertices), cryst})
    @enforce precompile(Tuple{typeof(CrystalNets.split_O_vertices), crystclust})
    @enforce precompile(Tuple{typeof(CrystalNets.identify_clustering), cryst, structure, clustering})
    @enforce precompile(Tuple{typeof(CrystalNets.identify_clustering), crystclust, structure, clustering})
    @enforce precompile(Tuple{typeof(CrystalNets.order_atomtype), Symbol})
    @enforce precompile(Tuple{typeof(CrystalNets._collapse_clusters), cryst, clust, Bool, Bool})
    @enforce precompile(Tuple{typeof(CrystalNets._collapse_clusters), cryst, clust, Bool})
    @enforce precompile(Tuple{typeof(CrystalNets._collapse_clusters), cryst, clust})
    @enforce precompile(Tuple{typeof(CrystalNets._find_clusters), cryst, Bool, Bool})
    @enforce precompile(Tuple{typeof(CrystalNets._find_clusters), crystclust, Bool, Bool})
    @enforce precompile(Tuple{typeof(CrystalNets.find_clusters), cryst})
    @enforce precompile(Tuple{typeof(CrystalNets.find_clusters), crystclust})
    @enforce precompile(Tuple{typeof(CrystalNets.collapse_clusters), cryst})
    @enforce precompile(Tuple{typeof(CrystalNets.collapse_clusters), crystclust})
    @enforce precompile(Tuple{typeof(CrystalNets.update_new_edgs!), Dict{PeriodicEdge3D,Bool}, neighs, BitMatrix, BitVector, Vector{Int}})
    @enforce precompile(Tuple{typeof(sort), SVector{4,Int}})
    @enforce precompile(Tuple{typeof(CrystalNets.edges_of_convex_hull), neighs, Int, fmat, p_pos, BitVector, Set{SVector{4,Int}}, BitMatrix})
    @enforce precompile(Tuple{typeof(CrystalNets.edges_of_convex_hull), neighs, Int, fmat, p_pos, BitVector, Set{SVector{4,Int}}})
    @enforce precompile(Tuple{typeof(CrystalNets.regroup_toremove), cryst, Vector{Int}, Vector{Vector{Int}}, String})
    @enforce precompile(Tuple{typeof(CrystalNets.pem_to_pe), cryst})
    @enforce precompile(Tuple{typeof(CrystalNets.regroup_vmap), cryst, Vector{Int}, Vector{Int}, String})
    @enforce precompile(Tuple{typeof(CrystalNets.allnodes_to_singlenodes), cryst})

    # guessbonds.jl
    @enforce precompile(Tuple{typeof(CrystalNets.guess_bonds), p_pos, types, fmat, opts})

    # types.jl, 2/2
    @enforce precompile(Tuple{typeof(CrystalNets.separate_components), cryst})
    @enforce precompile(Tuple{typeof(CrystalNets.separate_components), crystclust})
    for D in 1:3
        @enforce precompile(Tuple{typeof(CrystalNets._collect_net!), Vector{CrystalNet{D}}, Dict{graph,Int}, Int, cryst, clustering})
        @enforce precompile(Tuple{typeof(CrystalNets.collect_nets), Vector{cryst}, Val{D}})
    end
    @enforce precompile(Tuple{Type{unets}, Vector{Tuple{Vector{Int},Vector{CrystalNet1D}}},
                                  Vector{Tuple{Vector{Int},Vector{CrystalNet2D}}},
                                  Vector{Tuple{Vector{Int},Vector{CrystalNet3D}}}})
    @enforce precompile(Tuple{Type{unets}})
    @enforce precompile(Tuple{typeof(CrystalNets._repeatgroups!), Expr, Int})
    @enforce precompile(Tuple{Type{unets}, cryst})
    @enforce precompile(Tuple{Type{unets}, crystclust})
    @enforce precompile(Tuple{typeof(CrystalNets.__warn_nonunique), Int})
    @enforce precompile(Tuple{typeof(CrystalNets.__throw_interpenetrating), Int})
    @enforce precompile(Tuple{typeof(CrystalNets.__throw_multiplenets), Int})
    @enforce precompile(Tuple{Type{CrystalNet}, cryst})
    @enforce precompile(Tuple{Type{CrystalNet}, crystclust})
    for D in 1:3
        for D2 in 1:3
            @enforce precompile(Tuple{Type{CrystalNet{D}}, cell, types, PeriodicGraph{D2}, opts})
        end
        @enforce precompile(Tuple{Type{CrystalNet{D}}, unets})
        @enforce precompile(Tuple{Type{CrystalNet{D}}, cell, types, PeriodicGraph{D}, opts})
        @enforce precompile(Tuple{Type{unets}, PeriodicGraph{D}, opts})
        @enforce precompile(Tuple{Type{unets}, Vector{PeriodicEdge{D}}, opts})
        @enforce precompile(Tuple{Type{unets}, PeriodicGraph{D}})
        @enforce precompile(Tuple{Type{unets}, Vector{PeriodicEdge{D}}})
        @enforce precompile(Tuple{Type{CrystalNet}, PeriodicGraph{D}, opts})
        @enforce precompile(Tuple{Type{CrystalNet}, PeriodicGraph{D}})
    end
    @enforce precompile(Tuple{Type{CrystalNet}, unets})
    @enforce precompile(Tuple{Type{unets}, String, opts})
    for D in 1:3
        @enforce precompile(Tuple{Type{CrystalNet{D}}, PeriodicGraph{D}, opts})
        @enforce precompile(Tuple{Type{CrystalNet{D}}, PeriodicGraph{D}})
        @enforce precompile(Tuple{Type{genome}, PeriodicGraph{D}, Nothing, Bool, String})
        @enforce precompile(Tuple{Type{genome}, PeriodicGraph{D}, String, Bool, String})
        @enforce precompile(Tuple{Type{genome}, PeriodicGraph{D}, Nothing, Bool})
        @enforce precompile(Tuple{Type{genome}, PeriodicGraph{D}, String, Bool})
        @enforce precompile(Tuple{Type{genome}, PeriodicGraph{D}, Nothing})
        @enforce precompile(Tuple{Type{genome}, PeriodicGraph{D}, String})
    end
    @enforce precompile(Tuple{Type{genome}, String})
    @enforce precompile(Tuple{Type{genome}})
    @enforce precompile(Tuple{typeof(==), genome, genome})
    @enforce precompile(Tuple{typeof(hash), genome, UInt})
    @enforce precompile(Tuple{typeof(show), Base.TTY, genome})
    @enforce precompile(Tuple{typeof(show), Base.IOStream, genome})
    @enforce precompile(Tuple{typeof(parse), Type{genome}, String})
    @enforce precompile(Tuple{typeof(parse), Type{genome}, SubString{String}})
    @enforce precompile(Tuple{Type{result}, SizedVector{8,TopologicalGenome,Vector{TopologicalGenome}}, MVector{8,Int8}, Vector{Int8}})
    @enforce precompile(Tuple{Type{PeriodicGraph}, TopologicalGenome})
    for D in 1:3
        @enforce precompile(Tuple{Type{PeriodicGraph{D}}, TopologicalGenome})
    end
    @enforce precompile(Tuple{Type{result}})
    @enforce precompile(Tuple{Type{result}, Vector{Tuple{clustering,Union{clustering,genome}}}})
    @enforce precompile(Tuple{typeof(==), result, result})
    @enforce precompile(Tuple{typeof(hash), result, UInt})
    @enforce precompile(Tuple{typeof(get), CrystalNets.Returns{Nothing}, result, clustering})
    @enforce precompile(Tuple{typeof(get), result, clustering, Nothing})
    @enforce precompile(Tuple{typeof(getindex), result, clustering})
    @enforce precompile(Tuple{typeof(getindex), result, Symbol})
    @enforce precompile(Tuple{typeof(getindex), result, Int})
    @enforce precompile(Tuple{typeof(setindex!), result, genome, clustering})
    @enforce precompile(Tuple{typeof(setindex!), result, clustering, clustering})
    @enforce precompile(Tuple{Type{result}, String})
    @enforce precompile(Tuple{typeof(setindex!), result, Nothing, clustering})
    @enforce precompile(Tuple{typeof(show), Base.TTY, result})
    @enforce precompile(Tuple{typeof(show), Base.IOStream, result})
    @enforce precompile(Tuple{typeof(iterate), result, Int})
    @enforce precompile(Tuple{typeof(iterate), result})
    @enforce precompile(Tuple{typeof(eltype), Type{result}})
    @enforce precompile(Tuple{typeof(length), result})
    @enforce precompile(Tuple{typeof(parse), Type{result}, String})
    @enforce precompile(Tuple{typeof(parse), Type{result}, SubString{String}})
    @enforce precompile(Tuple{typeof(setindex!), result, Nothing, clustering})

    # input.jl
    @enforce precompile(Tuple{typeof(CrystalNets.parse_cif),String})
    @enforce precompile(Tuple{typeof(CrystalNets.parsestrip), Float64, String})
    @enforce precompile(Tuple{typeof(CrystalNets.parsestrip), Float32, String})
    @enforce precompile(Tuple{typeof(CrystalNets.parsestrip), BigFloat, String})
    @enforce precompile(Tuple{typeof(CrystalNets.parsestrip), String})
    @enforce precompile(Tuple{typeof(CrystalNets.popstring!), Dict{String, Union{Vector{String},String}}, String})
    @enforce precompile(Tuple{typeof(CrystalNets.popvecstring!), Dict{String, Union{Vector{String},String}}, String})
    @enforce precompile(Tuple{Type{CrystalNets.CIF}, Dict{String, Union{Vector{String},String}}})
    @enforce precompile(Tuple{Type{CrystalNets.CIF}, String})
    @enforce precompile(Tuple{typeof(CrystalNets.parse_arc), String})
    @enforce precompile(Tuple{typeof(CrystalNets.parse_arcs), String, Vector{String}})
    @enforce precompile(Tuple{typeof(CrystalNets.parse_arcs), String})
    @enforce precompile(Tuple{typeof(CrystalNets.parse_atom_name), String})
    @enforce precompile(Tuple{typeof(CrystalNets.parse_atom), String})
    @enforce precompile(Tuple{typeof(CrystalNets.chem_atoms), Chemfiles.Residue})
    @enforce precompile(Tuple{typeof(CrystalNets.attribute_residues), Vector{Chemfiles.Residue}, Int, Bool})
    @enforce precompile(Tuple{typeof(CrystalNets.check_collision), p_pos, fmat})
    @enforce precompile(Tuple{typeof(CrystalNets.fix_atom_on_a_bond!), graph, p_pos, fmat})
    @enforce precompile(Tuple{typeof(CrystalNets.least_plausible_neighbours), Vector{Float64}, Int})
    @enforce precompile(Tuple{typeof(CrystalNets.fix_valence!), PeriodicGraph3D, p_pos, types, Vector{Int}, Vector{Int}, Vector{Int}, fmat, Val{true}, opts})
    @enforce precompile(Tuple{typeof(CrystalNets.fix_valence!), PeriodicGraph3D, p_pos, types, Vector{Int}, Vector{Int}, Vector{Int}, fmat, Val{false}, opts})
    @enforce precompile(Tuple{typeof(CrystalNets.sanitize_removeatoms!), graph, p_pos, types, fmat, opts})
    @enforce precompile(Tuple{typeof(CrystalNets.remove_triangles!), graph, p_pos, types, fmat, edgs})
    @enforce precompile(Tuple{typeof(CrystalNets.remove_triangles!), graph, p_pos, types, fmat})
    @enforce precompile(Tuple{typeof(CrystalNets._detect_bent_bond), graph, p_pos, Int, Int, fmat})
    @enforce precompile(Tuple{typeof(CrystalNets.detect_bent_bond), graph, p_pos, Int, Int, fmat})
    @enforce precompile(Tuple{typeof(CrystalNets.sanity_checks!), graph, p_pos, types, fmat, opts})
    @enforce precompile(Tuple{typeof(CrystalNets._remove_homoatomic_bonds!), graph, types, Set{Symbol}, Bool})
    @enforce precompile(Tuple{typeof(CrystalNets._remove_homometallic_bonds!), graph, types, Vector{Int}})
    @enforce precompile(Tuple{typeof(CrystalNets.remove_homoatomic_bonds!), graph, types, Set{Symbol}, Bool})
    @enforce precompile(Tuple{typeof(CrystalNets.finalize_checks), cell, p_pos, types, Vector{Int}, bondlist, Bool, opts, String})
    @enforce precompile(Tuple{typeof(CrystalNets.parse_as_cif), CrystalNets.CIF, opts, String})
    @enforce precompile(Tuple{typeof(CrystalNets.parse_as_chemfile), Chemfiles.Frame, opts, String})
    @enforce precompile(Tuple{typeof(CrystalNets.parse_chemfile), String, opts})
    # precompile_kwarg(Tuple{typeof(parse_chemfile), String}, (kwargs,))

    # @enforce precompile(Tuple{typeof(CrystalNets.), })

    # arithmetics.jl
    for D in 1:3
        for T in inttypes
            @enforce precompile(Tuple{typeof(CrystalNets.find_ratbasis), ppos{D,T}})
            @enforce precompile(Tuple{typeof(CrystalNets.normal_basis_rational), ppos{D,T}})
        end
    end

    # stability.jl
    @enforce precompile(Tuple{Type{collision}, PeriodicGraph{0}, Int, Vector{Int}})
    @enforce precompile(Tuple{typeof(==), collision, collision})
    @enforce precompile(Tuple{typeof(length), collision})
    @enforce precompile(Tuple{typeof(hash), collision, UInt})
    for D in 1:3
        for T in inttypes
            @enforce precompile(Tuple{typeof(CrystalNets.shrink_collisions), cnet{D,T}, Vector{UnitRange{Int}}})
        end
    end
    @enforce precompile(Tuple{typeof(CrystalNets._order_collision), SimpleGraph{Int}, Vector{Vector{Int}}})
    for D in 1:3
        @enforce precompile(Tuple{typeof(CrystalNets.order_collision), PeriodicGraph{D}, UnitRange{Int}})
    end
    @enforce precompile(Tuple{typeof(CrystalNets.collision_utils), Vector{collision}, Int, Nothing})
    @enforce precompile(Tuple{typeof(CrystalNets.collision_utils), Vector{collision}, Int, Vector{Int}})
    @enforce precompile(Tuple{typeof(CrystalNets.collision_utils), Vector{collision}, Int})
    @enforce precompile(Tuple{typeof(CrystalNets.collision_utils), Vector{collision}, Vector{Int}})
    for D in 1:3
        @enforce precompile(Tuple{typeof(CrystalNets.expand_collisions), Vector{collision}, PeriodicGraph{D}, Vector{Int}})
        @enforce precompile(Tuple{typeof(CrystalNets.unsorted_node), PeriodicGraph{D}, UnitRange{Int}})
        @enforce precompile(Tuple{Type{collision}, PeriodicGraph{D}, UnitRange{Int}, Nothing})
        @enforce precompile(Tuple{Type{collision}, PeriodicGraph{D}, UnitRange{Int}})
    end
    @enforce precompile(Tuple{Type{collision}, collision, Vector{Int}})
    for D in 1:3
        for T in inttypes
            @enforce precompile(Tuple{typeof(CrystalNets.collect_collisions), cnet{D,T}})
            @enforce precompile(Tuple{typeof(CrystalNets.collision_nodes), cnet{D,T}})
        end
    end

    # topology.jl
    for D in 1:3
        @enforce precompile(Tuple{typeof(CrystalNets.check_symmetry_with_collisions), Vector{collision}})
        for T in reverse(inttypes)
            @enforce precompile(Tuple{typeof(CrystalNets.check_dimensionality), cnet{D,T}})
            @enforce precompile(Tuple{typeof(CrystalNets.possible_translations), cnet{D,T}})
            @enforce precompile(Tuple{typeof(CrystalNets.find_all_valid_translations), cnet{D,T}, Vector{collision}})
            @enforce precompile(Tuple{typeof(CrystalNets.minimal_volume_matrix), NTuple{D, Vector{Tuple{Int,Int,pos{D,T}}}}})
            @enforce precompile(Tuple{typeof(CrystalNets.reduce_with_matrix), cnet{D,T}, smat{D,T,D*D}, Vector{collision}})
            @enforce precompile(Tuple{typeof(CrystalNets.minimize), cnet{D,T}, Vector{collision}})
            @enforce precompile(Tuple{typeof(CrystalNets.findfirstbasis), ppos{D,T}})
            @enforce precompile(Tuple{typeof(CrystalNets.findbasis), Vector{Tuple{Int,Int,pos{D,T}}}})
            @enforce precompile(Tuple{typeof(CrystalNets.candidate_key), cnet{D,T}, Int, smat{D,T,D*D}, Vector{Tuple{Int,Int,pos{D,T}}}})
        end
        @enforce precompile(Tuple{typeof(CrystalNets.partition_by_coordination_sequence), PeriodicGraph{D}, NoSymmetryGroup})
        for T in inttypes
            @enforce precompile(Tuple{typeof(CrystalNets.partition_by_coordination_sequence), PeriodicGraph{D}, symmgroup{T}})
        end
        @enforce precompile(Tuple{typeof(CrystalNets.partition_by_coordination_sequence), PeriodicGraph{D}})
        for T in reverse(inttypes)
            @enforce precompile(Tuple{typeof(CrystalNets.find_initial_candidates), cnet{D,T}, Vector{Vector{Int}}, Vector{Int}})
            @enforce precompile(Tuple{typeof(CrystalNets.find_candidates_onlyneighbors), cnet{D,T}, Vector{Vector{Int}}, Vector{Int}})
        end
    end
    for T in inttypes
        @enforce precompile(Tuple{typeof(CrystalNets.find_candidates_fallback), cnet{3,T}, Vector{Int}, Vector{Vector{Int}}, Vector{Int}})
        @enforce precompile(Tuple{typeof(CrystalNets.extract_through_symmetry), Dict{Int,Vector{smat{3,T,9}}}, NoSymmetryGroup})
        @enforce precompile(Tuple{typeof(CrystalNets.extract_through_symmetry), Dict{Int,Vector{smat{3,T,9}}}, symmgroup{T}})
        for D in 1:3
            @enforce precompile(Tuple{typeof(CrystalNets.find_candidates), cnet{D,T}, Vector{collision}})
            @enforce precompile(Tuple{typeof(CrystalNets.topological_key), cnet{D,T}, Vector{collision}})
            @enforce precompile(Tuple{typeof(CrystalNets.topological_key), cnet{D,T}})
        end
    end

    # query.jl
    @enforce precompile(Tuple{typeof(CrystalNets.recognize_topology), String, Dict{String,String}})
    @enforce precompile(Tuple{typeof(CrystalNets.recognize_topology), String})
    for T in inttypes
        @enforce precompile(Tuple{typeof(CrystalNets.topological_genome), cnet{0,T}, Vector{collision}})
    end
    for D in 1:3
        @enforce precompile(Tuple{typeof(CrystalNets.recognize_topology), PeriodicGraph{D}, Dict{String,String}})
        @enforce precompile(Tuple{typeof(CrystalNets.recognize_topology), PeriodicGraph{D}})
        for T in reverse(inttypes)
            @enforce precompile(Tuple{typeof(CrystalNets.topological_genome), cnet{D,T}, Vector{collision}})
            @enforce precompile(Tuple{typeof(CrystalNets.topological_genome), cnet{D,T}})
        end
    end
    @enforce precompile(Tuple{typeof(CrystalNets.topological_genome), unets})
    for D in 1:3
        @enforce precompile(Tuple{typeof(CrystalNets.topological_genome), PeriodicGraph{D}, opts})
        @enforce precompile(Tuple{typeof(CrystalNets.topological_genome), PeriodicGraph{D}})
    end
    @enforce precompile(Tuple{typeof(CrystalNets.topological_genome), String, opts})
    @enforce precompile(Tuple{typeof(CrystalNets.topological_genome), String})
    @enforce precompile(Tuple{typeof(CrystalNets._loop_group!), Expr, Symbol, Symbol, Expr})
    @enforce precompile(Tuple{typeof(CrystalNets.determine_topology), String, opts})
    @enforce precompile(Tuple{typeof(CrystalNets.determine_topology), String})
    @enforce precompile(Tuple{typeof(CrystalNets.guess_topology), String, opts})
    @enforce precompile(Tuple{typeof(CrystalNets.guess_topology), String})

    # The following are not precompiled because execution time always dominates compilation time
    # @enforce precompile(Tuple{typeof(CrystalNets.determine_topology_dataset), String, Bool, Bool, opts})
    # @enforce precompile(Tuple{typeof(CrystalNets.determine_topology_dataset), String, Bool, Bool})
    # @enforce precompile(Tuple{typeof(CrystalNets.determine_topology_dataset), String, Bool})
    # @enforce precompile(Tuple{typeof(CrystalNets.determine_topology_dataset), String})
    # @enforce precompile(Tuple{typeof(CrystalNets.guess_topology_dataset), String, Bool, Bool, opts})
    # @enforce precompile(Tuple{typeof(CrystalNets.guess_topology_dataset), String, Bool, Bool})
    # @enforce precompile(Tuple{typeof(CrystalNets.guess_topology_dataset), String, Bool})
    # @enforce precompile(Tuple{typeof(CrystalNets.guess_topology_dataset), String})

    # executable.jl
    @enforce precompile(Tuple{typeof(CrystalNets.parse_commandline), Vector{String}})
    @enforce precompile(Tuple{typeof(CrystalNets.main), Vector{String}})
    @enforce precompile(Tuple{typeof(CrystalNets.main), Vector{SubString{String}}})
    @enforce precompile(Tuple{typeof(CrystalNets.main), String})
    @enforce precompile(Tuple{typeof(CrystalNets.julia_main)})
end

_precompile_()

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
