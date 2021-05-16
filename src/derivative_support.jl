using ChainRulesCore, Zygote

#=
# forward differentiation rule:
function ChainRulesCore.frule(
    (_, _, Δgen),
    ::typeof(IndexFunArray{T,N,F}),
    gen::F,
    sz::NTuple{N, Int},
    ) where {N,F}
    Ω = IndexFunArray{T,N,typeof(gen)}(gen, sz)
    ∂Ω = IndexFunArray{T,N,typeof(gen)}(Δgen, sz)
    return (Ω, ∂Ω)
end
=#


# this function provides support for derivatives in reverse mode as calculated 
# by automatic differentiation tools such as Zygote
function ChainRulesCore.rrule(::Type{IndexFunArray}, # ::typeof(IndexFunArray) 
    myT::Type{T},
    gen::F,
    sz::NTuple{N, Int},  # This is the standard constructor, abeit not the direct contructor
    ) where {T,F,N}
    @show "Hi! In Chair Rule!"
    val_grad(idx) = Zygote._pullback(gen, idx)[2](1.0) # [2](1.0)
    # @show  val_grad([1,1])[1]  # this is a named tuple
    # gradgen(idx) = val_grad(idx)[1][:a] # (val_grad(x)[2](1.0))[2] # mygrad(x)[1][1] # 
    @show mySymbols=keys(val_grad(sz)[1])
    @show gradgen(idx) = val_grad(idx)[1] # (val_grad(x)[2](1.0))[2] # mygrad(x)[1][1] # 
    function IFA_pullback(ΔΩ)   # Zygote._pullback(f,pi/2.0)[2](1.0) # 1.0 is only the seed
        Fcts = ((idx)-> val_grad(idx)[1][aSymbol] for aSymbol in mySymbols)
        print("oja:\n\n")
        TupleVals = (apply_tuple_list.(.*,ΔΩ,IndexFunArray(typeof(Fun(sz)), Fun, sz)) for Fun in Fcts) # {T,N,typeof(Fun)} yields a tuple of values
        # TupleVals = (collect(ΔΩ .* IndexFunArray(typeof(Fun(sz)), Fun, sz)) for Fun in Fcts) # {T,N,typeof(Fun)} yields a tuple of values
        @show TupleVals
        print("ojo:\n\n")
        ∂gen = NamedTuple{mySymbols}(TupleVals) # converts the value tuple into a named tuple with the appropriate type names
        @show typeof(∂gen)
        @show ∂gen
        print("oji:\n\n")
        return (NO_FIELDS, NO_FIELDS, ∂gen, NO_FIELDS) # myT, gen, sz, closures   NO_FIELDS  (offset=∂gen[:offset],)
    end
    Ω = IndexFunArray(T, gen, sz) # {T,N,typeof(gen)}
    return (Ω, IFA_pullback)
end
