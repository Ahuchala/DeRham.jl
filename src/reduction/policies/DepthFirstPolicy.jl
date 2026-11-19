struct DepthFirstPolicy <: AbstractReductionPolicy
    degree::Int
    S::Vector{Int}
end

function _depthfirst_chunk_length(u, direction, maximum)
    chunk = maximum
    for i = 1:length(u)
        if direction[i] > 0
            chunk = min(chunk, fld(u[i], direction[i]))
        end
    end
    return chunk
end

function run_reduction_policy!(
    policy::DepthFirstPolicy,
    engine::AbstractReductionEngine,
    keys_by_order,
    directions::ReductionDirectionSet,
)
    terminal_order = reduction_engine_geometry(engine).terminal_order
    terminal_chains = Any[]

    for (_, inputs) in keys_by_order
        for input in inputs
            chain = materialize(engine, ReductionChainKey(input.key[1], input.key[2]))
            while chain.key.order > terminal_order
                direction = rev_chooseV(collect(chain.key.u), policy.degree, policy.S)
                maximum_chunk = chain.key.order - terminal_order
                chunk = _depthfirst_chunk_length(chain.key.u, direction, maximum_chunk)
                chunk >= 1 || throw(ArgumentError("depthfirst selected an empty chunk"))
                directions(chain.key.u, chunk, direction) || throw(
                    ArgumentError("depthfirst selected an invalid reduction direction"),
                )
                reduce!(engine, chain, direction, chunk)
            end
            push!(terminal_chains, finish!(engine, chain))
        end
    end

    return terminal_chains
end

public DepthFirstPolicy
