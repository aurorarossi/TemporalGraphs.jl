# Binary min-heaps from DataStructures.jl. Elements are tuples `(key, id)`, so ties
# are broken by id and the settling order is deterministic.
const MinHeap{T} = BinaryMinHeap{T}

# Per-node vectors that are allocated on first use and, between two computations,
# reset only where they were used: a single-source computation costs O(reached)
# instead of O(n) allocations. `used[v]` (one byte per node, read on every edge of
# the scans) tells whether the list of v is in use in the current computation.
struct LazyLists{L}
    lists::Vector{Vector{L}}
    used::Vector{Bool}
    touched::Vector{Int}
end
LazyLists{L}(n::Integer) where {L} = LazyLists{L}(Vector{Vector{L}}(undef, n), zeros(Bool, n), Int[])

# the list of v, created if needed (the caller inserts into it)
@inline function _list!(z::LazyLists{L}, v::Int) where {L}
    @inbounds if !z.used[v]
        z.used[v] = true
        push!(z.touched, v)
        isassigned(z.lists, v) || (z.lists[v] = L[])
    end
    return @inbounds z.lists[v]
end

@inline _has_list(z::LazyLists, v::Int) = @inbounds z.used[v]

function _reset!(z::LazyLists)
    @inbounds for v in z.touched
        empty!(z.lists[v])
        z.used[v] = false
    end
    empty!(z.touched)
    return z
end

# Run f(), rethrowing the original exception instead of the TaskFailedException of a
# failed parallel task.
function _unwrap_task_errors(f)
    try
        return f()
    catch err
        while true
            if err isa TaskFailedException
                err = err.task.exception
            elseif err isa CompositeException && !isempty(err.exceptions)
                err = first(err.exceptions)
            else
                break
            end
        end
        throw(err)
    end
end
