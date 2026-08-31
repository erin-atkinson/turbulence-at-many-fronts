# Take a vector of nd arrays and turn it into a single n+1d array
function expand(v::V) where {T, n, A<:Array{T, n}, V<:Vector{A}}
    l = length(v)
    s = size(v[1])
    new_size = (l, s...)
    a = zeros(T, new_size)
    
    inds = [Colon() for x in s] # There must be some syntax I'm missing
    for i in 1:l
        a[i, inds...] .= v[i]
    end
    a
end

expand(v::V) where {T, V<:Vector{T}} = v

@inline function openifstring(func::Function, file) 
    func(file)
end

@inline function openifstring(func::Function, file::String)
    jldopen(func, file)
end

# Helper function for getting all iterations
@inline function timeseries_of(func::Function, file, field, iterations)
    data = openifstring(file) do file
        map(iteration->get_field(func, file, field, iteration), iterations)
    end
    return expand(data)
end

# With no function given
@inline function timeseries_of(file, field, iterations)
    timeseries_of(identity, file, field, iterations)
end

# Over a named tuple of fields
@inline function timeseries_of(func::Function, filename, iterations; kwargs...)
    return NamedTuple(k => timeseries_of(func, filename, kwargs[k], iterations) for k in keys(kwargs))
end

@inline function timeseries_of(filename, iterations; kwargs...)
    return NamedTuple(k => timeseries_of(filename, kwargs[k], iterations) for k in keys(kwargs))
end

