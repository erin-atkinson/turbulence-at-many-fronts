flux_mixed_layer_fit(a, z) = max(z / a[1] + one(z), zero(z))

@inline function flux_based_mixed_layer_func(i, j, k, grid, w′b′, zs, sp, scalar_model, u0)
    fluxes = @inbounds w′b′[i, j, 1:length(zs)] / sp.B
    
    # Solve problem
    prob = NonlinearCurveFitProblem(scalar_model, u0, zs, fluxes)
    sol = solve(prob)
    h = sol.u[1]
    
    return max(h, zero(h))
end

@doc raw"""
    FluxBasedMixedLayer(w′b′, sp)
Return a kernel function operation that calculates the mixed layer depth using the turbulent flux

See also [`MPEDensity`](@ref), [`COOLING`](@ref), [`MIXEDDensity`](@ref)
"""
function FluxBasedMixedLayer(w′b′, sp)
    grid = w′b′.grid
    loc = locationornothing((Center, Center, Nothing), w′b′)
    zs = znodes(w′b′)[1:end-1]
    
    scalar_model = ScalarModel(flux_mixed_layer_fit)
    u0 = [sp.H]
    
    return KernelFunctionOperation{loc...}(flux_based_mixed_layer_func, grid, w′b′, zs, sp, scalar_model, u0)
end
