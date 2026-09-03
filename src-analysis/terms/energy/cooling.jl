function cooling_func(i, j, k, grid, clock, h_ml, h_ml_prev, sp)
    t = clock.time
    B = b_flux_func(t, sp)
    return @inbounds a_avg(i, j, k, grid, h_ml, h_ml_prev) * B * sp.Lx 
end

@doc raw"""
    Cooling(clock, h_ml, h_ml_prev, sp)
Potential energy input by the surface cooling, referenced to the base of the mixed layer
```math
\text{COOLING} = Bh
```
"""
function COOLING(clock, h_ml, h_ml_prev, sp)
    grid = h_ml.grid
    return KernelFunctionOperation{Nothing, Nothing, Nothing}(cooling_func, grid, clock, h_ml, h_ml_prev, sp)
end
