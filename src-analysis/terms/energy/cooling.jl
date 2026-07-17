function cooling_func(i, j, k, grid, clock, h_ml, h_ml_next, sp)
    t = clock.time
    B = b_flux_func(t, sp)
    return @inbounds -a_avg(i, j, k, grid, h_ml, h_ml_next) * B * sp.Lx 
end

function COOLING(clock, h_ml, h_ml_next, sp)
    grid = h_ml.grid
    return KernelFunctionOperation{Nothing, Nothing, Nothing}(cooling_func, grid, clock, h_ml, h_ml_next, sp)
end
