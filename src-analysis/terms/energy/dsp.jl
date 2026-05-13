@inline function DSP_func(i, j, k, grid, u, v, u_next, v_next, U)

    α = ∂xᶜᵃᵃ(i, j, k, grid, U)

    uu = ℑxᶜᵃᵃ(i, j, k, grid, fGg, u, a_avg, u, u_next)
    vv = ℑyᵃᶜᵃ(i, j, k, grid, fGg, v, a_avg, v, v_next)
    
    return α * (uu - vv)
end

"""
    DSP(u, v, u_next, v_next, U)
Deformation shear production due to a background U
"""
function DSP(u, v, u_next, v_next, U)
    grid = u.grid
    loc = locationornothing((Center, Nothing, Center), u)
    return KernelFunctionOperation{loc...}(DSP_func, grid, u, v, u_next, v_next, U)
end
