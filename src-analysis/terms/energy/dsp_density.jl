@inline function dsp_density_func(i, j, k, grid, u, v, u_next, v_next, U)

    α = -∂xᶜᵃᵃ(i, j, k, grid, U)

    uu = ℑxᶜᵃᵃ(i, j, k, grid, fGg, u, a_avg, u, u_next)
    vv = ℑyᵃᶜᵃ(i, j, k, grid, fGg, v, a_avg, v, v_next)
    
    return α * (uu - vv)
end

"""
    DSPDensity(u, v, u_next, v_next, U)
Deformation shear production density due to a background U
"""
function DSPDensity(u, v, u_next, v_next, U)
    grid = u.grid
    loc = locationornothing((Center, Center, Center), u)
    return KernelFunctionOperation{loc...}(dsp_density_func, grid, u, v, u_next, v_next, U)
end
