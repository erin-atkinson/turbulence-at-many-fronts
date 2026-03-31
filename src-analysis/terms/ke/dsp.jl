@inline function DSP_func(i, j, k, grid, u, v, u_next, v_next, U)

    α = ∂xᶜᵃᵃ(i, j, k, grid, U)
    
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u)
    vᶜᶜᶜ = ℑyᵃᶜᵃ(i, j, k, grid, v)
    
    uᶜᶜᶜ_next = ℑxᶜᵃᵃ(i, j, k, grid, u_next)
    vᶜᶜᶜ_next = ℑyᵃᶜᵃ(i, j, k, grid, v_next)

    uu_dfm = ℑxᶜᵃᵃ(i, j, k, grid, fGg, u_dfm, a_avg, u_dfm, u_next_dfm)
    vv_dfm = ℑyᵃᶜᵃ(i, j, k, grid, fGg, v_dfm, a_avg, v_dfm, v_next_dfm)
    
    return α * (uᶜᶜᶜ_next * uᶜᶜᶜ - vᶜᶜᶜ_next * vᶜᶜᶜ)
end

"""
    DSP(u, v, u_next, v_next, U)
Deformation shear production due to a background U
"""
function DSP(u, v, u_next, v_next, U)
    grid = u.grid
    loc = location(u)
    return KernelFunctionOperation{Center, Nothing, Center}(DSP_func, grid, u, v, u_next, v_next, U)
end
