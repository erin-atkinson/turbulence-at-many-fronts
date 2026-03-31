@inline function LSP_func(i, j, k, grid, velocities, velocities_next, fluxes)

    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_next = velocities_next.u
    v_next = velocities_next.v
    w_next = velocities_next.w

    uu = @inbounds fluxes.uu[i, j, k]
    uv = @inbounds fluxes.uv[i, j, k]
    uw = @inbounds fluxes.uw[i, j, k]

    ux = ∂xᶜᶜᶜ(i, j, k, grid, a_avg, u, u_next)
    vx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, a_avg, v, v_next)
    wx = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, a_avg, w, w_next)
    
    return -(
          uu * ux
        + uv * vx
        + uw * wx
    )
end

"""
    LSP(velocities, velocities_next, fluxes)
Input coarse-grained velocities and fluxes for lateral shear production
"""
function LSP(velocities, velocities_next, fluxes)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(LSP_func, grid, velocities, velocities_next, fluxes)
end
