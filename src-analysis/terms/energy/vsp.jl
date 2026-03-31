@inline function VSP_func(i, j, k, grid, velocities, velocities_next, fluxes)

    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_next = velocities_next.u
    v_next = velocities_next.v
    w_next = velocities_next.w

    uu = fluxes.uu
    uv = fluxes.uv
    uw = fluxes.uw

    uz = ∂xzᶜᶜᶜ(i, j, k, grid, ∂zᶠᶜᶠ, a_avg, u, u_next)
    vz = ℑxᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, a_avg, v, v_next)
    wz = ∂xᶜᶜᶜ(i, j, k, grid, a_avg, w, w_next)
    
    return -(
          uu * uz
        + uv * vz
        + uw * wz
    )
end

"""
    VSP(velocities, velocities_next, fluxes)
Input coarse-grained velocities and fluxes for vertical shear production
"""
function VSP(velocities, velocities_next, fluxes)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(VSP_func, grid, velocities, velocities_next, fluxes)
end
