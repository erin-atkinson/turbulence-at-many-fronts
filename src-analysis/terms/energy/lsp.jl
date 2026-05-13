@inline function LSP_func(i, j, k, grid, velocities, velocities_next, fluxes)

    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_next = velocities_next.u
    v_next = velocities_next.v
    w_next = velocities_next.w

    uu = fluxes.uu
    uv = fluxes.uv
    uw = fluxes.uw

    uuux = ∂xᶜᶜᶜ(i, j, k, grid, fGg, uu, a_avg, u, u_next)
    uvvx = ℑxyᶜᶜᵃ(i, j, k, grid, fGg, uv, ∂xᶠᶠᶜ, a_avg, v, v_next)
    uwwx = ℑxzᶜᵃᶜ(i, j, k, grid, fGg, uw, ∂xᶠᶜᶠ, a_avg, w, w_next)
    
    return -(
          uuux
        + uvvx
        + uwwx
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
