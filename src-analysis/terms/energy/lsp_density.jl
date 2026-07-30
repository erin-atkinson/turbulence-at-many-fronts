@inline function lsp_density_func(i, j, k, grid, velocities, velocities_prev, turbulent_fluxes)

    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_prev = velocities_prev.u
    v_prev = velocities_prev.v
    w_prev = velocities_prev.w

    uu = turbulent_fluxes.uu
    uv = turbulent_fluxes.uv
    uw = turbulent_fluxes.uw

    uuux = ∂xᶜᶜᶜ(i, j, k, grid, fGg, uu, a_avg, u, u_prev)
    uvvx = ℑxyᶜᶜᵃ(i, j, k, grid, fGg, uv, ∂xᶠᶠᶜ, a_avg, v, v_prev)
    uwwx = ℑxzᶜᵃᶜ(i, j, k, grid, fGg, uw, ∂xᶠᶜᶠ, a_avg, w, w_prev)
    
    return -(
          uuux
        + uvvx
        + uwwx
    )
end

"""
    LSPDensity(velocities, velocities_prev, turbulent_fluxes)
Input coarse-grained velocities and fluxes for lateral shear production density
"""
function LSPDensity(velocities, velocities_prev, turbulent_fluxes)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(lsp_density_func, grid, velocities, velocities_prev, turbulent_fluxes)
end
