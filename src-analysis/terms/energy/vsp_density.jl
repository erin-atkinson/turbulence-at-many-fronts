@inline function vsp_density_func(i, j, k, grid, velocities, velocities_next, turbulent_fluxes)

    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_next = velocities_next.u
    v_next = velocities_next.v
    w_next = velocities_next.w

    wu = turbulent_fluxes.wu
    wv = turbulent_fluxes.wv
    ww = turbulent_fluxes.ww

    wuuz = ℑxzᶜᵃᶜ(i, j, k, grid, fGg, wu, ∂zᶠᶜᶠ, a_avg, u, u_next)
    wvvz = ℑyzᵃᶜᶜ(i, j, k, grid, fGg, wv, ∂zᶜᶠᶠ, a_avg, v, v_next)
    wwwz = ∂zᶜᶜᶜ(i, j, k, grid, fGg, ww, a_avg, w, w_next)
    
    return -(
          wuuz
        + wvvz
        + wwwz
    )
end

"""
    VSPDensity(velocities, velocities_next, fluxes)
Input coarse-grained velocities and fluxes for vertical shear production density
"""
function VSPDensity(velocities, velocities_next, turbulent_fluxes)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(vsp_density_func, grid, velocities, velocities_next, turbulent_fluxes)
end
