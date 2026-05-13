@inline function VSP_func(i, j, k, grid, velocities, velocities_next, fluxes)

    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_next = velocities_next.u
    v_next = velocities_next.v
    w_next = velocities_next.w

    wu = fluxes.wu
    wv = fluxes.wv
    ww = fluxes.ww

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
    VSP(velocities, velocities_next, fluxes)
Input coarse-grained velocities and fluxes for vertical shear production
"""
function VSP(velocities, velocities_next, fluxes)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(VSP_func, grid, velocities, velocities_next, fluxes)
end
