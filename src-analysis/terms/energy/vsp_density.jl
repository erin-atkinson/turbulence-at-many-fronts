@inline function vsp_density_func(i, j, k, grid, velocities, velocities_prev, turbulent_fluxes)

    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_prev = velocities_prev.u
    v_prev = velocities_prev.v
    w_prev = velocities_prev.w

    wu = turbulent_fluxes.wu
    wv = turbulent_fluxes.wv
    ww = turbulent_fluxes.ww

    wuuz = ℑxzᶜᵃᶜ(i, j, k, grid, fGg, wu, ∂zᶠᶜᶠ, a_avg, u, u_prev)
    wvvz = ℑyzᵃᶜᶜ(i, j, k, grid, fGg, wv, ∂zᶜᶠᶠ, a_avg, v, v_prev)
    wwwz = ∂zᶜᶜᶜ(i, j, k, grid, fGg, ww, a_avg, w, w_prev)
    
    return -(
          wuuz
        + wvvz
        + wwwz
    )
end

@doc raw"""
    VSPDensity(velocities, velocities_prev, turbulent_fluxes)
Return a kernel function operation that calculates the production of turbulent kinetic energy from vertical mixing

The vertical shear production is calcuated as
```math
\text{VSP} = \int \,\text{d}x\text{d}z \left [ \overline{w'\vec u'}\cdot {\color{red} \frac{\partial \overline{\vec u}}{\partial z}}\right]
```

This is a component of the mean kinetic energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MKE} = \text{DSP} + \text{WIND} + \text{BUOYANCY} + \text{SPONGE}_\text{MKE} - \text{LSP} - \text{VSP} + \text{STRAIN}_\text{MKE}
```

See also [`MKEDensity`](@ref), [`DSPDensity`](@ref), [`STRESS`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMKEDensity`](@ref), [`LSPDensity`](@ref), [`STRAINMKEDensity`](@ref)
"""
function VSPDensity(velocities, velocities_prev, turbulent_fluxes)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(vsp_density_func, grid, velocities, velocities_prev, turbulent_fluxes)
end
