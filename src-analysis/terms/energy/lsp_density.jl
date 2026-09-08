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

@doc raw"""
    LSPDensity(velocities, velocities_prev, turbulent_fluxes)
Return a kernel function operation that calculates the production of turbulent kinetic energy from lateral mixing

The lateral shear production is calcuated as
```math
\text{LSP} = \int \,\text{d}x\text{d}z \left [ \overline{u'\vec u'}\cdot {\color{red} \frac{\partial \overline{\vec u}}{\partial x}}\right]
```

This is a component of the mean kinetic energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MKE} = \text{DSP} + \text{WIND} + \text{BUOYANCY} + \text{SPONGE}_\text{MKE} - \text{LSP} - \text{VSP} + \text{STRAIN}_\text{MKE}
```

See also [`MKEDensity`](@ref), [`DSPDensity`](@ref), [`STRESS`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMKEDensity`](@ref), [`VSPDensity`](@ref), [`STRAINMKEDensity`](@ref)
"""
function LSPDensity(velocities, velocities_prev, turbulent_fluxes)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(lsp_density_func, grid, velocities, velocities_prev, turbulent_fluxes)
end
