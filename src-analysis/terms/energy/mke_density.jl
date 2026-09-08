@inline function mke_density_func(i, j, k, grid, velocities)
    u = velocities.u
    v = velocities.v
    w = velocities.w

    uu = ℑxᶜᵃᵃ(i, j, k, grid, fg, u, u)
    vv = ℑyᵃᶜᵃ(i, j, k, grid, fg, v, v)
    ww = ℑzᵃᵃᶜ(i, j, k, grid, fg, w, w)

    return (uu + vv + ww) / 2
end

@doc raw"""
    MKEDensity(velocities)
Return a kernel function operation that calculates the mean kinetic energy density

The mean kinetic energy is calculated as
```math
\text{MKE} = \int\text{d}x\text{d}z\, \frac{1}{2} (\overline u^2 + \overline v^2 + \overline w^2)
```

This is a component of the mean kinetic energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MKE} = \text{DSP} + \text{WIND} + \text{BUOYANCY} + \text{SPONGE}_\text{MKE} - \text{LSP} - \text{VSP} + \text{STRAIN}_\text{MKE}
```

See also [`DSPDensity`](@ref), [`STRESS`](@ref), [`BUOYANCYDensity`](@ref), [`SPONGEMKEDensity`](@ref), [`LSPDensity`](@ref), [`VSPDensity`](@ref), [`STRAINMKEDensity`](@ref)
"""
function MKEDensity(velocities)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(mke_density_func, grid, velocities)
end
