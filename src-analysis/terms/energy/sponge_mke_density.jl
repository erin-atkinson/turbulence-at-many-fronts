@inline function sponge_mke_density_func(i, j, k, grid, velocities, velocities_prev, sp)
    u = velocities.u
    v = velocities.v
    w = velocities.w
    
    u_prev = velocities_prev.u
    v_prev = velocities_prev.v
    w_prev = velocities_prev.w

    σuu = ℑxᶜᵃᵃ(i, j, k, grid, σff_avg, (Face(), Center(), Center()), u, u_prev, sp)
    σvv = ℑyᵃᶜᵃ(i, j, k, grid, σff_avg, (Center(), Face(), Center()), v, v_prev, sp)
    σww = ℑzᵃᵃᶜ(i, j, k, grid, σff_avg, (Center(), Center(), Face()), w, w_prev, sp)

    return σuu + σvv + σww
end

@doc raw"""
    SPONGEMKEDensity(velocities, velocities_prev, sp)
Return a kernel function operation that calculates the work done by the sponge layer on the resolved flow

This is calcuated as
```math
\text{SPONGE}_\text{MKE} = \int \,\text{d}x\text{d}z \left [ -\sigma \,\overline{\vec u}\cdot {\color{red} \overline{\vec u}}\right]
```

This is a component of the mean kinetic energy equation:
```math
\frac{\text{d}}{\text{d}t}\text{MKE} = \text{DSP} + \text{WIND} + \text{BUOYANCY} + \text{SPONGE}_\text{MKE} - \text{LSP} - \text{VSP} + \text{STRAIN}_\text{MKE}
```

See also [`MKEDensity`](@ref), [`DSPDensity`](@ref), [`STRESS`](@ref), [`BUOYANCYDensity`](@ref), [`LSPDensity`](@ref), [`VSPDensity`](@ref), [`STRAINMKEDensity`](@ref)
"""
function SPONGEMKEDensity(velocities, velocities_prev, sp)
    grid = velocities.u.grid
    loc = locationornothing((Center, Center, Center), velocities.u)
    return KernelFunctionOperation{loc...}(sponge_mke_density_func, grid, velocities, velocities_prev, sp)
end