@inline function J_adv_x_func(i, j, k, grid, sp, q, u, U)
    u_tot = ℑxzᶜᵃᶠ(i, j, k, grid, u) + ℑxzᶜᵃᶠ(i, j, k, grid, U)
    qᶜⁿᶠ = ℑxᶜᵃᵃ(i, j, k, grid, q)
    
    return u_tot * qᶜⁿᶠ
end

@inline function J_mix_x_func(i, j, k, grid, sp, v_avg, b_avg, Fy, 𝒟)
    Fyᶜⁿᶠ = ℑzᵃᵃᶠ(i, j, k, grid, Fy)
    𝒟ᶜⁿᶠ = ℑzᵃᵃᶠ(i, j, k, grid, 𝒟)

    vz = ∂zᶜᶜᶠ(i, j, k, grid, v_avg)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b_avg)
    
    return Fyᶜⁿᶠ * bz - 𝒟ᶜⁿᶠ * vz
end

function PV_J_adv_x(sp, q, u, U)
    grid = q.grid
    loc = (Center, Nothing, Face)
    return KernelFunctionOperation{loc...}(J_adv_x_func, grid, sp, q, u, U)
end

function PV_J_mix_x(sp, v_avg, b_avg, Fy, 𝒟)
    grid = v_avg.grid
    loc = (Center, Nothing, Face)
    return KernelFunctionOperation{loc...}(J_mix_x_func, grid, sp, v_avg, b_avg, Fy, 𝒟)
end