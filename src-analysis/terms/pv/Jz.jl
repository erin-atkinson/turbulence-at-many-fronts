@inline function J_adv_z_func(i, j, k, grid, sp, q, w)
    wᶠⁿᶜ = ℑxzᶠᵃᶜ(i, j, k, grid, w)
    qᶠⁿᶜ = ℑzᵃᵃᶜ(i, j, k, grid, q)

    return wᶠⁿᶜ * qᶠⁿᶜ
end

@inline function J_mix_z_func(i, j, k, grid, sp, v_avg, b_avg, Fy, 𝒟)
    Fyᶠⁿᶜ = ℑxᶠᵃᵃ(i, j, k, grid, Fy)
    𝒟ᶠⁿᶜ = ℑxᶠᵃᵃ(i, j, k, grid, 𝒟)

    vx = ∂xᶠᶜᶜ(i, j, k, grid, v_avg)
    bx = ∂xᶠᶜᶜ(i, j, k, grid, b_avg)
    
    return -Fyᶠⁿᶜ * bx + 𝒟ᶠⁿᶜ * (vx + sp.f)
end

function PV_J_adv_z(sp, q, w)
    grid = q.grid
    loc = (Face, Nothing, Center)
    return KernelFunctionOperation{loc...}(J_adv_z_func, grid, sp, q, w)
end

function PV_J_mix_z(sp, v_avg, b_avg, Fy, 𝒟)
    grid = v_avg.grid
    loc = (Face, Nothing, Center)
    return KernelFunctionOperation{loc...}(J_mix_z_func, grid, sp, v_avg, b_avg, Fy, 𝒟)
end
