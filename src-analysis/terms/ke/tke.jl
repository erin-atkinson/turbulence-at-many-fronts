@inline function TKE3D_func(i, j, k, grid, clock, input_fields, coarse_fields)

    u = ℑxᶜᵃᵃ(i, j, k, grid, input_fields.u)
    v = ℑyᵃᶜᵃ(i, j, k, grid, input_fields.v)
    w = ℑzᵃᵃᶜ(i, j, k, grid, input_fields.w)
    
    U = ℑxᶜᵃᵃ(i, j, k, grid, coarse_fields.u_coarse)
    V = ℑyᵃᶜᵃ(i, j, k, grid, coarse_fields.v_coarse)
    W = ℑzᵃᵃᶜ(i, j, k, grid, coarse_fields.w_coarse)

    KE = (u^2 + v^2 + w^2) / 2
    MKE = (U^2 + V^2 + W^2) / 2
    return KE - MKE
end
