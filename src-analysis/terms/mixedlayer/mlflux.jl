# Entrainment across mixed layer
@inline function ml_flux_func(i, j, k, grid, clock, ℓx, ℓy, ℓz, mld, mld_next, field, field_mld)
    field_at_ml = ml_interpolate_func(i, j, k, grid, ℓx, ℓy, ℓz, mld, field)
    dhdt = @inbounds (mld_next[i, j, k] - mld[i, j, k]) / clock.last_Δt

    
    return @inbounds dhdt / mld[i, j, k] * (field_at_ml - field_mld[i, j, k])
end

function ML_Flux(clock, field, field_mld, mld, mld_next)
    (ℓx, ℓy, ℓz) = location(field)
    grid = field.grid
    
    return KernelFunctionOperation{ℓx, ℓy, Nothing}(ml_flux_func, grid, clock, ℓx(), ℓy(), ℓz(), mld, mld_next, field, field_mld)
end