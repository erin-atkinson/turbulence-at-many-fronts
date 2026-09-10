include("terms/terms.jl")

fields = (:u, :v, :w, :b, :ub, :wb, :b_prev)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end

loc = (Center(), Nothing(), Center())

println("Mean flux densities")
flux_density_x = Field(UcFlux(centered, u_bar, b_bar))
flux_density_background = Field(UcFlux(centered, input_fields.U, b_bar))
flux_density_z = Field(WcFlux(centered, w_bar, b_bar))
flux_density = (; flux_density_x, flux_density_background, flux_density_z)

println("Mean advection")
advection_x = Field(@at loc -u_bar * ∂x(b_bar))
advection_background = Field(@at loc -input_fields.U * ∂x(b_bar))
advection_z = Field(@at loc -w_bar * ∂z(b_bar))
advection = (; advection_x, advection_background, advection_z)

println("Turbulent flux densities")
turbulent_flux_density_x = Field(ub_bar - flux_density_x)
turbulent_flux_density_z = Field(wb_bar - flux_density_z)
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_z)

println("Derivatives of turbulent flux densities")
mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_z = Field(-∂z(turbulent_flux_density_z))
mixing = (; mixing_x, mixing_z)

println("Tendency for a fluid parcel")
surface = Field(-SurfaceFluxB(grid, clock, sp))
parcel = (; surface)

dependency_fields = merge(flux_density, advection, turbulent_flux_density, mixing, parcel)
output_fields = dependency_fields

println("Quadratic balance equation")
balance_terms = (
    :advection_x, :advection_background, :advection_z,
    :mixing_x, :mixing_z,
    :surface
)

b_avg = Field((b_bar + b_prev_bar) / 2)
b_avg_surface = Field(ZSlice(b_avg, 0))
dependency_fields = (; dependency_fields..., b_avg, b_avg_surface)
quadratic = NamedTuple()

for ξ in balance_terms
    quadratic_ξ = Symbol(:quadratic_, ξ)
    field = ξ == :surface ? :b_avg_surface : :b_avg
    @eval begin
        $quadratic_ξ = Field(Integral($ξ * $field))
        quadratic = (; quadratic..., $quadratic_ξ)
    end
end
quadratic_total = Field(sum(quadratic))
quadratic = (; quadratic..., quadratic_total)

dependency_fields = merge(dependency_fields, quadratic)
output_fields = merge(output_fields, quadratic)

skip_update = filter(a->a ∉ fields, keys(input_fields))

