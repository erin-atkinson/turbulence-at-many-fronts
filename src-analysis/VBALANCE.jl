include("terms/terms.jl")

fields = (:u, :v, :w, :uv, :wv, :v_prev)

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
flux_density_x = Field(UvFlux(centered, u_bar, v_bar))
flux_density_background = Field(UvFlux(centered, input_fields.U, v_bar))
flux_density_z = Field(WvFlux(centered, w_bar, v_bar))
flux_density = (; flux_density_x, flux_density_background, flux_density_z)

println("Mean advection")
advection_x = Field(@at loc -u_bar * ∂x(v_bar))
advection_background = Field(@at loc -input_fields.U * ∂x(v_bar))
advection_z = Field(@at loc -w_bar * ∂z(v_bar))
advection = (; advection_x, advection_background, advection_z)

println("Turbulent flux densities")
turbulent_flux_density_x = Field(uv_bar - flux_density_x)
turbulent_flux_density_z = Field(wv_bar - flux_density_z)
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_z)

println("Derivatives of turbulent flux densities")
mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_z = Field(-∂z(turbulent_flux_density_z))
mixing = (; mixing_x, mixing_z)

println("Tendency for a fluid parcel")
coriolis = Field(@at loc -1 * sp.f * u_bar)
strain = Field(v_bar * ∂x(input_fields.U))
sponge = Field(SpongeLayer(v_bar, sp))
surface = Field(-SurfaceFluxV(grid, clock, sp))
parcel = (; coriolis, strain, sponge, surface)

dependency_fields = merge(flux_density, advection, turbulent_flux_density, mixing, parcel)
output_fields = dependency_fields

println("Quadratic balance equation")
balance_terms = (
    :advection_x, :advection_background, :advection_z,
    :mixing_x, :mixing_z,
    :coriolis, :strain, :sponge, :surface
)

v_avg = Field((v_bar + v_prev_bar) / 2)
v_avg_surface = Field(ZSlice(v_avg, 0))
dependency_fields = (; dependency_fields..., v_avg, v_avg_surface)
quadratic = NamedTuple()

for ξ in balance_terms
    quadratic_ξ = Symbol(:quadratic_, ξ)
    field = ξ == :surface ? :v_avg_surface : :v_avg
    @eval begin
        $quadratic_ξ = Field(Integral($ξ * $field))
        quadratic = (; quadratic..., $quadratic_ξ)
    end
end
quadratic_total = Field(sum(quadratic))
v² = Field(Integral(v * v))

quadratic = (; quadratic..., quadratic_total, v²)

dependency_fields = merge(dependency_fields, quadratic)
output_fields = merge(output_fields, quadratic)

skip_update = filter(a->a ∉ fields, keys(input_fields))
