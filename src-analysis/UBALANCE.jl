include("terms/terms.jl")

fields = (:u, :v, :w, :p, :uu, :wu, :u_prev)

for ξ in fields
    @eval $ξ = input_fields.$ξ
end

loc = (Face(), Nothing(), Center())

println("Mean flux densities")
flux_density_x = Field(UuFlux(centered, u, u))
flux_density_background = Field(UuFlux(centered, input_fields.U, u))
flux_density_z = Field(WuFlux(centered, w, u))
flux_density = (; flux_density_x, flux_density_background, flux_density_z)

println("Mean advection")
advection_x = Field(@at loc -u * ∂x(u))
advection_background = Field(@at loc -input_fields.U * ∂x(u))
advection_z = Field(@at loc -w * ∂z(u))
advection = (; advection_x, advection_background, advection_z)

println("Turbulent flux densities")
turbulent_flux_density_x = Field(uu - flux_density_x)
turbulent_flux_density_z = Field(wu - flux_density_z)
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_z)

println("Derivatives of turbulent flux densities")
mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_z = Field(-∂z(turbulent_flux_density_z))
mixing = (; mixing_x, mixing_z)

println("Tendency for a fluid parcel")
coriolis = Field(@at loc sp.f * v)
pressure = Field(-∂x(p))
strain = Field(-u * ∂x(input_fields.U))
sponge = Field(SpongeLayer(u, sp))
surface = Field(-SurfaceFluxU(grid, clock, sp))
ageostrophic = Field(pressure + coriolis)
parcel = (; coriolis, pressure, strain, sponge, surface, ageostrophic)

dependency_fields = merge(flux_density, advection, turbulent_flux_density, mixing, parcel)
output_fields = dependency_fields

println("Quadratic balance equation")
balance_terms = (
    :advection_x, :advection_background, :advection_z,
    :mixing_x, :mixing_z,
    :coriolis, :pressure, :strain, :sponge, :surface,
    :ageostrophic
)

u_avg = Field((u + u_prev) / 2)
u_avg_surface = Field(ZSlice(u_avg, 0))
dependency_fields = (; dependency_fields..., u_avg, u_avg_surface)
quadratic = NamedTuple()

for ξ in balance_terms
    quadratic_ξ = Symbol(:quadratic_, ξ)
    field = ξ == :surface ? :u_avg_surface : :u_avg
    @eval begin
        $quadratic_ξ = Field(Integral($ξ * $field))
        quadratic = (; quadratic..., $quadratic_ξ)
    end
end
quadratic_total = Field(
      quadratic_advection_x 
    + quadratic_advection_background
    + quadratic_advection_z
    + quadratic_mixing_x
    + quadratic_mixing_z
    + quadratic_strain
    + quadratic_sponge
    + quadratic_surface
    + quadratic_ageostrophic
)
u² = Field(Integral(u * u))

quadratic = (; quadratic..., quadratic_total, u²)

dependency_fields = merge(dependency_fields, quadratic)
output_fields = merge(output_fields, quadratic)

skip_update = filter(a->a ∉ fields, keys(input_fields))
