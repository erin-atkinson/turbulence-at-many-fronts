include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

fields = (:u, :v, :w, :p, :uu, :wu)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end

loc = (Face(), Nothing(), Center())

flux_density_x = Field(UuFlux(centered, u_bar, u_bar))
flux_density_background = Field(UuFlux(centered, input_fields.U, u_bar))
flux_density_z = Field(WuFlux(centered, w_bar, u_bar))
flux_density = (; flux_density_x, flux_density_background, flux_density_z)

advection_x = Field(@at loc -u_bar * ∂x(u_bar))
advection_background = Field(@at loc -input_fields.U * ∂x(u_bar))
advection_z = Field(@at loc -w_bar * ∂z(u_bar))
advection = (; advection_x, advection_background, advection_z)

turbulent_flux_density_x = Field(uu_bar - flux_density_x)
turbulent_flux_density_z = Field(wu_bar - flux_density_z)
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_z)

mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_z = Field(-∂z(turbulent_flux_density_z))
mixing = (; mixing_x, mixing_z)

coriolis_x = Field(@at loc sp.f * v_bar)
pressure_x = Field(-∂x(p_bar))
strain_x = Field(-u_bar * ∂x(input_fields.U))
sponge = Field(SpongeLayer(u_bar, sp))
other = (; coriolis_x, pressure_x, strain_x, sponge)

total = Field(-∂x(flux_density_x) - ∂z(flux_density_z) + mixing_x + mixing_z + advection_background + coriolis_x + pressure_x + strain_x)

skip_update = filter(a->a ∉ fields, keys(input_fields))
dependency_fields = merge(mean_fields, flux_density, advection, turbulent_flux_density, mixing, other, (; total))
output_fields = dependency_fields
