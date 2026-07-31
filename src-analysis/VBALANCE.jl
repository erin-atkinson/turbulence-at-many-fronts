include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

fields = (:u, :v, :w, :uv, :wv)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end

loc = (Center(), Nothing(), Center())

flux_density_x = Field(UvFlux(centered, u_bar, v_bar))
flux_density_background = Field(UvFlux(centered, input_fields.U, v_bar))
flux_density_z = Field(WvFlux(centered, w_bar, v_bar))
flux_density = (; flux_density_x, flux_density_background, flux_density_z)

advection_x = Field(@at loc -u_bar * ∂x(v_bar))
advection_background = Field(@at loc -input_fields.U * ∂x(v_bar))
advection_z = Field(@at loc -w_bar * ∂z(v_bar))
advection = (; advection_x, advection_background, advection_z)

turbulent_flux_density_x = Field(uv_bar - flux_density_x)
turbulent_flux_density_z = Field(wv_bar - flux_density_z)
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_z)

mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_z = Field(-∂z(turbulent_flux_density_z))
mixing = (; mixing_x, mixing_z)

coriolis_y = Field(@at loc -1 * sp.f * u_bar)
strain_y = Field(v_bar * ∂x(input_fields.U))
sponge = Field(SpongeLayer(v_bar, sp))

other = (; coriolis_y, strain_y, sponge)

total = Field(-∂x(flux_density_x) - ∂z(flux_density_z) + mixing_x + mixing_z + advection_background + coriolis_y + strain_y + sponge)

skip_update = filter(a->a ∉ fields, keys(input_fields))
dependency_fields = merge(mean_fields, flux_density, advection, turbulent_flux_density, mixing, other, (; total))
output_fields = dependency_fields
