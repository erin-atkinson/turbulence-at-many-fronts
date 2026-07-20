include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

fields = (:u, :v, :w, :b, :ub, :wb)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end

loc = (Center(), Nothing(), Center())

flux_density_x = Field(UcFlux(centered, u_bar, b_bar))
flux_density_background = Field(UcFlux(centered, input_fields.U, b_bar))
flux_density_z = Field(WcFlux(centered, w_bar, b_bar))
flux_density = (; flux_density_x, flux_density_background, flux_density_z)

advection_x = Field(@at loc -u_bar * ∂x(b_bar))
advection_background = Field(@at loc -input_fields.U * ∂x(b_bar))
advection_z = Field(@at loc -w_bar * ∂z(b_bar))
advection = (; advection_x, advection_background, advection_z)

turbulent_flux_density_x = Field(ub_bar - flux_density_x - flux_density_background)
turbulent_flux_density_z = Field(wb_bar - flux_density_z)
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_z)

mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_z = Field(-∂z(turbulent_flux_density_z))
mixing = (; mixing_x, mixing_z)

total = Field(-∂x(flux_density_x + flux_density_background) - ∂z(flux_density_z) + mixing_x + mixing_z)

skip_update = filter(a->a ∉ fields, keys(input_fields))
dependency_fields = merge(mean_fields, flux_density, advection, turbulent_flux_density, mixing, (; total))
output_fields = dependency_fields
