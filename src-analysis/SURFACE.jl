# Terms in the surface buoyancy balance equation
include("terms/terms.jl")

fields = (:u, :w, :b, :ub, :wb, :b_prev)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = afm(input_fields.$ξ)
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end

# Some definition of the surface layer
ε = sp.Δb/6
mld = Field(MLD(b_bar, ε))

constant_mld = Field(ConstantMLD(b_bar, ε))
constant_mld_prev = Field(ConstantMLD(b_prev_bar, ε))

growth_rate = Field(GrowthRate(clock, constant_mld, constant_mld_prev))
mld_fields = (; mld, constant_mld, constant_mld_prev, growth_rate)

for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    ξ_mld = Symbol(ξ, :_mld)
    @eval begin
        $ξ_bar = Field(ML_Average($ξ_bar, constant_mld))
        mld_fields = (; mld_fields..., $ξ_mld)
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

turbulent_flux_density_z = Field(wb_bar - flux_density_z)
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_z)

mixing_z = Field(-∂z(turbulent_flux_density_z))
mixing = (; mixing_x, mixing_z)

mld_flux_density_x = Field(UcFlux(centered, u_mld, b_mld))
mld_advection_x = Field(@at (loc[1], loc[2], Nothing()) -u_mld * ∂x(b_mld))
mld_advection_background = Field(ML_Average(flux_density_background, constant_mld))

shear_dispersion_flux = Field(flux_density - mld_flux_density)
mld_turbulent_flux = Field(ML_Average(turbulent_flux_density_x, constant_mld))

shear_dispersion = Field(-∂x(shear_dispersion_flux))
mixing = Field(-∂x(mld_turbulent_flux))


