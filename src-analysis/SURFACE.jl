# Terms in the surface buoyancy balance equation
include("terms/terms.jl")
include("terms/CoarseGraining.jl")
include("terms/advection/advection.jl")
include("terms/advection/operators.jl")
include("terms/mixedlayer/mld.jl")

fields = (:u, :v, :w, :b, :ub, :wb, :b_prev)
kernel = Gaussian(grid, 50, 0, 0)

mean_fields = NamedTuple()
for ξ in fields
    ξ_afm = Symbol(ξ, :_afm)
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_afm = afm(input_fields.$ξ)
        $ξ_bar = Field(Coarse($ξ_afm, kernel))
        mean_fields = (; mean_fields..., $ξ_afm, $ξ_bar)
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
        $ξ_mld = Field(ML_Average($ξ_bar, constant_mld))
        mld_fields = (; mld_fields..., $ξ_mld)
    end
end

loc = (Center(), Nothing(), Center())

flux_density_x = Field(UcFlux(centered, u_bar, b_bar))
flux_density_background = Field(UcFlux(centered, input_fields.U, b_bar))
flux_density_z = Field(WcFlux(centered, u_bar, b_bar))
flux_density = (; flux_density_x, flux_density_background)

advection_x = Field(@at loc -u_bar * ∂x(b_bar))
advection_background = Field(@at loc -input_fields.U * ∂x(b_bar))
advection = (; advection_x, advection_background)

turbulent_flux_density_x = Field(ub_bar - flux_density_x)
turbulent_flux_density_z = Field(wb_bar - flux_density_z)
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_z)

mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_z = Field(-∂z(turbulent_flux_density_z))
𝒟_mld = Field(ML_Average(mixing_x, constant_mld))
mixing = (; mixing_x, mixing_z, 𝒟_mld)

mld_flux_density_x = Field(UcFlux(centered, u_mld, b_mld))
flux_density_x_mld = Field(ML_Average(flux_density_x, constant_mld))

mld_advection_x = Field(@at (loc[1], nothing, nothing) -u_mld * ∂x(b_mld))
mld_advection_background = Field(ML_Average(advection_background, constant_mld))
mld_advection = (; mld_flux_density_x, flux_density_x_mld, mld_advection_x, mld_advection_background)

shear_dispersion_flux = Field(flux_density_x_mld - mld_flux_density_x)
𝒟_SD = Field(-∂x(shear_dispersion_flux))
shear_dispersion = (; shear_dispersion_flux, 𝒟_SD)

w_at_mld = Field(ML_Interpolate(w_bar, constant_mld))
b_at_mld = Field(ML_Interpolate(b_bar, constant_mld))
mixing_average =  Field(ML_Average(mixing_z, constant_mld))

F_h = Field(growth_rate * (b_at_mld - b_mld) - w_at_mld * b_mld / constant_mld + mixing_average)
entrainment = (; w_at_mld, mixing_average, b_at_mld, F_h)

Γ = Field(∂x(b_mld)^2 / 2)
divergence = Field(-∂x(u_mld) * ∂x(b_mld)^2)
strain = Field(-∂x(U) * ∂x(b_mld)^2)
dispersion = Field(∂x(𝒟_SD) * ∂x(b_mld))
horizontalmixing = Field(∂x(𝒟_mld) * ∂x(b_mld))
vertical = Field(∂x(F_h) * ∂x(b_mld))
gradient = (; Γ, divergence, strain, dispersion, horizontalmixing, vertical)

skip_update = filter(a->a ∉ fields, keys(input_fields))
dependency_fields = merge(
    mean_fields,
    mld_fields, 
    flux_density,
    advection,
    turbulent_flux_density, 
    mixing, 
    mld_advection, 
    shear_dispersion, 
    entrainment,
    gradient
)
output_fields = (; u_mld, v_mld, b_mld, mld, constant_mld, mld_advection_x, mld_advection_background, 𝒟_mld, 𝒟_SD, F_h, gradient...)
