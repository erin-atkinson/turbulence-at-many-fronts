include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/ubalance/turbulence.jl")

frames = 1:10:1000

# Average
u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)
pNHS_dfm = dfm(input_fields.pNHS)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm, pNHS_dfm)
dependency_fields = mean_fields
@info mean_fields

# Build weights
const kernel_size = (40 * sp.Nh) ÷ 768
const kernel_σ = sp.Lh / 75

# Coarse-grained fields
u_coarse = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size, kernel_σ, u_dfm))
v_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, v_dfm))
w_coarse = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, w_dfm))
b_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, b_dfm))
pNHS_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, pNHS_dfm))
coarse = (; u_coarse, v_coarse, w_coarse, b_coarse, pNHS_coarse)
dependency_fields = merge(dependency_fields, coarse)
@info coarse

# Fluxes
uu = Field(KernelFunctionOperation{Center, Nothing, Center}(FfGg, grid, ℑxᶜᵃᵃ, u_coarse, ℑxᶜᵃᵃ, u_coarse))
wu = Field(KernelFunctionOperation{Face, Nothing, Face}(FfGg, grid, ℑxᶠᵃᵃ, w_coarse, ℑzᵃᵃᶠ, u_coarse))
Uu = Field(KernelFunctionOperation{Center, Nothing, Center}(FfGg, grid, ℑxᶜᵃᵃ, U, ℑxᶜᵃᵃ, u_coarse))
fluxes = (; uu, wu, Uu)
dependency_fields = merge(dependency_fields, fluxes)
@info fluxes

# Turbulent fluxes
uu_full = Field(KernelFunctionOperation{Center, Center, Center}(uu_func, grid, clock, input_fields, dependency_fields, sp))
wu_full = Field(KernelFunctionOperation{Face, Center, Face}(wu_func, grid, clock, input_fields, dependency_fields, sp))
println("Full")

uu_dfm = dfm(uu_full)
wu_dfm = dfm(wu_full)
println("DFM")

uu_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, uu_dfm))
wu_coarse = Field(KernelFunctionOperation{Face, Nothing, Face}(coarse_grain_variable_x, grid, Face(), kernel_size, kernel_σ, wu_dfm))
println("Coarse")

Fx = Field(uu_coarse - (uu + Uu))
Fz = Field(wu_coarse - wu)

turbulent_fluxes = (; uu_full, wu_full, uu_dfm, wu_dfm, uu_coarse, wu_coarse, Fx, Fz)
dependency_fields = merge(dependency_fields, turbulent_fluxes)
@info turbulent_fluxes

# Gradients
∂u∂x = Field(∂x(u_coarse)) # CNC
∂u∂z = Field(∂z(u_coarse)) # FNF
gradients = (; ∂u∂x, ∂u∂z)
dependency_fields = merge(dependency_fields, gradients)
@info gradients

# Advection by the mean
∂uu∂x = Field(∂x(uu)) # FNC
∂wu∂z = Field(∂z(wu)) # FNC
U∂u∂x = Field(∂x(Uu) - ∂x(U) * u_coarse)
mean_adv = (; ∂uu∂x, ∂wu∂z, U∂u∂x)
dependency_fields = merge(dependency_fields, mean_adv)
@info mean_adv

# Gradients of turbulent fluxes
∂Fx∂x = Field(∂x(Fx))
∂Fz∂z = Field(∂z(Fz))
turbulent_adv = (; ∂Fx∂x, ∂Fz∂z)
dependency_fields = merge(dependency_fields, turbulent_adv)
@info turbulent_adv

# Diffusivity estimates
ν₁ₕ = Field(-Fx / ∂u∂x)
ν₁ᵥ = Field(-Fz / ∂u∂z)
diffusivities = (; ν₁ₕ, ν₁ᵥ)
dependency_fields = merge(dependency_fields, diffusivities)
@info diffusivities

# Pressure
pHS_coarse = Field(CumulativeIntegral(b_coarse; dims=3))
p_coarse = Field(@at (Center, Nothing, Center) pHS_coarse + pNHS_coarse)
pressure = (; pHS_coarse, p_coarse)
dependency_fields = merge(dependency_fields, pressure)
@info pressure 

px = Field(∂x(p_coarse))
others = (; px)
dependency_fields = merge(dependency_fields, others)
@info others

output_fields = merge(mean_adv, turbulent_adv, diffusivities, others)
skip_update = (:u_next, :v_next, :w_next, :b_next, :pNHS_next)
