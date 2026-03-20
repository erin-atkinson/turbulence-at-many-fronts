include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/vbalance/turbulence.jl")

frames = 1:10:1000

# Average
u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)

mean_fields = (; u_dfm, v_dfm, w_dfm)
dependency_fields = mean_fields
@info mean_fields

# Build weights
const kernel_size = (40 * sp.Nh) ÷ 768
const kernel_σ = sp.Lh / 75

# Coarse-grained fields
u_coarse = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size, kernel_σ, u_dfm))
v_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, v_dfm))
w_coarse = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, w_dfm))
coarse = (; u_coarse, v_coarse, w_coarse)
dependency_fields = merge(dependency_fields, coarse)
@info coarse

# Fluxes
uv = Field(KernelFunctionOperation{Face, Nothing, Center}(fGg, grid, u_coarse, ℑxᶠᵃᵃ, v_coarse))
wv = Field(KernelFunctionOperation{Center, Nothing, Face}(fGg, grid, w_coarse, ℑzᵃᵃᶠ, v_coarse))
Uv = Field(KernelFunctionOperation{Face, Nothing, Center}(fGg, grid, U, ℑxᶠᵃᵃ, v_coarse))
fluxes = (; uv, wv, Uv)
dependency_fields = merge(dependency_fields, fluxes)
@info fluxes

# Turbulent fluxes
uv_full = Field(KernelFunctionOperation{Face, Face, Center}(uv_func, grid, clock, input_fields, dependency_fields, sp))
wv_full = Field(KernelFunctionOperation{Center, Face, Face}(wv_func, grid, clock, input_fields, dependency_fields, sp))
println("Full")

uv_dfm = dfm(uv_full)
wv_dfm = dfm(wv_full)
println("DFM")

uv_coarse = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, uv_dfm))
wv_coarse = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Face(), kernel_size, kernel_σ, wv_dfm))
println("Coarse")

Fx = Field(uv_coarse - (uv + Uv))
Fz = Field(wv_coarse - wv)

turbulent_fluxes = (; uv_full, wv_full, uv_dfm, wv_dfm, uv_coarse, wv_coarse, Fx, Fz)
dependency_fields = merge(dependency_fields, turbulent_fluxes)
@info turbulent_fluxes

# Gradients
∂v∂x = Field(∂x(v_coarse)) # FNC
∂v∂z = Field(∂z(v_coarse)) # CNF
gradients = (; ∂v∂x, ∂v∂z)
dependency_fields = merge(dependency_fields, gradients)
@info gradients

# Advection by the mean
∂uv∂x = Field(∂x(uv)) # CNC
∂wv∂z = Field(∂z(wv)) # CNC
U∂v∂x = Field(∂x(Uv) - ∂x(U) * v_coarse)
mean_adv = (; ∂uv∂x, ∂wv∂z, U∂v∂x)
dependency_fields = merge(dependency_fields, mean_adv)
@info mean_adv

# Gradients of turbulent fluxes
∂Fx∂x = Field(∂x(Fx))
∂Fz∂z = Field(∂z(Fz))
turbulent_adv = (; ∂Fx∂x, ∂Fz∂z)
dependency_fields = merge(dependency_fields, turbulent_adv)
@info turbulent_adv

# Diffusivity estimates
ν₂ₕ = Field(-Fx / ∂v∂x)
ν₂ᵥ = Field(-Fz / ∂v∂z)
diffusivities = (; ν₂ₕ, ν₂ᵥ)
dependency_fields = merge(dependency_fields, diffusivities)
@info diffusivities

output_fields = merge(mean_adv, turbulent_adv, diffusivities)
skip_update = (:b, :pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
