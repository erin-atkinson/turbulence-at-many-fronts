include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/bbalance/turbulence.jl")

frames = 1:10:1000

# Average
u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)
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
coarse = (; u_coarse, v_coarse, w_coarse, b_coarse)
dependency_fields = merge(dependency_fields, coarse)
@info coarse

# Fluxes
ub = Field(KernelFunctionOperation{Face, Nothing, Center}(fGg, grid, u_coarse, ℑxᶠᵃᵃ, b_coarse))
wb = Field(KernelFunctionOperation{Center, Nothing, Face}(fGg, grid, w_coarse, ℑzᵃᵃᶠ, b_coarse))
Ub = Field(KernelFunctionOperation{Face, Nothing, Center}(fGg, grid, U, ℑxᶠᵃᵃ, b_coarse))
fluxes = (; ub, wb, Ub)
dependency_fields = merge(dependency_fields, fluxes)
@info fluxes

# Turbulent fluxes
ub_full = Field(KernelFunctionOperation{Face, Center, Center}(ub_func, grid, clock, input_fields, dependency_fields, sp))
wb_full = Field(KernelFunctionOperation{Center, Center, Face}(wb_func, grid, clock, input_fields, dependency_fields, sp))
println("Full")

ub_dfm = dfm(ub_full)
wb_dfm = dfm(wb_full)
println("DFM")

ub_coarse = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size, kernel_σ, ub_dfm))
wb_coarse = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, wb_dfm))
println("Coarse")

Fx = Field(ub_coarse - (ub + Ub))
Fz = Field(wb_coarse - wb)

turbulent_fluxes = (; ub_full, wb_full, ub_dfm, wb_dfm, ub_coarse, wb_coarse, Fx, Fz)
dependency_fields = merge(dependency_fields, turbulent_fluxes)
@info turbulent_fluxes

# Gradients
∂b∂x = Field(∂x(b_coarse)) # FNC
∂b∂z = Field(∂z(b_coarse)) # CNF
gradients = (; ∂b∂x, ∂b∂z)
dependency_fields = merge(dependency_fields, gradients)
@info gradients

# Advection by the mean
∂ub∂x = Field(∂x(ub)) # CNC
∂wb∂z = Field(∂z(wb)) # CNC
U∂b∂x = Field(∂x(Ub) - ∂x(U) * b_coarse)
mean_adv = (; ∂ub∂x, ∂wb∂z, U∂b∂x)
dependency_fields = merge(dependency_fields, mean_adv)
@info mean_adv

# Gradients of turbulent fluxes
∂Fx∂x = Field(∂x(Fx))
∂Fz∂z = Field(∂z(Fz))
turbulent_adv = (; ∂Fx∂x, ∂Fz∂z)
dependency_fields = merge(dependency_fields, turbulent_adv)
@info turbulent_adv

# Diffusivity estimates
κₕ = Field(-Fx / ∂b∂x)
κᵥ = Field(-Fz / ∂b∂z)
diffusivities = (; κₕ, κᵥ)
dependency_fields = merge(dependency_fields, diffusivities)
@info diffusivities

output_fields = merge(mean_adv, turbulent_adv, diffusivities)
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
