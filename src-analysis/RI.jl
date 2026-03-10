include("terms/terms.jl")
include("terms/richardson/ri.jl")

frames = frames[1:10:1000]

u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
b_dfm = dfm(input_fields.b)
mean_fields = (; u_dfm, v_dfm, b_dfm)

const kernel_size = (40 * sp.Nh) ÷ 768
const kernel_σ = sp.Lh / 75

u_coarse = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size, kernel_σ, u_dfm))
v_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, v_dfm))
b_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, b_dfm))
coarse_fields = (; u_coarse, v_coarse, b_coarse)

Ri_dfm = Field(KernelFunctionOperation{Center, Nothing, Center}(Ri_func, grid, clock, sp, u_dfm, v_dfm, b_dfm))
Ri_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(Ri_func, grid, clock, sp, u_coarse, v_coarse, b_coarse))

Rib_dfm = Field(KernelFunctionOperation{Center, Nothing, Center}(Rib_func, grid, clock, sp, b_dfm))
Rib_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(Rib_func, grid, clock, sp, b_coarse))

output_fields = (; Ri_dfm, Rib_dfm, Ri_coarse, Rib_coarse)
dependency_fields = merge(mean_fields, coarse_fields, output_fields)

skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next, :w)
