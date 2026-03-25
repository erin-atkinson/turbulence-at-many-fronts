include("terms/terms.jl")

#frames = 1:10:1000

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

output_fields = coarse
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
