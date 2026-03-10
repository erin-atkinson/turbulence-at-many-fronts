include("terms/terms.jl")
include("terms/vorticity/vorticity.jl")
include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/pv/q.jl")
include("terms/pv/uq.jl")
include("terms/pv/Vq.jl")

u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)
mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)
dependency_fields = mean_fields

const kernel_size = (40 * sp.Nh) ÷ 768
const kernel_σ = sp.Lh / 75

u_coarse = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size, kernel_σ, u_dfm))
v_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, v_dfm))
w_coarse = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, w_dfm))
b_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, kernel_σ, b_dfm))
coarse_fields = (; u_coarse, v_coarse, w_coarse, b_coarse)


coarse_q = Field(KernelFunctionOperation{Face, Nothing, Face}(q_func, grid, clock, sp, u_coarse, v_coarse, w_coarse, b_coarse))
dfm_q = Field(KernelFunctionOperation{Face, Nothing, Face}(q_func, grid, clock, sp, u_dfm, v_dfm, w_dfm, b_dfm))
q = Field(KernelFunctionOperation{Face, Face, Face}(q_func, grid, clock, sp, u, v, w, b))
q_dfm = dfm(q)
q_fields = (; coarse_q, dfm_q, q, q_dfm)
dependency_fields = merge(dependency_fields, q_fields)

output_fields = q_fields
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
