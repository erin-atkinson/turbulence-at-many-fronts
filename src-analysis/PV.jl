include("terms/terms.jl")
include("terms/vorticity/vorticity.jl")
include("terms/advection/advection.jl")

include("terms/pv/q.jl")
include("terms/pv/uq.jl")
include("terms/pv/Vq.jl")

# Average
u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)
dependency_fields = mean_fields
@info mean_fields

# Filter widths
σx = coarse_σx
σz = coarse_σz

# Coarse-grained fields
u_coarse = Field(Coarse(u_dfm, Gaussian(grid, σx, 1.0, σz)))
v_coarse = Field(Coarse(v_dfm, Gaussian(grid, σx, 1.0, σz)))
w_coarse = Field(Coarse(w_dfm, Gaussian(grid, σx, 1.0, σz)))
b_coarse = Field(Coarse(b_dfm, Gaussian(grid, σx, 1.0, σz)))
coarse = (; u_coarse, v_coarse, w_coarse, b_coarse)
dependency_fields = merge(dependency_fields, coarse)
@info coarse

coarse_q = Field(KernelFunctionOperation{Face, Nothing, Face}(q_func, grid, clock, sp, u_coarse, v_coarse, w_coarse, b_coarse))
dfm_q = Field(KernelFunctionOperation{Face, Nothing, Face}(q_func, grid, clock, sp, u_dfm, v_dfm, w_dfm, b_dfm))
q = Field(KernelFunctionOperation{Face, Face, Face}(q_func, grid, clock, sp, input_fields.u, input_fields.v, input_fields.w, input_fields.b))
q_dfm = dfm(q)
q_fields = (; coarse_q, dfm_q, q, q_dfm)
dependency_fields = merge(dependency_fields, q_fields)

output_fields = (; coarse_q, dfm_q, q_dfm)
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
