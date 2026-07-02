include("terms/terms.jl")

#frames = 1:1:1000

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
kernel = Gaussian(grid, σx, 0, σz)
u_coarse = Field(Coarse(u_dfm, kernel))
v_coarse = Field(Coarse(v_dfm, kernel))
w_coarse = Field(Coarse(w_dfm, kernel))
b_coarse = Field(Coarse(b_dfm, kernel))

ψ_coarse = Field(CumulativeIntegral(-u_coarse; dims=3))
coarse = (; u_coarse, v_coarse, w_coarse, b_coarse, ψ_coarse)
dependency_fields = merge(dependency_fields, coarse)
@info coarse

output_fields = coarse
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
