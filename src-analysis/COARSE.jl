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
u_coarse = Field(Coarse(u_dfm, Gaussian(grid, σx, 1.0, σz)))
v_coarse = Field(Coarse(v_dfm, Gaussian(grid, σx, 1.0, σz)))
w_coarse = Field(Coarse(w_dfm, Gaussian(grid, σx, 1.0, σz)))
b_coarse = Field(Coarse(b_dfm, Gaussian(grid, σx, 1.0, σz)))
coarse = (; u_coarse, v_coarse, w_coarse, b_coarse)
dependency_fields = merge(dependency_fields, coarse)
@info coarse

output_fields = coarse
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
