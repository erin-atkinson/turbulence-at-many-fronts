include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

include("terms/mixedlayer/mld.jl")

# Average
@info "Down-front averaged fields"
u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
b_dfm = dfm(input_fields.b)

mean_fields = (; u_dfm, v_dfm, b_dfm)
dependency_fields = mean_fields
@info mean_fields

# Filter widths
σx = coarse_σx
σz = coarse_σz

# Coarse-grained fields
@info "Coarse-grained fields"
kernel = Gaussian(grid, σx, 1.0, σz)
u_coarse = Field(Coarse(u_dfm, kernel))
v_coarse = Field(Coarse(v_dfm, kernel))
b_coarse = Field(Coarse(b_dfm, kernel))

coarse_fields = (; u_coarse, v_coarse, b_coarse)
dependency_fields = merge(dependency_fields, coarse_fields)
@info coarse_fields

@info "MLD averaged"
ε = sp.Δb/6
mld = Field(MLD(b_coarse, ε))

u_coarse_mld = Field(ML_Average(u_coarse, mld))
v_coarse_mld = Field(ML_Average(v_coarse, mld))
b_coarse_mld = Field(ML_Average(b_coarse, mld))

mld_fields = (; mld, u_coarse_mld, v_coarse_mld, b_coarse_mld)
dependency_fields = merge(dependency_fields, mld_fields)
@info mld_fields

output_fields = mld_fields
skip_update = (:w, :pNHS, :u_next, :w_next, :b_next, :pNHS_next)
