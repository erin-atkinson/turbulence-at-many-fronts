include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

include("terms/mixedlayer/mld.jl")
frames = frames[1:10:end]
@info "Down-front averaged fields"
u_dfm = dfm(input_fields.u)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)
b_next_dfm = dfm(input_fields.b_next)

mean_fields = (; u_dfm, w_dfm, b_dfm, b_next_dfm)
dependency_fields = mean_fields
@info mean_fields

# Filter widths
σx = coarse_σx
σz = coarse_σz

@info "Coarse-grained fields"
kernel = Gaussian(grid, σx, 1.0, σz)
u_coarse = Field(Coarse(u_dfm, kernel))
w_coarse = Field(Coarse(w_dfm, kernel))
b_coarse = Field(Coarse(b_dfm, kernel))
b_next_coarse = Field(Coarse(b_next_dfm, kernel))

coarse_fields = (; u_coarse, w_coarse, b_coarse, b_next_coarse)
dependency_fields = merge(dependency_fields, coarse_fields)
@info coarse_fields

@info "MLD averaged"
ε = sp.Δb/6
mld = Field(ConstantMLD(b_coarse, ε))
mld_next = Field(ConstantMLD(b_next_coarse, ε))

u_coarse_mld = Field(ML_Average(u_coarse, mld))
b_coarse_mld = Field(ML_Average(b_coarse, mld))
Γ = Field(∂x(b_coarse_mld))

mld_fields = (; mld, mld_next, u_coarse_mld, b_coarse_mld, Γ)
dependency_fields = merge(dependency_fields, mld_fields)
@info mld_fields

@info "Horizontal total flux"
ub = Field(UcFlux(weno, input_fields.u, input_fields.b; background=input_fields.U))
ub_dfm = dfm(ub)
ub_coarse = Field(Coarse(ub_dfm, kernel))
ub_coarse_mld = Field(ML_Average(ub_coarse, mld))

ub_total = (; ub, ub_dfm, ub_coarse, ub_coarse_mld)
dependency_fields = merge(dependency_fields, ub_total)
@info ub_total

@info "Vertical total flux"
wb = Field(WcFlux(weno, input_fields.w, input_fields.b))
wb_dfm = dfm(wb)
wb_coarse = Field(Coarse(wb_dfm, kernel))

wb_total = (; wb, wb_dfm, wb_coarse)
dependency_fields = merge(dependency_fields, wb_total)
@info wb_total

@info "Advective fluxes due to coarse-grained fields"
coarse_ub = Field(UcFlux(centered, u_coarse, b_coarse; background=input_fields.U))
coarse_wb = Field(WcFlux(centered, w_coarse, b_coarse))

coarse_ub_mld = Field(ML_Average(coarse_ub, mld))

coarse_adv = (; coarse_ub, coarse_wb, coarse_ub_mld)
dependency_fields = merge(dependency_fields, coarse_adv)
@info coarse_adv

@info "Divergence and strain"
divergence = Field(-Γ^2 * ∂x(u_coarse_mld) )
strain = Field(-Γ^2 * ∂x(input_fields.U))

frontogenesis = (; divergence, strain)
dependency_fields = merge(dependency_fields, frontogenesis)
@info frontogenesis

@info "Horizontal mixing"
u′b′_mld = Field(ub_coarse_mld - coarse_ub_mld)
horizontalmixing = Field(-Γ * ∂x(∂x(u′b′_mld)))

horizontal = (; u′b′_mld, horizontalmixing)
dependency_fields = merge(dependency_fields, horizontal)
@info horizontal

@info "Advection due to MLD fields"
advection = Field(Γ * u_coarse_mld * ∂x(Γ))
background = Field(Γ * input_fields.U * ∂x(Γ))

mld_advection = (; advection, background)
dependency_fields = merge(dependency_fields, mld_advection)
@info mld_advection

@info "Shear dispersion"
coarse_mld_ub = Field(UcFlux(centered, u_coarse_mld, b_coarse_mld; background=input_fields.U))
u_dagu_dag_mld = Field(coarse_ub_mld - coarse_mld_ub)
sheardispersion = Field(-Γ * ∂x(∂x(u_dagu_dag_mld)))

dispersion = (; coarse_mld_ub, u_dagu_dag_mld, sheardispersion)
dependency_fields = merge(dependency_fields, dispersion)
@info dispersion

@info "Bottom mixing"
w′b′_at_mld = Field(ML_Interpolate(wb_coarse - coarse_wb, mld))
bottommixing = Field(Γ * ∂x(w′b′_at_mld) / mld)

bottom = (; w′b′_at_mld, bottommixing)
dependency_fields = merge(dependency_fields, bottom)
@info bottom

@info "Entrainment"
b_at_mld = Field(ML_Interpolate(b_coarse, mld))
entrainment = Field(Γ * GrowthRate(clock, mld, mld_next) * ∂x(b_at_mld - b_coarse_mld))

entrain = (; b_at_mld, entrainment)
dependency_fields = merge(dependency_fields, entrain)
@info entrain

output_fields = (; Γ, advection, background, divergence, strain, horizontalmixing, sheardispersion, bottommixing, entrainment)
skip_update = (:v, :pNHS, :u_next, :w_next, :v_next, :pNHS_next)
