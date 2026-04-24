include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

include("terms/vorticity/vorticity.jl")

include("terms/pv/q.jl")
include("terms/pv/Jx.jl")
include("terms/pv/Jz.jl")

frames = frames[1:10:end]

# Average
@info "Down-front averaged fields"
u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)

v_next_dfm = dfm(input_fields.v_next)
b_next_dfm = dfm(input_fields.b_next)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm, v_next_dfm, b_next_dfm)
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
w_coarse = Field(Coarse(w_dfm, kernel))
b_coarse = Field(Coarse(b_dfm, kernel))

v_next_coarse = Field(Coarse(v_next_dfm, kernel))
b_next_coarse = Field(Coarse(b_next_dfm, kernel))

v_avg_coarse = Field((v_coarse + v_next_coarse) / 2)
b_avg_coarse = Field((b_coarse + b_next_coarse) / 2)

coarse_fields = (; u_coarse, v_coarse, w_coarse, b_coarse, v_next_coarse, b_next_coarse, v_avg_coarse, b_avg_coarse)
dependency_fields = merge(dependency_fields, coarse_fields)
@info coarse_fields

# Total fluxes
@info "Total advective fluxes"
(uv, wv, ub, wb) = let u = input_fields.u, 
    v = input_fields.v, w = input_fields.w,
    U = input_fields.U, b = input_fields.b
    
    uv = Field(UvFlux(weno, u, v; background=U))
    wv = Field(WvFlux(weno, w, v))
    
    ub = Field(UcFlux(weno, u, b; background=U))
    wb = Field(WcFlux(weno, w, b))
    (uv, wv, ub, wb)
end
total_fluxes = (; uv, wv, ub, wb)
dependency_fields = merge(dependency_fields, total_fluxes)
@info total_fluxes

# Averaged
@info "Averaged total advective fluxes"
uv_dfm = dfm(uv)
wv_dfm = dfm(wv)

ub_dfm = dfm(ub)
wb_dfm = dfm(wb)
mean_fluxes = (; uv_dfm, wv_dfm, ub_dfm, wb_dfm)
dependency_fields = merge(dependency_fields, mean_fluxes)
@info mean_fluxes

# Coarse-grained fluxes
@info "Coarse-grained total advective fluxes"
uv_coarse = Field(Coarse(uv_dfm, kernel))
wv_coarse = Field(Coarse(wv_dfm, kernel))

ub_coarse = Field(Coarse(ub_dfm, kernel))
wb_coarse = Field(Coarse(wb_dfm, kernel))
coarse_fluxes = (; uv_coarse, wv_coarse, ub_coarse, wb_coarse)
dependency_fields = merge(dependency_fields, coarse_fluxes)
@info coarse_fluxes

# Advective fluxes
@info "Advective fluxes due to coarse-grained fields"
coarse_uv = Field(UvFlux(centered, u_coarse, v_coarse; background=U))
coarse_wv = Field(WvFlux(centered, w_coarse, v_coarse))

coarse_ub = Field(UcFlux(centered, u_coarse, b_coarse; background=U))
coarse_wb = Field(WcFlux(centered, w_coarse, b_coarse))

fluxes = (; coarse_uv, coarse_wv, coarse_ub, coarse_wb)
dependency_fields = merge(dependency_fields, fluxes)
@info fluxes

# Turbulent fluxes
@info "Turbulent fluxes"
u′v′ = Field(uv_coarse - coarse_uv)
w′v′ = Field(wv_coarse - coarse_wv)

u′b′ = Field(ub_coarse - coarse_ub)
w′b′ = Field(wb_coarse - coarse_wb)

turbulent_fluxes = (; u′v′, w′v′, u′b′, w′b′)
dependency_fields = merge(dependency_fields, turbulent_fluxes)
@info turbulent_fluxes

# Momentum and buoyancy forcing
@info "Momentum and buoyancy forcing"
Fy = Field(∂x(u′v′) + ∂z(w′v′))
𝒟 = Field(∂x(u′b′) + ∂z(w′b′))

forcing = (; Fy, 𝒟)
dependency_fields = merge(dependency_fields, forcing)
@info forcing

# Potential vorticity
@info "Potential vorticity"
q = Field(PV(sp, u_coarse, v_coarse, w_coarse, b_coarse))
J_adv_x = Field(PV_J_adv_x(sp, q, u_coarse, input_fields.U))
J_mix_x = Field(PV_J_mix_x(sp, v_avg_coarse, b_avg_coarse, Fy, 𝒟))

J_adv_z = Field(PV_J_adv_z(sp, q, w_coarse))
J_mix_z = Field(PV_J_mix_z(sp, v_avg_coarse, b_avg_coarse, Fy, 𝒟))

Q_adv = Field(∂x(J_adv_x) + ∂z(J_adv_z))
Q_mix = Field(∂x(J_mix_x) + ∂z(J_mix_z))

potential_vorticity = (; q, J_adv_x, J_mix_x, J_adv_z, J_mix_z, Q_adv, Q_mix)
dependency_fields = merge(dependency_fields, potential_vorticity)
@info potential_vorticity

output_fields = potential_vorticity
skip_update = (:pNHS, :u_next, :w_next, :pNHS_next)
