include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

include("terms/energy/ke.jl")
include("terms/energy/pe.jl")

include("terms/energy/dsp.jl")
include("terms/energy/lsp.jl")
include("terms/energy/vsp.jl")

frames = frames[1:10:end]

# Average
@info "Down-front averaged fields"
u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)

u_next_dfm = dfm(input_fields.u_next)
v_next_dfm = dfm(input_fields.v_next)
w_next_dfm = dfm(input_fields.w_next)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm, u_next_dfm, v_next_dfm, w_next_dfm)
dependency_fields = mean_fields
@info mean_fields

# Filter widths
σx = sp.L
σz = coarse_σz

# Coarse-grained fields
@info "Coarse-grained fields"
kernel = Gaussian(grid, σx, 1.0, σz)
u_coarse = Field(Coarse(u_dfm, kernel))
v_coarse = Field(Coarse(v_dfm, kernel))
w_coarse = Field(Coarse(w_dfm, kernel))
b_coarse = Field(Coarse(b_dfm, kernel))

u_next_coarse = Field(Coarse(u_next_dfm, kernel))
v_next_coarse = Field(Coarse(v_next_dfm, kernel))
w_next_coarse = Field(Coarse(w_next_dfm, kernel))

coarse_fields = (; u_coarse, v_coarse, w_coarse, b_coarse, u_next_coarse, v_next_coarse, w_next_coarse)
dependency_fields = merge(dependency_fields, coarse_fields)
@info coarse_fields

# Total fluxes
@info "Total advective fluxes"
(uu, uv, uw, wu, wv, ww, wb) = let u = input_fields.u, 
    v = input_fields.v, w = input_fields.w,
    U = input_fields.U, b = input_fields.b
    
    uu = Field(UuFlux(weno, u, u; background=U))
    uv = Field(UvFlux(weno, u, v; background=U))
    uw = Field(UwFlux(weno, u, w; background=U))
    
    wu = Field(WuFlux(weno, w, u))
    wv = Field(WvFlux(weno, w, v))
    ww = Field(WwFlux(weno, w, w))
    
    wb = Field(WcFlux(weno, w, b))
    (uu, uv, uw, wu, wv, ww, wb)
end
total_fluxes = (; uu, uv, uw, wu, wv, ww, wb)
dependency_fields = merge(dependency_fields, total_fluxes)
@info total_fluxes

# Averaged
@info "Averaged total advective fluxes"
uu_dfm = dfm(uu)
uv_dfm = dfm(uv)
uw_dfm = dfm(uw)

wu_dfm = dfm(wu)
wv_dfm = dfm(wv)
ww_dfm = dfm(ww)
wb_dfm = dfm(wb)
mean_fluxes = (; uu_dfm, uv_dfm, uw_dfm, wu_dfm, wv_dfm, ww_dfm, wb_dfm)
dependency_fields = merge(dependency_fields, mean_fluxes)
@info mean_fluxes

# Coarse-grained fluxes
@info "Coarse-grained total advective fluxes"
uu_coarse = Field(Coarse(uu_dfm, kernel))
uv_coarse = Field(Coarse(uv_dfm, kernel))
uw_coarse = Field(Coarse(uw_dfm, kernel))

wu_coarse = Field(Coarse(wu_dfm, kernel))
wv_coarse = Field(Coarse(wv_dfm, kernel))
ww_coarse = Field(Coarse(ww_dfm, kernel))
wb_coarse = Field(Coarse(wb_dfm, kernel))
coarse_fluxes = (; uu_coarse, uv_coarse, uw_coarse, wu_coarse, wv_coarse, ww_coarse, wb_coarse)
dependency_fields = merge(dependency_fields, coarse_fluxes)
@info coarse_fluxes

# Advective fluxes
@info "Advective fluxes due to coarse-grained fields"
coarse_uu = Field(UuFlux(centered, u_coarse, u_coarse; background=U))
coarse_uv = Field(UvFlux(centered, u_coarse, v_coarse; background=U))
coarse_uw = Field(UwFlux(centered, u_coarse, w_coarse; background=U))

coarse_wu = Field(WuFlux(centered, w_coarse, u_coarse))
coarse_wv = Field(WvFlux(centered, w_coarse, v_coarse))
coarse_ww = Field(WwFlux(centered, w_coarse, w_coarse))

coarse_wb = Field(WcFlux(centered, w_coarse, b_coarse))
fluxes = (; coarse_uu, coarse_uv, coarse_uw, coarse_wu, coarse_wv, coarse_ww, coarse_wb)
dependency_fields = merge(dependency_fields, fluxes)
@info fluxes

# Turbulent fluxes
@info "Turbulent fluxes"
u′u′ = Field(uu_coarse - coarse_uu)
u′v′ = Field(uv_coarse - coarse_uv)
u′w′ = Field(uw_coarse - coarse_uw)

w′u′ = Field(wu_coarse - coarse_wu)
w′v′ = Field(wv_coarse - coarse_wv)
w′w′ = Field(ww_coarse - coarse_ww)

bflux = Field(@at (Center, Nothing, Center) wb_coarse - coarse_wb)
turbulent_fluxes = (; u′u′, u′v′, u′w′, w′u′, w′v′, w′w′, bflux)
dependency_fields = merge(dependency_fields, turbulent_fluxes)
@info turbulent_fluxes

# Energies
@info "Energies and energy fluxes"
ke = Field(KE(input_fields.u, input_fields.v, input_fields.w))
mke = Field(KE(u_coarse, v_coarse, w_coarse))
ke_dfm = dfm(ke)
ke_coarse = Field(Coarse(ke_dfm, kernel))
tke = Field(ke_coarse - mke)
pe = Field(PE(b_coarse))

lsp, vsp = let velocities = (; u=u_coarse, v=v_coarse, w=w_coarse),
    velocities_next = (; u=u_next_coarse, v=v_next_coarse, w=w_next_coarse),
    fluxes = (; uu=u′u′, uv=u′v′, uw=u′w′, wu=w′u′, wv=w′v′, ww=w′w′)
    lsp = Field(LSP(velocities, velocities_next, fluxes))
    vsp = Field(VSP(velocities, velocities_next, fluxes))

    lsp, vsp
end
dsp = Field(DSP(u_coarse, v_coarse, u_next_coarse, v_next_coarse, input_fields.U))

energies = (; ke, mke, ke_dfm, ke_coarse, tke, pe, lsp, vsp, dsp)
dependency_fields = merge(dependency_fields, energies)
@info energies

output_fields = (; mke, pe, tke, lsp, vsp, dsp, bflux)
skip_update = (:pNHS, :b_next, :pNHS_next)
