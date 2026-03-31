include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

include("terms/energy/tke.jl")
include("terms/energy/dsp.jl")
include("terms/energy/lsp.jl")
include("terms/energy/vsp.jl")

# Average
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
σx = coarse_σx
σz = coarse_σz

# Coarse-grained fields
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
(uu, uv, uw, wu, wv, ww, wb) = let u = input_fields.u, 
    v = input_fields.u, w = input_fields.w,
    U = input_fields.U, b = input_fields.b
    
    uu = Field(UuFlux(u, u; background=U))
    uv = Field(UvFlux(u, v; background=U))
    uw = Field(UwFlux(u, w; background=U))
    
    wu = Field(WuFlux(w, u))
    wv = Field(WvFlux(w, v))
    ww = Field(WwFlux(w, w))
    
    wb = Field(WcFlux(w, b))
    (uu, uv, uw, wu, wv, ww, wb)
end
total_fluxes = (; uu, uv, uw, wu, wv, ww, wb)
@info total_fluxes

# Averaged
uu_dfm = dfm(uu)
uv_dfm = dfm(uv)
uw_dfm = dfm(uw)

wu_dfm = dfm(wu)
wv_dfm = dfm(wv)
ww_dfm = dfm(ww)
wb_dfm = dfm(wb)
mean_fluxes = (; uu_dfm, uv_dfm, uw_dfm, wu_dfm, wv_dfm, ww_dfm, wb_dfm)
@info mean_fluxes

# Coarse-grained fluxes
uu_coarse = Field(Coarse(uu, kernel))
uv_coarse = Field(Coarse(uv, kernel))
uw_coarse = Field(Coarse(uw, kernel))

wu_coarse = Field(Coarse(wu, kernel))
wv_coarse = Field(Coarse(wv, kernel))
ww_coarse = Field(Coarse(ww, kernel))
bflux = Field(Coarse(wb, kernel))
coarse_fluxes = (; uu_coarse, uv_coarse, uw_coarse, wu_coarse, wv_coarse, ww_coarse, bflux)
@info coarse_fluxes



LSP3D = Field(KernelFunctionOperation{Center, Center, Center}(LSP3D_func, grid, clock, coarse_fields, coarse_fluxes, sp))
VSP3D = Field(KernelFunctionOperation{Center, Center, Center}(VSP3D_func, grid, clock, input_fields, coarse_fields, sp))
BFLUX3D = Field(KernelFunctionOperation{Center, Center, Center}(BFLUX3D_func, grid, clock, input_fields, coarse_fields, sp))
DSP3D = Field(KernelFunctionOperation{Center, Center, Center}(DSP3D_func, grid, clock, input_fields, coarse_fields, sp))
TKE3D = Field(KernelFunctionOperation{Center, Center, Center}(TKE3D_func, grid, clock, input_fields, coarse_fields, sp))
DTKEDt3D = Field(KernelFunctionOperation{Center, Center, Center}(DTKEDt3D_func, grid, clock, input_fields, coarse_fields, sp))

TKE3D_fields = (; LSP3D, VSP3D, BFLUX3D, DSP3D, TKE3D, DTKEDt3D)

(LSP, VSP, BFLUX, DSP, TKE, DTKEDt) = map(dfm, TKE3D_fields)

TKE_fields = (; LSP, VSP, BFLUX, DSP, TKE, DTKEDt)

ε = Field(KernelFunctionOperation{Center, Nothing, Center}(ε_func, grid, clock, input_fields, TKE_fields, sp))

TKE_fields = merge(TKE_fields, (; ε))

dependency_fields = merge(mean_fields, TKE3D_fields, TKE_fields)
output_fields = TKE_fields

skip_update = (:pNHS, :b_next, :pNHS_next)