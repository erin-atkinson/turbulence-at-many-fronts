include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

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

@info "Coarse-grained fields"
kernel = Gaussian(grid, σx, 1.0, σz)
u_coarse = Field(Coarse(u_dfm, kernel))
w_coarse = Field(Coarse(w_dfm, kernel))
b_coarse = Field(Coarse(b_dfm, kernel))

@info "Horizontal total flux"
utotb = Field(UcFlux(weno, input_fields.u, input_fields.b; background=input_fields.U))
utotb_dfm = dfm(utotb)
utotb_coarse = Field(Coarse(utotb_dfm, kernel))

ub_total = (; utotb, utotb_dfm, utotb_coarse)
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
coarse_utotb = Field(UcFlux(centered, u_coarse, b_coarse; background=input_fields.U))
coarse_ub = Field(UcFlux(centered, u_coarse, b_coarse))
coarse_wb = Field(WcFlux(centered, w_coarse, b_coarse))

coarse_adv = (; coarse_utotb, coarse_ub, coarse_wb)
dependency_fields = merge(dependency_fields, coarse_adv)
@info coarse_adv

@info "Turbulent fluxes"
u′b′_coarse = Field(utotb_coarse - coarse_utotb)
w′b′_coarse = Field(wb_coarse - coarse_wb)

turbulent_fluxes = (; u′b′_coarse, w′b′_coarse)
dependency_fields = merge(dependency_fields, turbulent_fluxes)
@info turbulent_fluxes

@info "Advection by the mean"
adv_x = Field(-∂x(coarse_ub))
adv_z = Field(-∂z(coarse_wb))
strain = Field(-∂x(coarse_utotb - coarse_ub) + ∂x(input_fields.U) * b_coarse)
adv = (; adv_x, adv_z, strain)
dependency_fields = merge(dependency_fields, adv)
@info adv

@info "Gradients of turbulent fluxes"
mix_x = Field(-∂x(u′b′_coarse))
mix_z = Field(-∂z(w′b′_coarse))
mix = (; mix_x, mix_z)
dependency_fields = merge(dependency_fields, mix)
@info mix

#=
@info "Sponge layer"
sponge = Field(KernelFunctionOperation{Center, Center, Center}(b_forcing_func, grid, clock, input_fields.b, sp))
sponge_dfm = dfm(sponge)
sponge_coarse = Field(Coarse(sponge_dfm, kernel))

sponge_layer = (; sponge, sponge_dfm, sponge_coarse)
dependency_fields = merge(dependency_fields, sponge_layer)
@info sponge_layer
=#

output_fields = merge(adv, mix, (; coarse_ub, u′b′_coarse, coarse_wb, w′b′_coarse))
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
