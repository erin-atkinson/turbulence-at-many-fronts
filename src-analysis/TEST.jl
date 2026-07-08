include("terms/terms.jl")
include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

frames = frames[1:10:end]

# Filter widths
σx = coarse_σx
σz = coarse_σz

kernel = Gaussian(grid, σx, 1.0, σz)

@info "Sponge layer"
sponge = Field(KernelFunctionOperation{Center, Center, Center}(b_forcing_func, grid, clock, input_fields.b, sp))
sponge_slice = Field(KernelFunctionOperation{Center, Nothing, Center}(b_forcing_func, grid, clock, input_fields.b, sp))
sponge_dfm = dfm(sponge)
sponge_coarse = Field(Coarse(sponge_dfm, kernel))

sponge_layer = (; sponge, sponge_slice, sponge_dfm, sponge_coarse)
dependency_fields = sponge_layer
@info sponge_layer

output_fields = (; sponge_slice, sponge_dfm, sponge_coarse)
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
