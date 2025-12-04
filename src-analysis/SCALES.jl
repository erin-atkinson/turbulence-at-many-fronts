include("terms/terms.jl")

# Average
u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)

# Build weights
kernel_size = 1024

Δx = -kernel_size:kernel_size .* minimum_xspacing(grid)
Δz = -4:4 .* minimum_zspacing(grid)

weights_x = gaussian_weights(Δx, sp.L)
weights_z = gaussian_weights(Δz, 2)

@info weights_x
@info weights_z

uz = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_z, grid, weights_z, u_dfm))
vz = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_z, grid, weights_z, v_dfm))
wz = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_z, grid, weights_z, w_dfm))
bz = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_z, grid, weights_z, b_dfm))

coarsez = (; uz, vz, wz, bz)

# 0.05LD
u01 = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size, 0.1sp.L, uz))
v01 = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.1sp.L, vz))
w01 = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.1sp.L, wz))
b01 = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.1sp.L, bz))

coarse01 = (; u01, v01, w01, b01)

# 0.25LD
u05 = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size, 0.5sp.L, uz))
v05 = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.5sp.L, vz))
w05 = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.5sp.L, wz))
b05 = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.5sp.L, bz))

coarse05 = (; u05, v05, w05, b05)

# 1.0LD
u10 = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size, 0.5sp.L, uz))
v10 = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.5sp.L, vz))
w10 = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.5sp.L, wz))
b10 = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.5sp.L, bz))

coarse10 = (; u10, v10, w10, b10)

dependency_fields = merge(mean_fields, coarsez, coarse01, coarse05, coarse10)
output_fields = merge(coarse01, coarse05, coarse10)
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
