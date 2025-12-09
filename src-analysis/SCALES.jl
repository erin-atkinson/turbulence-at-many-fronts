include("terms/terms.jl")

# Average
u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)

# Build weights
kernel_size = 1024
kernel_size_005 = 128

#Δx = -kernel_size:kernel_size .* minimum_xspacing(grid)
Δz = -4:4 .* minimum_zspacing(grid)

#weights_x = gaussian_weights(Δx, sp.L)
weights_z = gaussian_weights(Δz, 2)

#@info weights_x
@info weights_z

uz = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_z, grid, weights_z, u_dfm))
vz = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_z, grid, weights_z, v_dfm))
wz = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_z, grid, weights_z, w_dfm))
bz = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_z, grid, weights_z, b_dfm))

coarsez = (; uz, vz, wz, bz)

# 0.1LD
u_coarse = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size, 0.1sp.L, uz))
v_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.1sp.L, vz))
w_coarse = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.1sp.L, wz))
b_coarse = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size, 0.1sp.L, bz))

coarse = (; u_coarse, v_coarse, w_coarse, b_coarse)

# 0.01LD
u005 = Field(KernelFunctionOperation{Face, Nothing, Center}(coarse_grain_variable_x, grid, Face(), kernel_size_005, 0.01sp.L, uz))
v005 = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size_005, 0.01sp.L, vz))
w005 = Field(KernelFunctionOperation{Center, Nothing, Face}(coarse_grain_variable_x, grid, Center(), kernel_size_005, 0.01sp.L, wz))
b005 = Field(KernelFunctionOperation{Center, Nothing, Center}(coarse_grain_variable_x, grid, Center(), kernel_size_005, 0.01sp.L, bz))

medfine = (; u005, v005, w005, b005)

# Going to call this the submesoscale 5m < L < 100m
u_sms = Field(u005 - u_coarse)
v_sms = Field(v005 - v_coarse)
w_sms = Field(w005 - w_coarse)
b_sms = Field(b005 - b_coarse)

sms = (; u_sms, v_sms, w_sms, b_sms)

dependency_fields = merge(mean_fields, coarsez, coarse, medfine, sms)
output_fields = merge(coarse, sms)
skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
