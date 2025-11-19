include("terms/terms.jl")
include("terms/richardson/ri.jl")

u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
b_dfm = dfm(input_fields.b)
mean_fields = (; u_dfm, v_dfm, b_dfm)


N² = Field(KernelFunctionOperation{Center, Nothing, Face}(N²_func, grid, clock, input_fields, mean_fields, sp))
M² = Field(KernelFunctionOperation{Face, Nothing, Center}(M²_func, grid, clock, input_fields, mean_fields, sp))

Sx = Field(KernelFunctionOperation{Center, Nothing, Face}(Sx_func, grid, clock, input_fields, mean_fields, sp))
Sy = Field(KernelFunctionOperation{Center, Nothing, Face}(Sy_func, grid, clock, input_fields, mean_fields, sp))

Ri = Field(KernelFunctionOperation{Center, Nothing, Center}(Ri_func, grid, clock, input_fields, mean_fields, sp))
Rib = Field(KernelFunctionOperation{Center, Nothing, Center}(Rib_func, grid, clock, input_fields, mean_fields, sp))

dependency_fields = (; mean_fields..., N², M², Sx, Sy, Ri, Rib)
output_fields = (; N², M², Sx, Sy, Ri, Rib)

skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next, :w)
