@inline function uv_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    return advective_momentum_flux_density_Uv(i, j, k, grid, weno, total_velocities.u, v)
end

@inline function wv_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    return advective_momentum_flux_density_Wv(i, j, k, grid, weno, total_velocities.w, v)
end

