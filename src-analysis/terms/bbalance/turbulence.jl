@inline function ub_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    U = fields.U
    V = fields.V
    W = fields.W

#    ub = @inbounds dependency_fields.ub[i, j, k]
#    Ub = @inbounds dependency_fields.Ub[i, j, k]

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    return advective_tracer_flux_density_x(i, j, k, grid, weno, total_velocities.u, b)# - (ub + Ub)
end

@inline function wb_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    U = fields.U
    V = fields.V
    W = fields.W

#    wb = @inbounds dependency_fields.wb[i, j, k]

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    return advective_tracer_flux_density_z(i, j, k, grid, weno, total_velocities.w, b)# - wb
end

