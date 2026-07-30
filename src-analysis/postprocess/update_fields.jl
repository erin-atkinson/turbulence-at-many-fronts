function update_fields!(fields, fds, clock, frame; skip_update)

    fieldnames = Symbol.(keys(fds.fields))
    for fieldname in fieldnames
        fieldname ∈ skip_update && continue
        set!(fields[fieldname], fds[fieldname][frame])
    end
    
    # Set previous state
    f2 = max(frame - 1, 1)
    :u_prev in skip_update || set!(fields.u_prev, fds.u[f2])
    :v_prev in skip_update || set!(fields.v_prev, fds.v[f2])
    :w_prev in skip_update || set!(fields.w_prev, fds.w[f2])
    :b_prev in skip_update || set!(fields.b_prev, fds.b[f2])

    compute_background!(fields.U, fields.V, fields.W, clock)

    return nothing
end
