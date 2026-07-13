function update_fields!(fields, fds, clock, frame; skip_update)

    fieldnames = Symbol.(keys(fds.fields))
    for fieldname in fieldnames
        fieldname ∈ skip_update && continue
        set!(fields[fieldname], fds[fieldname][frame])
    end
    
    # Set next state
    l = length(fds.u)
    f2 = min(frame + 1, l)
    :u_next in skip_update || set!(fields.u_next, fds.u[f2])
    :v_next in skip_update || set!(fields.v_next, fds.v[f2])
    :w_next in skip_update || set!(fields.w_next, fds.w[f2])
    :b_next in skip_update || set!(fields.b_next, fds.b[f2])

    compute_background!(fields.U, fields.V, fields.W, clock)

    return nothing
end
