function update_fields!(fields, fds, clock, frame; skip_update)

    fieldnames = Symbol.(keys(fds.fields))
    for fieldname in fieldnames
        fieldname ∈ skip_update && continue
        set!(fields[fieldname], fds[fieldname][frame])
        fill_halo_regions!(fields[fieldname])
    end
    
    # Set previous state
    f2 = max(frame - 1, 1)
    for fieldname in (:u, :v, :w, :b)
        prevfieldname = Symbol(fieldname, :_prev)
        prevfieldname ∈ skip_update && continue
        set!(fields[prevfieldname], fds[fieldname][f2])
        fill_halo_regions!(fields[prevfieldname])
    end

    compute_background!(fields.U, fields.V, fields.W, clock)

    return nothing
end
