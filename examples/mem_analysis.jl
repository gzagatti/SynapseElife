using Coverage

coevolve = [f for f in analyze_malloc(".") if endswith(f.filename, "636555.mem")]
coevolve = append!(
    coevolve,
    [
        f for f in analyze_malloc("/home/guilherme/.julia") if
        endswith(f.filename, "636555.mem")
    ],
)
coevolve = sort(coevolve; lt = (x, y) -> x.bytes < y.bytes)
@info "Total memory coevolve" sum([f.bytes for f in coevolve]) / 1e6
Base.print_matrix(
    stdout,
    [(f.bytes / 1e6, f.filename, f.linenumber) for f in coevolve if f.bytes > 0],
)

pdmp = [f for f in analyze_malloc(".") if endswith(f.filename, "618695.mem")]
pdmp = append!(
    pdmp,
    [
        f for f in analyze_malloc("/home/guilherme/.julia") if
        endswith(f.filename, "618695.mem")
    ],
)
pdmp = sort(pdmp; lt = (x, y) -> x.bytes < y.bytes)
@info "Total memory pdmp" sum([f.bytes for f in pdmp]) / 1e6
Base.print_matrix(
    stdout,
    [(f.bytes / 1e6, f.filename, f.linenumber) for f in pdmp if f.bytes > 0],
)
