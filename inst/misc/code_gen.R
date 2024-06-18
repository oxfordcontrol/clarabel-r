## TypedSexp::String(s) => settings.direct_solve_method = if let Some(result) = s.to_vec().get(0) {
##     result.to_string()
## } else {
##     settings.direct_solve_method
## }


library(clarabel)
# Function to match on types using savvy
f <- function(name, type = c("integer", "double", "logical", "character")) {
  type <- match.arg(type)
  c(
    sprintf('"%s" => ', name),
    sprintf('   match typed_value {'),
    switch(type,
           integer = sprintf('       TypedSexp::Integer(i) => settings.%s = i.as_slice()[0],', name),
           double = sprintf('       TypedSexp::Real(f) => settings.%s = f.as_slice()[0],', name),
           logical = sprintf('       TypedSexp::Logical(b) => settings.%s = b.as_slice_raw()[0] != 0,', name),
           character = sprintf('       TypedSexp::String(s) => if let Some(result) = s.to_vec().get(0) {\\n
			settings.%s = result.to_string();\\n		    },\\n		    _ => (),', name)
           ),
    "_ => (),"
    "},"
    )
}


s <- clarabel_control()
sn <- c(names(s),
        "chordal_decomposition_enable",
        "chordal_decomposition_merge_method",
        "chordal_decomposition_compact",
        "chordal_decomposition_complete_dual")
st <- c(sapply(s, typeof),
        "logical", "character", "logical", "logical")
out <- lapply(seq_along(st), function(i) f(sn[i], st[i]))
writeLines(unlist(out), "update_settings.rs")
