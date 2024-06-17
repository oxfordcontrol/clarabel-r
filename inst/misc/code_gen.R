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
           character = sprintf('       TypedSexp::String(s) => settings.%s = String::from(s.to_vec().get(0)),', name)
           ),
    sprintf("_ => settings.%s = settings.%s,", name, name),
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
writeLines(unlist(out), "junk.rs")


function(name, type) {
  c(
    sprintf('"%s" => ', name),
    sprintf('   match typed_value() {'),
    if(type == "integer") {
      c(sprintf('       TypedSexp::Integer(i) => settings.%s = i.as_slice()[0],', name),
        sprintf('       _ => savvy::io::r_warn("%s should be scalar integer")', name),
        "},")
    } else {
      c(sprintf('       TypedSexp::Real(f) => settings.%s = f.as_slice()[0],', name),
        sprintf('       _ => savvy::io::r_warn("%s should be scalar double")', name),
        "},")
    }
  )
}
