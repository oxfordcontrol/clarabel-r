#!/usr/bin/env Rscript

hash_file <- "./vendor/lapack-sys/.cargo-checksum.json"
js <- jsonlite::fromJSON(hash_file)
files <- readLines("./modified_files.txt")
new_hash <- old_hash <- rep("", length(files))
names(new_hash) <- names(old_hash) <- files
for (f in files) {
  # Split file path into components
  path_comps <- rev(gtools::split_path(x = f))
  # Make path relative to package dir by dropping first two
  rel_path <- do.call(file.path, as.list(path_comps[-(1:2)]))
  # Save old hash for troubleshooting and reference
  old_hash[f] <- js$files[[rel_path]]
  new_hash[f] <- as.character(openssl::sha256(file(f)))
  # Set new hash
  js$files[[rel_path]] <- new_hash
}
jsonlite::write_json(x = js, path = hash_file)

## Write csv output
readr::write_csv(data.frame(file = files, old_hash = old_hash, new_hash = new_hash),
                 file = "modified_files.csv")

                 
