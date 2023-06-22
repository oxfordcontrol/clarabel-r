#!/bin/sh -e

## Author: Hiroaki Yutani in the string2path package

cargo vendor

## Some files in Lapack are missing newlines resulting in notes.
## So fix them.
## define the directory to start from
startdir=./vendor/lapack-sys/lapack

## -type f means 'find files only', not directories
## -print0 and -0 handle file names with spaces
## -name "*.f" -o -name "*.c" matches files ending with .f or .c
find "$startdir" -type f \( -name "*.f" -o -name "*.c" \) -print0 | while IFS= read -r -d '' file; do
    # if the file does not end with a newline
    if [ "$(tail -c 1 "$file" | wc -l)" -eq 0 ]; then
        # append a newline to the end of the file
        printf '\n' >> "$file"
	printf "Appended new line to $file\n"	
    fi
done

# c.f. https://reproducible-builds.org/docs/archives/
gtar \
  --sort=name \
  --mtime='1970-01-01 00:00:00Z' \
  --owner=0 \
  --group=0 \
  --numeric-owner \
  --xz \
  --create \
  --file=vendor.tar.xz \
  vendor

echo
echo
echo "#############################################"
echo "#                                           #"
echo "#  UPDATE src/cargo_vendor_config.toml !!!  #"
echo "#                                           #"
echo "#############################################"
