#!/bin/sh -e

## Author: Hiroaki Yutani in the string2path package
## Modifications: naras@stanford

cargo vendor

## Clarabel has a CITATION.bib file that CRAN flags
# c.f. https://reproducible-builds.org/docs/archives/
echo
echo "Tarring up files"
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
echo "##############################################"
echo "#                                            #"
echo "#  UPDATE src/cargo_vendor_config.toml !!!   #"
echo "#                                            #"
echo "##############################################"
echo
