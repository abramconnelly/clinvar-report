#!/bin/bash

u="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
echo "Fetching the latest clinvar.vcf.gz file from $u"

wget "$u" -O clinvar.vcf.gz

echo "decompressing..."
gunzip clinvar.vcf.gz

echo "done."
