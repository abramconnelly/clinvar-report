#!/bin/bash

echo ""
echo "********************************************************************************"
echo "This script will now download the latest clinvar.vcf.gz data, decompress it,"
echo "and put it in the data/ directory...."
echo "********************************************************************************"
echo ""

pushd data
./grab_clinvar.sh
popd

echo ""
echo "********************************************************************************"
echo "Assuming the clinvar.vcf file is downloaded and ready, you can generate a test"
echo "report on the sample data provided by issuing the following command:"
echo ""
echo "  ./src/clinvar-report.py -C data/clinvar.vcf -i data/abe.vcf"
echo ""
echo "********************************************************************************"


