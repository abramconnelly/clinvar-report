clinvar_report
===

Generate a report based on the ClinVar variant database.

This code is a slightly modified and stand alone version of the report generation from the [Genevieve](https://github.com/PersonalGenomesOrg/genevieve) project.

---

Getting Started
----

You'll need a local copy of the `clinvar.vcf` file from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz .
You can either download this yourself or run the `setup.sh` script provided.

```
$ ./setup.sh

********************************************************************************
This script will now download the latest clinvar.vcf.gz data, decompress it,
and put it in the data/ directory....
********************************************************************************
...

Fetching the latest clinvar.vcf.gz file from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
...

decompressing...

...

********************************************************************************
Assuming the clinvar.vcf file is downloaded and ready, you can generate a test
report on the sample data provided by issuing the following command:

  ./src/clinvar_report.py -C data/clinvar.vcf -i data/abe.vcf

********************************************************************************

```

---

Generating a CSV report
----

A sample VCF file is provided, `abe.vcf`, that you can do a test run on:

```
$ ./src/clinvar_report.py -C data/clinvar.vcf -i data/abe.vcf
```

Which will produce output like the following:

```
Chromosome,Position,Name,Significance,Frequency,Zygosity,ACC URL
1,45444038,not_specified,other,0.2965,Het,http://www.ncbi.nlm.nih.gov/clinvar/RCV000081999.4
1,46660295,not_specified,other,0.3335,Het,http://www.ncbi.nlm.nih.gov/clinvar/RCV000081805.4
1,46870761,Drug_addiction\x2c_susceptibility_to,other,0.2616,Het,http://www.ncbi.nlm.nih.gov/clinvar/RCV000007116.2
1,70904800,Homocysteine\x2c_total_plasma\x2c_elevated,pathogenic,0.2053,Het,http://www.ncbi.nlm.nih.gov/clinvar/RCV000003075.1
1,154426970,Serum_level_of_interleukin-6_soluble_receptor,other,0.2931,Het,http://www.ncbi.nlm.nih.gov/clinvar/RCV000015767.1
1,156105028,not_specified,other,0.1931,Het,http://www.ncbi.nlm.nih.gov/clinvar/RCV000041374.5
1,156106185,not_specified,other,0.2492,Het,http://www.ncbi.nlm.nih.gov/clinvar/RCV000041315.5
1,156107534,not_specified,other,0.2202,Het,http://www.ncbi.nlm.nih.gov/clinvar/RCV000041327.4
1,161009523,Hyperlipidemia\x2c_familial_combined\x2c_susceptibility_to,other,0.1703,Het,http://www.ncbi.nlm.nih.gov/clinvar/RCV000013088.1
...
```


