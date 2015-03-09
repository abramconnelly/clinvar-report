#!/bin/bash
#
#
descr="ClinVar Report"

pt="clinvar-report.pipeline"

basedir="crunch_scripts"

clinvar="2423f05145f28120f8ff1bcff260932f+67/clinvar.vcf"
vcf="71509d7fdb9f52b389b2efd84815be2f+63/abe.vcf"

if [ "$clinvar" == "" ] || [ "$vcf" == "" ]
then
  echo "provide clinvar and vcf files"
  exit 1
fi

template=` cat "$pt" | jq ".components[].repository=\"$USER\"" `

opt=""
opt=" $opt ClinVarReport::CLINVAR=$clinvar"
opt=" $opt ClinVarReport::VCF=$vcf"

echo "$template" | jq .


echo arv-run-pipeline-instance --description "$descr" --submit --template <( echo "$template" ) $opt
arv-run-pipeline-instance --description "$descr" --submit --template <( echo "$template" ) $opt
