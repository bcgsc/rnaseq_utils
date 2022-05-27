set -euo pipefail

# A script for extracting ground truth transcript IDs from a Trans-NanoSim simulated reads FASTQ file
# written by Ka Ming Nip @kmnip

cat=cat

if [[ ${1} == *.gz ]]
then
    cat=zcat
fi

${cat} ${1} | awk 'NR % 4 == 1 {print}' | sed -e 's/^@//g' -e 's/_.*//g' | sort | uniq
