# Run all-vs-all BLASTP for selected sequences and generate matrix file for heatmap visualization of pairwise identities

set -eux

# Input arguments
SEQS_FILE="$1"

# Check consistency between phylogenetic groups

OUTPUT_DIR=~/piezo1_pipeline/all_vs_all_blastp/$(basename "${SEQS_FILE}" .fasta)_results
mkdir -p "${OUTPUT_DIR}"

TEMP_DIR=$(mktemp -d)
trap "rm -rf ${TEMP_DIR}" EXIT

# Extract species from headers
extract_species() {
    local header="$1"
    if [[ "${header}" =~ \[([^]]+)\] ]]; then
        echo "${BASH_REMATCH[1]}" | awk '{print $1"_"$2}'
    fi
}

while IFS= read -r line; do
    if [[ "${line}" =~ ^\> ]]; then
        header="${line#>}"
        seqid=$(echo "${header}" | cut -d' ' -f1)
        species=$(extract_species "${header}")
        echo -e "${seqid}\t${species}\t${header}"
    fi
done < "${SEQS_FILE}" > "${TEMP_DIR}/seq_species_list.tsv"


# All-vs-all BLASTP
makeblastdb \
    -in "${SEQS_FILE}" \
    -dbtype prot \
    -out "${TEMP_DIR}/all_db" > /dev/null

blastp \
    -query "${SEQS_FILE}" \
    -db "${TEMP_DIR}/all_db" \
    -max_target_seqs 100 \
    -max_hsps 1 \
    -evalue 1e-20 \
    -outfmt "6 qseqid sseqid pident qcovs length qlen slen gaps bitscore" \
    > "${OUTPUT_DIR}/all_vs_all_blast.tsv"

# Build pairwise identity tsv file
# Map sequence IDs to species names, then create species-level matrix
echo -e "Query_Species\tSubject_Species\tQuery_ID\tSubject_ID\tIdentity\tQuery_Cov\tGaps" \
    > "${OUTPUT_DIR}/pairwise_identity.tsv"

while IFS=$'\t' read -r qseqid sseqid pident qcovs alen qlen slen gaps bitscore; do
    # Skip self-hits
    [[ "${qseqid}" == "${sseqid}" ]] && continue

    qspecies=$(grep -w "^${qseqid}" "${TEMP_DIR}/seq_species_list.tsv" | cut -f2)
    sspecies=$(grep -w "^${sseqid}" "${TEMP_DIR}/seq_species_list.tsv" | cut -f2)

    echo -e "${qspecies}\t${sspecies}\t${qseqid}\t${sseqid}\t${pident}\t${qcovs}\t${gaps}"
done < "${OUTPUT_DIR}/all_vs_all_blast.tsv" \
    >> "${OUTPUT_DIR}/pairwise_identity.tsv"

# Species-level identity matrix (for heatmap visualization)
# For each species pair, take the average of reciprocal identities
tail -n +2 "${OUTPUT_DIR}/pairwise_identity.tsv" \
    | awk -F'\t' '{
        if ($1 < $2) key = $1"\t"$2;
        else key = $2"\t"$1;
        sum[key] += $5;
        count[key]++;
    }
    END {
        print "Species_1\tSpecies_2\tMean_Identity\tN_comparisons";
        for (k in sum) {
            printf "%s\t%.2f\t%d\n", k, sum[k]/count[k], count[k];
        }
    }' | sort -t$'\t' -k3 -n \
    > "${OUTPUT_DIR}/species_pairwise_matrix.tsv"

# Clean up associative array
unset sum count
echo ""
echo "DONE. Results in  ${OUTPUT_DIR}"