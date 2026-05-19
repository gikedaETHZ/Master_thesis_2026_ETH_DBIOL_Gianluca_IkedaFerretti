# Pipeline for Quality Assessment of Protein Sequences against a Reference
# Uses pident, qcov, and length ratio
set -eux

# Usage to show help message
usage() {
    echo "Usage: $0 <input_fasta> <reference_accession>"
    echo ""
    echo "Run quality assessment of protein sequences against a selected reference with BLASTP."
    echo ""
    echo "Arguments:"
    echo " input_fasta:             Multi-FASTA file containing protein sequences to be assessed, including the reference sequence."
    echo " reference_accession:     Reference protein accession ID (e.g., NP_001136336.2 Human PIEZO1 for primates)"
    echo " output_dir:              Optional output directory for results (default: ./<input_fasta>_QA_results)"
    exit 1
}   

if [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
    usage
fi

if [ "$#" -lt 2 ]; then
    usage
fi

# Input arguments
SEQS_FILE="$1"

REF_SPECIES="$2"

OUTPUT_DIR="${3:-./$(basename "${SEQS_FILE}")_QA_results}"
mkdir -p "${OUTPUT_DIR}"

# Make temporary directory for intermediate files
TEMP_DIR=$(mktemp -d)
trap "rm -rf ${TEMP_DIR}" EXIT

# Threshold file format:
# THRESHOLD  MIN_ID  MIN_QUERY_COV  MIN_LEN_RATIO  MAX_LEN_RATIO
cat << EOF > "${OUTPUT_DIR}/quality_thresholds.tsv"
LIGHT           70  85  90  110
INTERMEDIATE    75  90  95  105
STRICT          80  95  98  102
EOF


# Getting all species from headers
extract_species() {
    local header="$1"
    if [[ "${header}" =~ \[([^]]+)\] ]]; then
        echo "${BASH_REMATCH[1]}" | awk '{print $1"_"$2}'
    fi
}

while IFS= read -r line; do
    if [[ "${line}" =~ ^\> ]]; then
        header="${line#>}"
        seqid=$(echo "${header}" | awk '{print $1}')
        species=$(extract_species "${header}")
        echo -e "${seqid}\t${species}\t${header}"
    fi
done < "${SEQS_FILE}" > "${TEMP_DIR}/seq_species_list.tsv"

# Get reference ID and Name
REF_ID=$(grep -i "${REF_SPECIES}" "${SEQS_FILE}" | head -n 1 | sed 's/^>//' | cut -d' ' -f1)
REF_NAME=$(grep -w "${REF_ID}" "${TEMP_DIR}/seq_species_list.tsv" | cut -f2)

# Get reference and query sequences into separate temporary files
seqkit grep -p "${REF_ID}" "${SEQS_FILE}" > "${TEMP_DIR}/reference.fasta" 


# BLASTP to given reference
makeblastdb \
    -in "${TEMP_DIR}/reference.fasta" \
    -dbtype prot \
    -out "${TEMP_DIR}/ref_db" > /dev/null

blastp \
    -query "${SEQS_FILE}" \
    -db "${TEMP_DIR}/ref_db" \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -evalue 1e-20 \
    -outfmt "6 qseqid qlen pident qcovs length slen bitscore" \
    > "${OUTPUT_DIR}/blastp_raw_results.tsv"

# Starting sensitivity summary file
echo -e "Threshold\tSpecies_with_Sequences\tSpecies_Missing\tSpecies_Total" \
> "${OUTPUT_DIR}/sensitivity_summary.tsv"


# Looping for each threshold set
while read -r THRESHOLD MIN_ID MIN_QUERY_COV MIN_LEN_RATIO MAX_LEN_RATIO; do

    THRESH_DIR="${OUTPUT_DIR}/${THRESHOLD}_results"
    mkdir -p "${THRESH_DIR}"

    PASSED_QA="${TEMP_DIR}/passed_qa.tsv"

    echo -e "Sequence_ID\tSpecies\tIdentity\tQuery_Cov\tLen_Ratio\tStatus\tReason" \
        > "${THRESH_DIR}/qa_results.tsv"

    > "${PASSED_QA}"

    while IFS=$'\t' read -r qseqid qlen pident qcovs alen slen bitscore; do

        species=$(grep -w "^${qseqid}" "${TEMP_DIR}/seq_species_list.tsv" | cut -f2)
        len_ratio=$(awk -v q="${qlen}" -v s="${slen}" 'BEGIN{printf "%.1f", q/s*100}')
        
        reason=""

        m_start=$(seqkit grep -p "${qseqid}" "${SEQS_FILE}" | seqkit seq -s | head -c1)
        if [[ "${m_start}" != "M" ]]; then
            reason+="non_canonical_start;"
        fi

        (( $(awk -v x="${pident}" -v y="${MIN_ID}" 'BEGIN{print x<y}') )) && reason+="low_identity;"
        (( $(awk -v x="${qcovs}" -v y="${MIN_QUERY_COV}" 'BEGIN{print x<y}') )) && reason+="low_query_coverage;"
        (( $(awk -v x="${len_ratio}" -v y="${MIN_LEN_RATIO}" 'BEGIN{print x<y}') )) && reason+="too_short;"
        (( $(awk -v x="${len_ratio}" -v y="${MAX_LEN_RATIO}" 'BEGIN{print x>y}') )) && reason+="too_long;"



        if [[ -z "${reason}" ]]; then
            status="PASS"
            echo -e "${qseqid}\t${species}\t${qlen}\t${pident}\t${qcovs}\t${len_ratio}" >> "${PASSED_QA}"
        else
            status="FAIL"
        fi

        echo -e "${qseqid}\t${species}\t${pident}\t${qcovs}\t${len_ratio}\t${status}\t${reason%%;}" \
            >> "${THRESH_DIR}/qa_results.tsv"

    done < "${OUTPUT_DIR}/blastp_raw_results.tsv"

    # Selecting best isoform per species from passed sequences
    cut -f2 "${TEMP_DIR}/passed_qa.tsv" | sort -u | while read -r species; do
        grep -w "${species}" "${TEMP_DIR}/passed_qa.tsv" > "${TEMP_DIR}/sp.tsv"

        best_id=""
        best_score=-1

        # scoring formula: score = (pident * 10) + (qcovs * 10) - lenght_penalty
        while IFS=$'\t' read -r qseqid sp qlen pident qcovs len_ratio; do
            len_penalty=$(awk -v x="${len_ratio}" 'BEGIN{print int(sqrt((x-100)^2)*10)}')
            score=$(( $(awk -v x="${qcovs}" 'BEGIN{print int(x*10)}') + $(awk -v x="${pident}" 'BEGIN{print int(x*10)}') - len_penalty ))

            if [[ "${score}" -gt "${best_score}" ]]; then
                best_score="${score}"
                best_id="${qseqid}"

            fi
        done < "${TEMP_DIR}/sp.tsv"

        echo "${best_id}"
    done  > "${TEMP_DIR}/good_ids.txt"

    # Generate FASTA file with selected sequences
    seqkit grep -f "${TEMP_DIR}/good_ids.txt" "${SEQS_FILE}" > "${THRESH_DIR}/final_selected_sequences.fasta"

    # Species with sequences
    cut -f2 "${TEMP_DIR}/passed_qa.tsv" | sort -u > "${TEMP_DIR}/good_species.txt"

    # All species in dataset
    cut -f2 "${TEMP_DIR}/seq_species_list.tsv" | grep -v "${REF_NAME}" | sort -u > "${TEMP_DIR}/all_species.txt"

    # Species without sequences (set difference)
    comm -23 "${TEMP_DIR}/all_species.txt" "${TEMP_DIR}/good_species.txt" > "${THRESH_DIR}/failed_species.txt"

    # Sensitivity counts
    SP_WITH_SEQS=$(wc -l < "${TEMP_DIR}/good_species.txt")
    SP_MISS=$(wc -l < "${THRESH_DIR}/failed_species.txt")
    TOTAL_SP=$(wc -l < "${TEMP_DIR}/all_species.txt")

    echo -e "${THRESHOLD}\t${SP_WITH_SEQS}\t${SP_MISS}\t${TOTAL_SP}" >> "${OUTPUT_DIR}/sensitivity_summary.tsv"

done < "${OUTPUT_DIR}/quality_thresholds.tsv"