# Script for running GeMoMa on multiple genomes
# GeMoMa uses as input the annotation and genome of a or multiple reference species,
# and  the genome of the target species.
set -eux

# Set arguments
TAXON="$1"

export TAXON

# Define directories, files and output path

## Target genomes directory
GENOME_DIR="${HOME}/piezo1_pipeline/genomes/${TAXON}"

## Output directory
OUTPUT_DIR="${HOME}/piezo1_pipeline/GeMoMa/results_${TAXON}"
mkdir -p "$OUTPUT_DIR"

## Make manifest file for parallelization
MANIFEST_FILE="${HOME}/piezo1_pipeline/GeMoMa/gemoma_manifest_${TAXON}.tsv"

find "$GENOME_DIR" \( -name "*.fna" -o -name "*.fasta" \) | parallel '
    CLADE=$(basename $(dirname {}))
    REF_GENOME=$(ls "'${HOME}'/piezo1_pipeline/GeMoMa/reference/${TAXON}/${CLADE}"/*.fna 2>/dev/null)
    REF_ANNOT=$(ls "'${HOME}'/piezo1_pipeline/GeMoMa/reference/${TAXON}/${CLADE}"/*.gff 2>/dev/null)
    if [[ -n "$REF_GENOME" && -n "$REF_ANNOT" ]]; then
        echo -e "{}\t${REF_GENOME}\t${REF_ANNOT}\t${CLADE}"
    fi
' > "$MANIFEST_FILE"

run_gemoma() {
    local target="$1"
    local ref_g="$2"
    local ref_a="$3"
    local clade="$4"
    local out_base="$5"
    local dir_name=$(basename "$target" | sed 's/\.[^.]*$//')
    local log_file="${out_base}/${dir_name}_error.log"

    echo "Starting $dir_name..."

    GeMoMa GeMoMaPipeline threads=4 \
        t="$target" \
        outdir="${out_base}/${dir_name}_GeMoMa_results" \
        s=own i="$clade" a="$ref_a" g="$ref_g" \
        p=true \
        GeMoMa.m=50000 \
        AnnotationFinalizer.r=NO > "$log_file" 2>&1 || { echo "Job $dir_name failed. Check $log_file"; return 1; }
}
export -f run_gemoma

# Loop through each genome file in the directory. 
# ! Set high intron max length since Piezo1 regions is quite large and contains long introns
cat "$MANIFEST_FILE" | parallel -j 1 --colsep '\t' run_gemoma {1} {2} {3} {4} "$OUTPUT_DIR"