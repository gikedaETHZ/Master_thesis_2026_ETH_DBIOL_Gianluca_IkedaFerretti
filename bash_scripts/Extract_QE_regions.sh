# Based on an MSA of PIEZO1 sequences from which Q/E repeats regions have been visually identified,
# start and end positions were selected as the first position where the amino acid is consereved in at least
# around 90% of the sequences.
# The number of Q and E residues, the length of the region, the ration Q+E/length is calculated for each region
# plus a sum of these variable Q/E regions is also calculated
# Input: MSA file of PIEZO1 sequences and the identified interval of Q/E repeat in the MSA.
# Interval file format: "start:end" since it uses seqkit to extract the region from the MSA.

# Inputs
MSA_FILE="$1"
INTERVAL_FILE="$2" # One region per line e.g: 930:970:QE_region1

{
    echo -e "species\tregion\tQ_count\tE_count\tQE_sum\tregion_length\tE_ratio\tNeg_charge\tNet_charge\tNet_charge_no_H"

    # Loop using the intervals from the interval file and gather sequences with seqkit
    # then use awk to count desired values mentioned in the beginning
    cat "$INTERVAL_FILE" | while IFS=':' read start end label; do
    seqkit subseq -r ${start}:${end} "$MSA_FILE" | seqkit seq -g | seqkit fx2tab | \
        awk -F'\t' -v reg="$label" '{
        name = $1; seq = $2;
        match(name, /\[([^\]]+)\]/, arr);
        species = arr[1];
        sub(/ \(.*/, "", species);
        len = length(seq);
        q = gsub(/Q/, "Q", seq);
        e = gsub(/E/, "E", seq);
        d = gsub(/D/, "D", seq);
        h = gsub(/H/, "H", seq);
        r = gsub(/R/, "R", seq);
        k = gsub(/K/, "K", seq);
        neg_charge = e + d;
        pos_charge = h + r + k;
        pos_charge_no_h = r + k;
        net_charge = pos_charge - neg_charge;
        net_charge_no_h = pos_charge_no_h - neg_charge;
        qe = q + e;
        e_ratio = e/len;
        qe_ratio = qe/len;
        printf "%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%d\n",
        species, reg, q, e, qe, len, e_ratio, neg_charge, net_charge, net_charge_no_h;
        }'
    done
} > QE_regions.tsv

sed -i 's/ /_/g' QE_regions.tsv




