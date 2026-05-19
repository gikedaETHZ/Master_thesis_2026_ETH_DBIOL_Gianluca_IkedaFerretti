#!/usr/bin/env python
# coding: utf-8

# In[1]:


import requests # allows to send HTTP requests easily (querying online databases)
import re # allows for regular expressions usage and many functions related to it
import os # allows to use operating system dependent functions
import argparse # allows to set variables when calling from command line
from Bio import SeqIO 
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import xml.etree.ElementTree as ET
import gzip
import time
import csv

# Define the taxonomic group to look for in databases and a reference length for that group.

def parse_args():
    parser = argparse.ArgumentParser(description="Collect PIEZO1 sequences from UniProt, NCBI, and DNAZoo")
    parser.add_argument("--taxon", required=True, help="Taxonomic group (e.g., primates, rodentia)")
    parser.add_argument("--gene", default="PIEZO1", help="Gene name (default:PIEZO1)")
    parser.add_argument("--output", default=None, help="Output FASTA filename")
    parser.add_argument("--csv", action="store_true", help="Also generate summary CSV")
    return parser.parse_args()


# In[2]:


# Make species name all the same just to be sure
def normalize_organism_name(organism):
    parts = organism.split()
    if len(parts) >= 2:
        return f"{parts[0].lower()} {parts[1].lower()}"
    return organism.lower()


# In[3]:


# Extracting sequences from UniProt
def collect_uniprot_sequences(
    gene_name,
    protein_name,
    taxonomy,
    annotation_score_min=2,
):

    base_url = "https://rest.uniprot.org/uniprotkb/search"

    query = f'gene:{gene_name} AND protein_name:Piezo type mechanosensitive ion channel component 1 AND taxonomy_name:{taxonomy} AND existence:3'

    params = {
        'query': query,
        'format': 'tsv',
        'fields': 'accession,gene_names,organism_name,annotation_score,protein_existence,sequence',
        'size': 500
    }

    print(f"Querying UniProt with: {query}")

    response = requests.get(base_url, params=params)

    # Print errors related to the request for troubleshooting 
    if response.status_code != 200:
        print(f"Error: {response.status_code}")
        print(response.text)
        return []

    # Prepares tsv for parsing
    lines = response.text.strip().split('\n')

    # Check if there are any
    if len(lines) <= 1:
        print("No sequences fond in UniProt")
        return []

    # Save header of sequence
    header = lines[0].split('\t')

    records = []
    for line in lines[1:]:
        fields = line.split('\t')
        row = dict(zip(header, fields))

        accession = row.get('Entry', 'unknown')
        organism = row.get('Organism', 'unknown')
        sequence = row.get('Sequence', '')

        try:
            score = float(row.get('Annotation',0))
        except ValueError:
            score = 0

        if score < annotation_score_min:
            continue

        if sequence:
            rec = SeqRecord(
                Seq(sequence),
                id=f"UniProt|{accession}",
                description=f"[{organism}]"
            )
            records.append(rec)

    print(f"UniProt: Retrieved {len(records)} sequences with annotation score >= {annotation_score_min}")

    return records


# In[4]:


# Extracting sequences from NCBI sequences
def collect_NCBI_sequences(gene, taxonomy, retmax=1000):

    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params={
        'db':'protein',
        'term': f'{gene}[gene] AND {taxonomy}[organism] AND refseq[filter]',
        'retmax':retmax,
        'retmode':'xml'
    }

    response = requests.get(search_url, params=params)

    if response.status_code != 200:
        print(f'Error: {response.status_code}')
        return []

    root = ET.fromstring(response.text)
    ids = [id_elem.text for id_elem in root.findall('.//Id')]


    print(f'Found {len(ids)} protein IDs')

    if not ids:
        return []

    fetch_url =  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    all_records= []

    batch_size = 50

    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]

        params = {
            'db':'protein',
            'id': ','.join(batch),
            'rettype':'fasta',
            'retmode':'text'
        }

        response = requests.get(fetch_url, params=params)

        if response.status_code == 200:
            records = list(SeqIO.parse(StringIO(response.text), 'fasta'))

            # Tag records with NCBI source
            for rec in records:
                rec.id = f"NCBI|{rec.id}"
            all_records.extend(records)
            print(f'Fetched batch {i // batch_size+1}: {len(records)} sequences')
        else:
            print(f'Error: {response.status_code}')

        time.sleep(0.35)

    print(f"NCBI: Retrieved {len(all_records)} sequences (all, no filtering)")
    return all_records


# In[5]:


# Extracting sequences from DNAZoo 

## Generate a list with all species for a certain taxon. For all species I mean all those available in NCBI.

def get_all_species_names(taxon, retmax=15000, filter_unidentified=True):

    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {
        'db':'taxonomy',
        'term':f'{taxon}[Subtree] AND species[Rank]',
        'retmax':retmax,
        'retmode':'xml'
    }

    response = requests.get(search_url, params=params)

    if response.status_code != 200:
        print(f'Error {response.status_code}')
        return []

    root = ET.fromstring(response.content)

    tax_ids = [id_elem.text for id_elem in root.findall('.//Id')]

    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    species_list = []
    batch_size = 200

    for i in range(0, len(tax_ids), batch_size):
        batch = tax_ids[i:i+batch_size]
        print(f'Fetching batch {i+1}')

        params = {
            'db':'taxonomy',
            'id':','.join(batch),
            'retmode':'xml' 
        }

        response = requests.get(fetch_url, params=params)

        if response.status_code != 200:
            print(f'Error {response.status_code}')
            return []

        root_fetch = ET.fromstring(response.content)

        for taxon_elem in root_fetch.findall('./Taxon'):
            sci_name_elem = taxon_elem.find('ScientificName')
            if sci_name_elem is not None:
                name = sci_name_elem.text.replace(' ', '_')

                if filter_unidentified:
                    if '_sp.' in name:
                        continue
                    if '_x_' in name:
                        continue
                    if '_cf._' in name:
                        continue

                species_list.append(name)

        time.sleep(0.35)

    print(f'Total species {len(species_list)}')
    return sorted(species_list)

## Find protein file URL for a species in DNAZoo

def get_dnazoo_protein_file_url(species_list):

    url = "https://dnazoo.s3.wasabisys.com/"
    params = {'prefix':f'{species_list}/'}

    response = requests.get(url, params=params) 
    if response.status_code != 200:
        return []

    root = ET.fromstring(response.content)
    namespace = {'s3': 'http://s3.amazonaws.com/doc/2006-03-01/'}

    for contents in root.findall('.//s3:Contents', namespace):
        key = contents.find('s3:Key', namespace)
        if key is not None and 'HiC.fasta_v2.functional.proteins.fasta.gz' in key.text:
            return f"https://dnazoo.s3.wasabisys.com/{key.text}"

    return None

## Download zip files of functional proteins fasta files from DNAZoo of available species and create a directory with them

def collect_DNAZoo_sequences(GENE, TAXON, base_dir="./dnazoo_downloads"):

    output_dir = os.path.join(base_dir, TAXON)

    os.makedirs(output_dir, exist_ok=True)

    species_list = get_all_species_names(TAXON)

    all_records = []
    downloaded_species = []

    for i, species in enumerate(species_list):
        print(f"[{i+1}/{len(species_list)}] {species}...", end=" ")

        url = get_dnazoo_protein_file_url(species)
        if url is None:
            print("No data")
            continue

        print("↓", end=" ")

        try:
            response = requests.get(url, stream=True, timeout=60)
            filename = os.path.join(output_dir, f"{species}_functional_proteins.fasta.gz")

            with open(filename, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            piezo1_count = 0
            with gzip.open(filename, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    if GENE in record.description.upper():
                        original_id = record.id
                        record.id = f"DNAZoo|{species}|{original_id}"
                        record.description = f"[{species.replace('_',' ')}]"
                        all_records.append(record)
                        piezo1_count += 1

            if piezo1_count > 0:
                print(f"Downloaded ({piezo1_count} PIEZO1)")
                downloaded_species.append(species)
            else:
                print("No piezo1 found")

            time.sleep(0.5)

        except Exception as e:
            print(f"Error: {e})")
            continue

    print(f"\nDNAZoo: Retrieved {len(all_records)} sequences from {len(downloaded_species)} species (all, no filtering)")
    return all_records


# In[6]:


# Gather all sequences and create a fasta file with all of those

def collect_all_piezo1_sequences(GENE, TAXON, output_file=None):
    print("=" * 80)
    print(f"Collecting all available PIEZO1 sequences for: {TAXON}")
    print("=" * 80)

    all_sequences = []

    # Run UniProt function
    print("--- UniProt ---")
    uniprot_seqs = collect_uniprot_sequences(
        gene_name=GENE,
        protein_name=GENE,
        taxonomy=TAXON
    )

    all_sequences.extend(uniprot_seqs)

    # Run NCBI function
    print("--- NCBI ---")
    ncbi_seqs = collect_NCBI_sequences(
        gene=GENE,
        taxonomy=TAXON
    )
    all_sequences.extend(ncbi_seqs)

    # Run DNAZoo
    print("--- DNAZoo ---")
    dnazoo_seqs = collect_DNAZoo_sequences(GENE, TAXON)
    all_sequences.extend(dnazoo_seqs)

    # Quick Summary
    print("=" * 80)
    print("COLLECTION SUMMARY")
    print("=" * 80)
    print(f"UniProt: {len(uniprot_seqs)} sequences")
    print(f"NCBI: {len(ncbi_seqs)} sequences")
    print(f"DNAZoo: {len(dnazoo_seqs)} sequences")
    print("=" * 80)
    print(f"TOTAL {len(all_sequences)} sequences")

    if output_file is None:
        output_file = f"piezo1_{TAXON}_available_sequences.fasta"

    with open(output_file, 'w') as f:
        SeqIO.write(sorted(all_sequences, key=lambda r: r.description), f, 'fasta')

    print(f"\nSaved to: {output_file}")

    return all_sequences


# In[7]:


if __name__ == "__main__":
    args = parse_args()

    output_file = args.output or f"piezo1_{args.taxon}_available_sequences.fasta"

    collect_all_piezo1_sequences(args.gene, args.taxon, output_file=output_file)

    if args.csv:
        csv_file=f"list_{args.taxon}_all_available_piezo1.csv"

        with open(output_file) as handle:
            piezo1_list = []

            for piezo1 in SeqIO.parse(handle, 'fasta'):
                description = piezo1.description
                accession = piezo1.id

                match = re.search(r'\[([^\]]+)\]', description)
                db_query = re.search(r'_DNAZoo', description)

                if match:
                    name = match.group(1).strip()
                    database = 'DNAZoo' if db_query else 'NCBI'
                else:
                    match = re.search(r'OS=([^(]+)', description)
                    if match:
                        name = match.group(1).strip()
                        database = 'UniProt'
                piezo1_list.append({
                    'species': name,
                    'accession': accession,
                    'database': database
                })

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['species', 'accession', 'database'])
            writer.writeheader()
            writer.writerows(piezo1_list)

        print(f"Summary CSV saved to: {csv_file}")

