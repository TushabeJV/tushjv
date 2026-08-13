#!/usr/bin/env python3

import subprocess
import csv
import os
import sys
import re
from Bio import SeqIO


# ============================================================
# INPUT FILES
# ============================================================

# Query contigs generated from the pit latrine
QUERY_FASTA = "mansoni_sequences.fasta"

# WBPS19 S. mansoni genome
REFERENCE_FASTA = "schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa"

# WBPS19 annotation
GFF3_FILE = "schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"


# ============================================================
# OUTPUT FILES expected
# ============================================================

BLAST_DB = "smansoni_WBPS19_genome"

BLAST_OUTPUT = "mansoni_contig_blast.tsv"

FINAL_OUTPUT = "mansoni_contig_gene_ids.csv"


# ============================================================
# BLAST PARAMETERS
# ============================================================

MIN_IDENTITY = 80.0
MIN_COVERAGE = 70.0
MAX_EVALUE = 1e-10

THREADS = 4


# ============================================================
# Checking Input files
# ============================================================

print("\nChecking input files...")

for filename in [
    QUERY_FASTA,
    REFERENCE_FASTA,
    GFF3_FILE
]:

    if not os.path.exists(filename):

        sys.exit(
            f"\nERROR: File not found:\n{filename}\n"
        )

    print(f"  OK: {filename}")


# ============================================================
# READ QUERY CONTIGS
# ============================================================

print("\nReading query FASTA...")

contigs = list(
    SeqIO.parse(QUERY_FASTA, "fasta")
)

print(
    f"Number of query contigs: {len(contigs)}"
)


# ============================================================
# Creating a  BLAST database
# ============================================================

print("\nCreating BLAST database...")

subprocess.run(
    [
        "makeblastdb",
        "-in",
        REFERENCE_FASTA,
        "-dbtype",
        "nucl",
        "-out",
        BLAST_DB
    ],
    check=True
)


# ============================================================
# Running BLASTN
# ============================================================

print("\nRunning BLASTN...")

blast_fields = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen"
]

subprocess.run(
    [
        "blastn",
        "-query",
        QUERY_FASTA,
        "-db",
        BLAST_DB,
        "-out",
        BLAST_OUTPUT,
        "-outfmt",
        "6 " + " ".join(blast_fields),
        "-evalue",
        "1e-10",
        "-max_target_seqs",
        "20",
        "-num_threads",
        str(THREADS)
    ],
    check=True
)


# ============================================================
# Choosing the best BLAST hits
# ============================================================

print("\nSelecting best BLAST hit for each contig...")

best_hits = {}


with open(BLAST_OUTPUT) as handle:

    for line in handle:

        if not line.strip():
            continue

        fields = line.rstrip().split("\t")

        (
            qseqid,
            sseqid,
            pident,
            length,
            mismatch,
            gapopen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
            qlen,
            slen
        ) = fields

        pident = float(pident)
        length = int(length)
        evalue = float(evalue)
        bitscore = float(bitscore)
        qlen = int(qlen)

        # Calculate query coverage
        coverage = (
                           length / qlen
                   ) * 100

        # Apply filters
        if pident < MIN_IDENTITY:
            continue

        if coverage < MIN_COVERAGE:
            continue

        if evalue > MAX_EVALUE:
            continue

        start = min(
            int(sstart),
            int(send)
        )

        end = max(
            int(sstart),
            int(send)
        )

        strand = (
            "+"
            if int(sstart) <= int(send)
            else "-"
        )

        hit = {
            "reference": sseqid,
            "start": start,
            "end": end,
            "strand": strand,
            "identity": pident,
            "coverage": coverage,
            "evalue": evalue,
            "bitscore": bitscore
        }

        # Keep highest bitscore
        if (
                qseqid not in best_hits
                or bitscore >
                best_hits[qseqid]["bitscore"]
        ):

            best_hits[qseqid] = hit


print(
    f"Contigs with acceptable BLAST hits: "
    f"{len(best_hits)}"
)

# ============================================================
# Parse GFF3 attributes to get sm gene IDs
# ============================================================

def parse_attributes(attribute_string):

    attributes = {}

    for item in attribute_string.split(";"):

        item = item.strip()

        if not item:
            continue

        if "=" in item:

            key, value = item.split(
                "=",
                1
            )

            attributes[key] = value

        elif " " in item:

            key, value = item.split(
                None,
                1
            )

            attributes[key] = value.strip('"')

    return attributes


# ============================================================
# Read in the gff3 file
# ============================================================

print("\nReading GFF3 annotation...")

genes = []


with open(
        GFF3_FILE,
        encoding="utf-8"
) as gff:

    for line in gff:

        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")

        if len(fields) != 9:
            continue

        seqid = fields[0]
        feature_type = fields[2]

        start = int(fields[3])
        end = int(fields[4])

        strand = fields[6]

        attributes = parse_attributes(
            fields[8]
        )

        # We only need gene features
        if feature_type.lower() != "gene":
            continue

        gene_id = (
                attributes.get("ID")
                or attributes.get("gene_id")
                or attributes.get("gene")
        )

        if not gene_id:
            continue

        genes.append(
            {
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
                "gene_id": gene_id
            }
        )


print(
    f"Number of gene annotations loaded: "
    f"{len(genes)}"
)


# ============================================================
# FIND GENE OVERLAPPING BLAST HIT
# ============================================================

def find_overlapping_genes(
        chromosome,
        start,
        end
):

    matches = []

    for gene in genes:

        if gene["seqid"] != chromosome:
            continue

        # Check overlap
        if (
                start <= gene["end"]
                and end >= gene["start"]
        ):

            overlap_start = max(
                start,
                gene["start"]
            )

            overlap_end = min(
                end,
                gene["end"]
            )

            overlap_length = (
                    overlap_end
                    - overlap_start
                    + 1
            )

            matches.append(
                (
                    gene,
                    overlap_length
                )
            )

    # Sort by largest overlap
    matches.sort(
        key=lambda x: x[1],
        reverse=True
    )

    return matches


# ============================================================
# CLEAN GENE ID
# ============================================================

def clean_gene_id(gene_id):

    """
    Extract standard Smp_XXXXXX ID.

    Example:

        gene:Smp_123450
        Smp_123450
        Smp_123450.1

    becomes:

        Smp_123450
    """

    match = re.search(
        r"(Smp_\d+)",
        gene_id
    )

    if match:

        return match.group(1)

    return gene_id


# ============================================================
# GENERATE FINAL RESULTS
# ============================================================

print("\nMapping BLAST hits to GFF3 genes...")

results = []


for record in contigs:

    contig_id = record.id

    if contig_id not in best_hits:

        results.append(
            {
                "contig_id": contig_id,
                "gene_id": "NO_BLAST_HIT",
                "reference": "",
                "start": "",
                "end": "",
                "strand": "",
                "identity": "",
                "coverage": "",
                "evalue": "",
                "bitscore": "",
                "gene_overlap": ""
            }
        )

        continue


    hit = best_hits[contig_id]


    overlapping_genes = (
        find_overlapping_genes(
            hit["reference"],
            hit["start"],
            hit["end"]
        )
    )


    if overlapping_genes:

        gene, overlap_length = (
            overlapping_genes[0]
        )

        gene_id = clean_gene_id(
            gene["gene_id"]
        )

        results.append(
            {
                "contig_id": contig_id,
                "gene_id": gene_id,
                "reference": hit["reference"],
                "start": hit["start"],
                "end": hit["end"],
                "strand": hit["strand"],
                "identity": round(
                    hit["identity"], 2
                ),
                "coverage": round(
                    hit["coverage"], 2
                ),
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "gene_overlap": overlap_length
            }
        )

    else:

        results.append(
            {
                "contig_id": contig_id,
                "gene_id": "NO_GENE_OVERLAP",
                "reference": hit["reference"],
                "start": hit["start"],
                "end": hit["end"],
                "strand": hit["strand"],
                "identity": round(
                    hit["identity"], 2
                ),
                "coverage": round(
                    hit["coverage"], 2
                ),
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "gene_overlap": 0
            }
        )


# ============================================================
# WRITE CSV
# ============================================================

print("\nWriting final output...")

with open(
        FINAL_OUTPUT,
        "w",
        newline=""
) as output:

    writer = csv.DictWriter(
        output,
        fieldnames=[
            "contig_id",
            "gene_id",
            "reference",
            "start",
            "end",
            "strand",
            "identity",
            "coverage",
            "evalue",
            "bitscore",
            "gene_overlap"
        ]
    )

    writer.writeheader()

    writer.writerows(results)


# ============================================================
# SUMMARY
# ============================================================

total = len(results)

identified = sum(
    1
    for r in results
    if r["gene_id"].startswith("Smp_")
)

no_blast = sum(
    1
    for r in results
    if r["gene_id"] == "NO_BLAST_HIT"
)

no_gene = sum(
    1
    for r in results
    if r["gene_id"] == "NO_GENE_OVERLAP"
)


print("\n============================================")
print("          ANALYSIS COMPLETE")
print("============================================")

print(f"Total contigs:          {total}")
print(f"Gene IDs identified:    {identified}")
print(f"No BLAST hit:           {no_blast}")
print(f"No gene overlap:        {no_gene}")

print("\nFinal output:")
print(f"  {FINAL_OUTPUT}")

print("============================================")