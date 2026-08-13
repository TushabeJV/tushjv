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

QUERY_FASTA = "mansoni_sequences.fasta"

REFERENCE_FASTA = (
    "schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa"
)

GFF3_FILE = (
    "schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
)


# ============================================================
# OUTPUT FILES
# ============================================================

BLAST_DB = "smansoni_WBPS19_genome"

BLAST_OUTPUT = "mansoni_contig_blast.tsv"

FINAL_OUTPUT = "mansoni_contig_gene_ids.csv"


# ============================================================
# PARAMETERS
# ============================================================

MIN_IDENTITY = 80.0
MIN_COVERAGE = 70.0
MAX_EVALUE = 1e-10

THREADS = 4

# Search this distance around a BLAST hit if no gene overlaps
NEARBY_DISTANCE = 10000


# ============================================================
# CHECK INPUT FILES
# ============================================================

print("\nChecking input files...")

for filename in [
    QUERY_FASTA,
    REFERENCE_FASTA,
    GFF3_FILE
]:

    if not os.path.isfile(filename):

        sys.exit(
            f"\nERROR: Cannot find:\n{filename}\n"
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
    f"Query contigs: {len(contigs)}"
)


# ============================================================
# PARSE GFF3 GENE FEATURES
# ============================================================

print("\nReading GFF3 gene annotations...")

genes_by_seq = {}


def parse_gff_attributes(attribute_string):

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

    return attributes


with open(
    GFF3_FILE,
    "r",
    encoding="utf-8"
) as gff:

    for line in gff:

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        seqid = fields[0]
        feature_type = fields[2]

        if feature_type != "gene":
            continue

        start = int(fields[3])
        end = int(fields[4])

        strand = fields[6]

        attributes = parse_gff_attributes(
            fields[8]
        )

        gene_identifier = attributes.get(
            "ID",
            ""
        )

        # Expected:
        # gene:Smp_000020

        match = re.search(
            r"(Smp_\d+)",
            gene_identifier
        )

        if match:

            gene_id = match.group(1)

        else:

            # Try Name if ID did not work
            name = attributes.get(
                "Name",
                ""
            )

            match = re.search(
                r"(Smp_\d+)",
                name
            )

            if match:

                gene_id = match.group(1)

            else:

                continue

        gene_record = {
            "gene_id": gene_id,
            "seqid": seqid,
            "start": start,
            "end": end,
            "strand": strand
        }

        genes_by_seq.setdefault(
            seqid,
            []
        ).append(
            gene_record
        )


# Sort genes by genomic position

for seqid in genes_by_seq:

    genes_by_seq[seqid].sort(
        key=lambda x: x["start"]
    )


total_genes = sum(
    len(x)
    for x in genes_by_seq.values()
)

print(
    f"Genes loaded from GFF3: {total_genes}"
)


# ============================================================
# CREATE BLAST DATABASE
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
# RUN BLASTN
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
        str(MAX_EVALUE),
        "-max_target_seqs",
        "20",
        "-num_threads",
        str(THREADS)
    ],
    check=True
)


# ============================================================
# SELECT BEST BLAST HIT
# ============================================================

print("\nSelecting best BLAST hit per contig...")

best_hits = {}


with open(BLAST_OUTPUT) as blast:

    for line in blast:

        if not line.strip():
            continue

        fields = line.rstrip().split("\t")

        (
            qseqid,
            sseqid,
            pident,
            alignment_length,
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
        alignment_length = int(
            alignment_length
        )

        qstart = int(qstart)
        qend = int(qend)

        sstart = int(sstart)
        send = int(send)

        evalue = float(evalue)
        bitscore = float(bitscore)

        qlen = int(qlen)

        # ----------------------------------------------------
        # CORRECT QUERY COVERAGE
        # ----------------------------------------------------

        aligned_query_length = (
            abs(qend - qstart) + 1
        )

        coverage = (
            aligned_query_length
            / qlen
        ) * 100


        if pident < MIN_IDENTITY:
            continue

        if coverage < MIN_COVERAGE:
            continue

        if evalue > MAX_EVALUE:
            continue


        genomic_start = min(
            sstart,
            send
        )

        genomic_end = max(
            sstart,
            send
        )

        strand = (
            "+"
            if sstart <= send
            else "-"
        )


        hit = {
            "reference": sseqid,
            "start": genomic_start,
            "end": genomic_end,
            "strand": strand,
            "identity": pident,
            "coverage": coverage,
            "evalue": evalue,
            "bitscore": bitscore
        }


        if (
            qseqid not in best_hits
            or bitscore >
               best_hits[qseqid]["bitscore"]
        ):

            best_hits[qseqid] = hit


print(
    f"Contigs with BLAST hits: "
    f"{len(best_hits)}"
)


# ============================================================
# FIND OVERLAPPING GENES
# ============================================================

def find_overlapping_genes(
    chromosome,
    hit_start,
    hit_end
):

    matches = []

    if chromosome not in genes_by_seq:
        return matches

    for gene in genes_by_seq[chromosome]:

        # No overlap possible
        if gene["end"] < hit_start:
            continue

        # Since genes are sorted, stop when beyond hit
        if gene["start"] > hit_end:
            break

        overlap_start = max(
            hit_start,
            gene["start"]
        )

        overlap_end = min(
            hit_end,
            gene["end"]
        )

        if overlap_start <= overlap_end:

            overlap_length = (
                overlap_end
                - overlap_start
                + 1
            )

            matches.append(
                {
                    "gene": gene,
                    "overlap": overlap_length
                }
            )

    matches.sort(
        key=lambda x: x["overlap"],
        reverse=True
    )

    return matches


# ============================================================
# FIND NEARBY GENES
# ============================================================

def find_nearby_genes(
    chromosome,
    hit_start,
    hit_end
):

    nearby = []

    if chromosome not in genes_by_seq:
        return nearby

    search_start = (
        hit_start
        - NEARBY_DISTANCE
    )

    search_end = (
        hit_end
        + NEARBY_DISTANCE
    )

    for gene in genes_by_seq[chromosome]:

        if gene["end"] < search_start:
            continue

        if gene["start"] > search_end:
            break

        # Calculate distance
        if gene["end"] < hit_start:

            distance = (
                hit_start
                - gene["end"]
            )

        elif gene["start"] > hit_end:

            distance = (
                gene["start"]
                - hit_end
            )

        else:

            distance = 0

        nearby.append(
            {
                "gene": gene,
                "distance": distance
            }
        )

    nearby.sort(
        key=lambda x: x["distance"]
    )

    return nearby


# ============================================================
# MAP CONTIGS TO GENES
# ============================================================

print("\nMapping BLAST hits to genes...")

results = []


for record in contigs:

    contig_id = record.id


    # --------------------------------------------------------
    # NO BLAST HIT
    # --------------------------------------------------------

    if contig_id not in best_hits:

        results.append(
            {
                "contig_id": contig_id,
                "gene_id": "NO_BLAST_HIT",
                "assignment": "NO_BLAST_HIT",
                "reference": "",
                "start": "",
                "end": "",
                "strand": "",
                "identity": "",
                "coverage": "",
                "evalue": "",
                "bitscore": "",
                "gene_overlap": "",
                "nearest_gene": "",
                "nearest_gene_distance": ""
            }
        )

        continue


    hit = best_hits[contig_id]


    # --------------------------------------------------------
    # FIND DIRECT GENE OVERLAP
    # --------------------------------------------------------

    overlapping = find_overlapping_genes(
        hit["reference"],
        hit["start"],
        hit["end"]
    )


    if overlapping:

        best_gene = overlapping[0]

        gene = best_gene["gene"]

        results.append(
            {
                "contig_id": contig_id,
                "gene_id": gene["gene_id"],
                "assignment": "DIRECT_GENE_OVERLAP",
                "reference": hit["reference"],
                "start": hit["start"],
                "end": hit["end"],
                "strand": hit["strand"],
                "identity": round(
                    hit["identity"],
                    2
                ),
                "coverage": round(
                    hit["coverage"],
                    2
                ),
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "gene_overlap": best_gene[
                    "overlap"
                ],
                "nearest_gene": "",
                "nearest_gene_distance": ""
            }
        )

        continue


    # --------------------------------------------------------
    # NO DIRECT OVERLAP
    # LOOK FOR NEAREST GENE
    # --------------------------------------------------------

    nearby = find_nearby_genes(
        hit["reference"],
        hit["start"],
        hit["end"]
    )


    if nearby:

        nearest = nearby[0]

        nearest_gene = nearest["gene"]

        results.append(
            {
                "contig_id": contig_id,
                "gene_id": "NO_DIRECT_OVERLAP",
                "assignment": "NEARBY_GENE",
                "reference": hit["reference"],
                "start": hit["start"],
                "end": hit["end"],
                "strand": hit["strand"],
                "identity": round(
                    hit["identity"],
                    2
                ),
                "coverage": round(
                    hit["coverage"],
                    2
                ),
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "gene_overlap": 0,
                "nearest_gene": nearest_gene[
                    "gene_id"
                ],
                "nearest_gene_distance": nearest[
                    "distance"
                ]
            }
        )

    else:

        results.append(
            {
                "contig_id": contig_id,
                "gene_id": "NO_GENE_FOUND",
                "assignment": "INTERGENIC_OR_UNANNOTATED",
                "reference": hit["reference"],
                "start": hit["start"],
                "end": hit["end"],
                "strand": hit["strand"],
                "identity": round(
                    hit["identity"],
                    2
                ),
                "coverage": round(
                    hit["coverage"],
                    2
                ),
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "gene_overlap": 0,
                "nearest_gene": "",
                "nearest_gene_distance": ""
            }
        )


# ============================================================
# WRITE OUTPUT
# ============================================================

print("\nWriting final CSV...")

fieldnames = [
    "contig_id",
    "gene_id",
    "assignment",
    "reference",
    "start",
    "end",
    "strand",
    "identity",
    "coverage",
    "evalue",
    "bitscore",
    "gene_overlap",
    "nearest_gene",
    "nearest_gene_distance"
]


with open(
    FINAL_OUTPUT,
    "w",
    newline=""
) as output:

    writer = csv.DictWriter(
        output,
        fieldnames=fieldnames
    )

    writer.writeheader()

    writer.writerows(results)


# ============================================================
# SUMMARY
# ============================================================

total = len(results)

direct = sum(
    1
    for r in results
    if r["assignment"]
    == "DIRECT_GENE_OVERLAP"
)

nearby = sum(
    1
    for r in results
    if r["assignment"]
    == "NEARBY_GENE"
)

no_blast = sum(
    1
    for r in results
    if r["assignment"]
    == "NO_BLAST_HIT"
)

intergenic = sum(
    1
    for r in results
    if r["assignment"]
    == "INTERGENIC_OR_UNANNOTATED"
)


print("\n============================================")
print("           ANALYSIS COMPLETE")
print("============================================")

print(
    f"Total contigs:              {total}"
)

print(
    f"Direct gene assignments:    {direct}"
)

print(
    f"Nearby gene candidates:     {nearby}"
)

print(
    f"No BLAST hit:               {no_blast}"
)

print(
    f"Intergenic/unannotated:     {intergenic}"
)

print("\nOutput:")
print(
    f"  {FINAL_OUTPUT}"
)

print("============================================")
