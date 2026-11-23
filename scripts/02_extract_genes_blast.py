import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

REF_GENES = "data/reference_genes/mpxv_target_genes.fasta"
GENOMES_FASTA = "data/raw/mpxv_complete_genomes.fasta"
DB_PREFIX = "data/blast_db/mpxv_cohort"
BLAST_TSV = "results/blast_hits.tsv"
OUT_DIR = "data/extracted_cds"


def run_local_blast():
    os.makedirs("data/blast_db", exist_ok=True)
    os.makedirs("results", exist_ok=True)

    print("[BLAST+] Building local nucleotide database from MPXV genomes...")
    subprocess.run(
        ["makeblastdb", "-in", GENOMES_FASTA, "-dbtype", "nucl", "-out", DB_PREFIX],
        stdout=subprocess.DEVNULL,
        check=True
    )

    print("[BLAST+] Searching target genes (OPG057, OPG210, OPG032, OPG085) across cohort...")
    cmd = [
        "blastn",
        "-query", REF_GENES,
        "-db", DB_PREFIX,
        "-out", BLAST_TSV,
        "-outfmt", "6 qseqid sseqid pident length qlen sstart send evalue",
        "-evalue", "1e-20",
        "-max_target_seqs", "50"
    ]
    subprocess.run(cmd, check=True)
    print(f"[BLAST+] Alignment coordinates saved to {BLAST_TSV}")


def extract_orthologs():
    os.makedirs(OUT_DIR, exist_ok=True)
    genomes = {rec.id: rec.seq for rec in SeqIO.parse(GENOMES_FASTA, "fasta")}

    cols = ["gene", "genome_id", "pident", "length", "qlen", "sstart", "send", "evalue"]
    df = pd.read_csv(BLAST_TSV, sep="\t", names=cols)

    # Keep the top hit per (gene, genome_id) pair
    df = df.sort_values(by=["gene", "genome_id", "length", "pident"], ascending=[True, True, False, False])
    best_hits = df.drop_duplicates(subset=["gene", "genome_id"])

    for gene, group in best_hits.groupby("gene"):
        extracted_records = []
        for _, row in group.iterrows():
            g_id = row["genome_id"]
            full_seq = genomes[g_id]
            start, end = int(row["sstart"]), int(row["send"])

            if start <= end:
                subseq = full_seq[start - 1 : end]
                strand = "+"
            else:
                subseq = full_seq[end - 1 : start].reverse_complement()
                strand = "-"

            clean_acc = g_id.split(".")
            rec_id = f"{clean_acc}_{gene}"
            desc = f"strand={strand} coords={start}..{end} pident={row['pident']}%"
            extracted_records.append(SeqRecord(subseq, id=rec_id, description=desc))

        out_fasta = os.path.join(OUT_DIR, f"{gene}_orthologs.fasta")
        SeqIO.write(extracted_records, out_fasta, "fasta")
        print(f"  -> {gene}: Extracted {len(extracted_records)} orthologs to {out_fasta}")


if __name__ == "__main__":
    run_local_blast()
    extract_orthologs()
