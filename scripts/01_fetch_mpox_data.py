import os
import time
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "prakruthis192@gmail.com"

REF_ACCESSION = "NC_063383.1"
TARGET_GENES = ["OPG057", "OPG210", "OPG032", "OPG085"]

GENOME_QUERY = (
    '(Monkeypox virus[Organism] OR Mpox virus[Organism]) AND '
    '"complete genome"[Title] AND 190000:205000[SLEN]'
)
MAX_GENOMES = 15

REF_GENES_OUT = "data/reference_genes/mpxv_target_genes.fasta"
GENOMES_OUT = "data/raw/mpxv_complete_genomes.fasta"


def fetch_reference_genes():
    """Download annotated MPXV reference GenBank record and extract target CDS."""
    print(f"[NCBI] Fetching reference annotation {REF_ACCESSION}...")
    handle = Entrez.efetch(db="nuccore", id=REF_ACCESSION, rettype="gbwithparts", retmode="text")
    ref_record = SeqIO.read(handle, "genbank")
    handle.close()

    extracted = []
    for feature in ref_record.features:
        if feature.type == "CDS":
            gene_names = feature.qualifiers.get("gene", []) + feature.qualifiers.get("locus_tag", [])
            for target in TARGET_GENES:
                if any(target in g for g in gene_names):
                    cds_seq = feature.extract(ref_record.seq)
                    product = feature.qualifiers.get("product", ["unknown_product"])[0]
                    rec = SeqRecord(cds_seq, id=target, description=product)
                    extracted.append(rec)
                    print(f"  -> Extracted reference {target} ({len(cds_seq)} bp): {product}")
                    break

    os.makedirs(os.path.dirname(REF_GENES_OUT), exist_ok=True)
    SeqIO.write(extracted, REF_GENES_OUT, "fasta")
    print(f"[Saved] {len(extracted)} reference genes saved to {REF_GENES_OUT}")


def fetch_cohort_genomes():
    """Download complete MPXV genomes from NCBI across outbreak cohorts."""
    print(f"[NCBI] Searching for {MAX_GENOMES} complete MPXV genomes...")
    handle = Entrez.esearch(db="nuccore", term=GENOME_QUERY, retmax=MAX_GENOMES, sort="relevance")
    search_res = Entrez.read(handle)
    handle.close()

    id_list = search_res["IdList"]
    print(f"[NCBI] Found {len(id_list)} genomes. Downloading FASTA (this takes ~10 seconds)...")
    time.sleep(0.5)

    fetch_handle = Entrez.efetch(db="nuccore", id=",".join(id_list), rettype="fasta", retmode="text")
    raw_fasta = fetch_handle.read()
    fetch_handle.close()

    os.makedirs(os.path.dirname(GENOMES_OUT), exist_ok=True)
    with open(GENOMES_OUT, "w") as f:
        f.write(raw_fasta)

    genomes = list(SeqIO.parse(GENOMES_OUT, "fasta"))
    print(f"[Saved] {len(genomes)} complete MPXV genomes saved to {GENOMES_OUT}")


if __name__ == "__main__":
    fetch_reference_genes()
    fetch_cohort_genomes()
