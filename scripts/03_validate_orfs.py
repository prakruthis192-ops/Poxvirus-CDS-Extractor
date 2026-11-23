import os
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

REF_GENES = "data/reference_genes/mpxv_target_genes.fasta"
EXTRACTED_DIR = "data/extracted_cds"
GENOMES_FASTA = "data/raw/mpxv_complete_genomes.fasta"
CSV_OUT = "results/orf_validation_summary.csv"

GENE_FUNCTIONS = {
    "OPG057": "F13L (Antiviral Target - Tecovirimat)",
    "OPG210": "B21R (Surface Glycoprotein / Immune Evasion)",
    "OPG085": "O2L (Core Structural / Housekeeping)",
    "OPG032": "C3L (Complement Control / Virulence)"
}


def clean_accession(raw_id):
    return raw_id.replace("[", "").replace("]", "").replace("'", "").replace(",", "").split(".")[0].strip()


def validate_extracted_orfs():
    ref_seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(REF_GENES, "fasta")}
    genome_ids = [clean_accession(rec.id) for rec in SeqIO.parse(GENOMES_FASTA, "fasta")]

    extracted_map = {}
    for fasta_path in glob.glob(os.path.join(EXTRACTED_DIR, "*_orthologs.fasta")):
        gene = os.path.basename(fasta_path).replace("_orthologs.fasta", "")
        for rec in SeqIO.parse(fasta_path, "fasta"):
            g_acc = clean_accession(rec.id)
            extracted_map[(g_acc, gene)] = str(rec.seq).upper()

    rows = []
    for gene in GENE_FUNCTIONS:
        ref_str = ref_seqs.get(gene, "")
        for g_acc in genome_ids:
            seq_str = extracted_map.get((g_acc, gene), "")
            if not seq_str:
                rows.append({
                    "Genome": g_acc,
                    "Gene": gene,
                    "Functional_Role": GENE_FUNCTIONS[gene],
                    "Status": "Absent / Deleted",
                    "Length_bp": 0,
                    "Intact_ORF": False,
                    "SNP_Count_vs_Ref": None,
                    "Pct_Identity": 0.0
                })
                continue

            starts_atg = seq_str.startswith("ATG")
            div_by_3 = (len(seq_str) % 3 == 0)
            body = seq_str[:-3] if seq_str[-3:] in {"TAA", "TAG", "TGA"} else seq_str
            no_internal_stop = ("*" not in str(Seq(body).translate())) if div_by_3 else False
            intact = starts_atg and div_by_3 and no_internal_stop

            if ref_str and len(seq_str) == len(ref_str):
                snps = sum(1 for a, b in zip(seq_str, ref_str) if a != b)
                pident = round(100.0 * (1 - snps / len(ref_str)), 3)
            else:
                snps = None
                pident = None

            rows.append({
                "Genome": g_acc,
                "Gene": gene,
                "Functional_Role": GENE_FUNCTIONS[gene],
                "Status": "Intact CDS" if intact else "Truncated / Frameshift",
                "Length_bp": len(seq_str),
                "Intact_ORF": intact,
                "SNP_Count_vs_Ref": snps,
                "Pct_Identity": pident
            })

    df = pd.DataFrame(rows)
    os.makedirs("results", exist_ok=True)
    df.to_csv(CSV_OUT, index=False)
    print(f"[Validation] Evaluated {len(df)} gene-genome pairs -> {CSV_OUT}")
    print(df.groupby(["Gene", "Status"]).size())


if __name__ == "__main__":
    validate_extracted_orfs()
