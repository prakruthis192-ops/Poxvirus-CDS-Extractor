import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

CSV_IN = "results/orf_validation_summary.csv"
FIG_OUT = "results/figures/mpxv_cds_conservation.png"


def plot_conservation():
    df = pd.read_csv(CSV_IN)
    os.makedirs("results/figures", exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.3, 1]})

    df["Status_Code"] = df["Intact_ORF"].astype(int)
    pivot_df = df.pivot(index="Genome", columns="Gene", values="Status_Code")
    gene_order = ["OPG057", "OPG210", "OPG085", "OPG032"]
    pivot_df = pivot_df[gene_order]

    cmap = sns.color_palette(["#d95f02", "#1b9e77"])
    sns.heatmap(
        pivot_df,
        cmap=cmap,
        cbar=False,
        linewidths=0.8,
        linecolor="white",
        ax=ax1
    )
    ax1.set_title("MPXV Target CDS Presence & Reading-Frame Integrity", fontweight="bold", pad=12)
    ax1.set_xlabel("Target Ortholog (OPG Locus)", fontsize=10)
    ax1.set_ylabel("MPXV Genome Accession", fontsize=10)
    ax1.set_xticklabels(
        ["OPG057\n(F13L)", "OPG210\n(B21R)", "OPG085\n(O2L)", "OPG032\n(C3L)"],
        fontsize=9
    )
    ax1.tick_params(axis="y", labelsize=8)

    legend_elements = [
        Patch(facecolor="#1b9e77", label="Intact CDS (100% ORF)"),
        Patch(facecolor="#d95f02", label="Absent / Clade II Deletion")
    ]
    ax1.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)

    intact_df = df[df["Intact_ORF"] == True].copy()
    sns.boxplot(
        data=intact_df,
        x="Gene",
        y="SNP_Count_vs_Ref",
        order=["OPG057", "OPG210", "OPG085"],
        hue="Gene",
        palette=["#377eb8", "#984ea3", "#4daf4a"],
        legend=False,
        width=0.45,
        ax=ax2
    )
    sns.stripplot(
        data=intact_df,
        x="Gene",
        y="SNP_Count_vs_Ref",
        order=["OPG057", "OPG210", "OPG085"],
        color="black",
        size=6,
        alpha=0.7,
        ax=ax2
    )
    ax2.set_title("Nucleotide Divergence vs. Reference (NC_063383.1)", fontweight="bold", pad=12)
    ax2.set_xlabel("Retained MPXV Gene", fontsize=10)
    ax2.set_ylabel("SNP Count Relative to Reference", fontsize=10)
    ax2.set_xticklabels(["OPG057\n(F13L)", "OPG210\n(B21R)", "OPG085\n(O2L)"], fontsize=9)
    ax2.grid(axis="y", linestyle="--", alpha=0.4)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(FIG_OUT, dpi=300, bbox_inches="tight")
    print(f"[Plot] Saved CDS conservation figure to {FIG_OUT}")


if __name__ == "__main__":
    plot_conservation()
