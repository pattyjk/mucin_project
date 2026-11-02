#!/usr/bin/env python
"""Dot-plot alpha diversity, Bray-Curtis PCoA and beta-dispersion figures.
Usage: python figures_alpha_pcoa_dispersion.py results/   (directory written by analysis_20C_mucin.R)"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = sys.argv[1] if len(sys.argv) > 1 else "results"
rng = np.random.default_rng(1)
plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 8, "axes.labelsize": 8})
cm = {"No": "#B0B0B0", "Yes": "#1f6fb4"}
def letter(ax, s): ax.text(-0.12, 1.08, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
abo = ["NoAB","Ampicillin","Chloramphenicol","Gentamycin","Neomycin","PSN","Streptomycin","Sulfamethazine","Trimethoprim"]
cso = ["Fucose","Galactose","Glucose","Glycerol","Mannitol","Methanol","N-acetyl glucosamine","N-acetylneuraminic acid","Tryptone","VitaminMix","YE","Mucin"]
sh = {"NoAB":"No AB","Ampicillin":"Amp","Chloramphenicol":"Chlor","Gentamycin":"Gent","Neomycin":"Neo","PSN":"PSN","Streptomycin":"Strep","Sulfamethazine":"Sulfa",
      "Trimethoprim":"Trim","N-acetyl glucosamine":"GlcNAc","N-acetylneuraminic acid":"Neu5Ac","VitaminMix":"Vit. mix","YE":"Yeast ext.","Mucin":"Mucin only"}
EXN = {"Antibiotics": "Antibiotics\n(Synthetic Amphibian Skin Media)", "CarbonTemp": "Single carbon\nsources"}

# ---- Figure 1: alpha + PCoA
al = pd.read_csv(os.path.join(R, "alpha_rarefied1000.csv"))
fig, axs = plt.subplots(2, 3, figsize=(7.4, 5.2), gridspec_kw=dict(width_ratios=[1.6,1.6,1.1], hspace=0.55, wspace=0.35))
for col, (ex, order) in enumerate([("Antibiotics", abo), ("CarbonTemp", cso)]):
    for row, v in enumerate(["Observed", "Shannon"]):
        ax = axs[row, col]
        for i, t in enumerate(order):
            for mu in ["No", "Yes"]:
                y = al[(al.Experiment==ex)&(al.Treatment==t)&(al.Mucin==mu)][v].values
                if len(y) == 0: continue
                x = i + (-0.18 if mu=="No" else 0.18)
                ax.scatter(x + rng.uniform(-0.06,0.06,len(y)), y, s=9, color=cm[mu], edgecolor="k", lw=0.3, zorder=3,
                           label=("Mucin" if mu=="Yes" else "No mucin") if (i==0 and row==0 and col==0) else None)
                ax.hlines(np.median(y), x-0.13, x+0.13, color="k", lw=1, zorder=4)
        ax.set_xticks(range(len(order))); ax.set_xticklabels([sh.get(o,o) for o in order] if row==1 else [], rotation=55, ha="right", fontsize=6)
        if col == 0: ax.set_ylabel("Observed ASVs (1,000 reads)" if v=="Observed" else "Shannon index")
        if row == 0: ax.set_title(EXN[ex], loc="left")
        ax.margins(x=0.03)
axs[0,0].legend(frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(0, 1.32))
for row, ex in enumerate(["Antibiotics", "CarbonTemp"]):
    ax = axs[row, 2]; p = pd.read_csv(os.path.join(R, f"pcoa_{ex}.csv"))
    for mu in ["No", "Yes"]:
        q = p[p.Mucin==mu]; ax.scatter(q.PC1, q.PC2, s=10, color=cm[mu], edgecolor="k", lw=0.3)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(f"PCoA 1 ({p.ve1[0]:.0f}%)"); ax.set_ylabel(f"PCoA 2 ({p.ve2[0]:.0f}%)")
    ax.set_title(("Antibiotics" if ex=="Antibiotics" else "Carbon sources") + ", Bray–Curtis", loc="left", fontsize=7)
for ax, l in zip([axs[0,0],axs[0,1],axs[0,2],axs[1,0],axs[1,1],axs[1,2]], "acebdf"): letter(ax, l)
fig.savefig(os.path.join(R, "mucin_20C_summary.png"), dpi=300, bbox_inches="tight")

# ---- Figure 2: beta dispersion
bd = pd.read_csv(os.path.join(R, "betadisper_distances.csv")); bd = bd[bd.distance=="bray"]
tests = pd.read_csv(os.path.join(R, "betadisper_tests.csv")).query("distance=='bray'").set_index("Experiment")
fig, axs = plt.subplots(1, 2, figsize=(5.6, 2.8), sharey=True, gridspec_kw=dict(wspace=0.08))
titles = {"dist_mucin_centroid": "Spread around mucin-level centroid\n(all treatments pooled)",
          "dist_group_centroid": "Spread among replicates\n(centroid per treatment × mucin)"}
pcol = {"dist_mucin_centroid": "p_mucinCentroid", "dist_group_centroid": "p_groupCentroid_mucin"}
for ax, col in zip(axs, titles):
    for i, ex in enumerate(["Antibiotics", "CarbonTemp"]):
        for mu, off in [("No",-0.18), ("Yes",0.18)]:
            y = bd[(bd.Experiment==ex)&(bd.Mucin==mu)][col].values
            ax.scatter(i+off+rng.uniform(-0.07,0.07,len(y)), y, s=8, color=cm[mu], edgecolor="k", lw=0.3, zorder=3,
                       label=("Mucin" if mu=="Yes" else "No mucin") if (i==0 and col=="dist_mucin_centroid") else None)
            ax.hlines(np.mean(y), i+off-0.14, i+off+0.14, color="k", lw=1.2, zorder=4)
        ax.text(i, 1.0, f"p = {tests.loc[ex, pcol[col]]:.2g}", ha="center", fontsize=6)
    ax.set_xticks([0,1]); ax.set_xticklabels(["Antibiotics","Carbon sources"]); ax.set_title(titles[col], loc="left", fontsize=7); ax.set_ylim(0, 1.06)
axs[0].set_ylabel("Bray–Curtis distance to centroid"); axs[0].legend(frameon=False, loc="lower left", fontsize=6)
letter(axs[0], "a"); letter(axs[1], "b")
fig.savefig(os.path.join(R, "betadisper_mucin.png"), dpi=300, bbox_inches="tight")
print("figures written to", R)
