"""
================================================================
Day 17 — Genome Sequence Analysis (REAL DATA)
Author  : Subhadip Jana
Organism: Escherichia coli K-12 substr. MG1655
Accession: NC_000913.3 | 4,641,652 bp | Complete genome

Analyses:
  1. Basic genome statistics (length, GC%, N50 etc.)
  2. GC content sliding window (GC skew landscape)
  3. GC skew analysis [(G-C)/(G+C)] — identifies ori/ter
  4. Dinucleotide & trinucleotide frequency (k-mer analysis)
  5. Codon usage bias (all 64 codons, RSCU)
  6. Cumulative GC skew — pinpoints replication origin
  7. AT/GC isochore-like regions (compositional complexity)
  8. Sequence complexity (Shannon entropy sliding window)
================================================================
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import Counter
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────
# SECTION 1: PARSE FASTA
# ─────────────────────────────────────────────────────────────

print("🔬 Loading E. coli K-12 MG1655 genome...")

with open("data/ecoli_k12_genome.fasta") as f:
    lines  = f.readlines()

header = lines[0].strip().lstrip(">")
genome = "".join(l.strip().upper() for l in lines[1:])
L      = len(genome)

print(f"✅ {header[:60]}")
print(f"   Length : {L:,} bp")
print(f"   GC     : {(genome.count('G')+genome.count('C'))/L*100:.3f}%")
print(f"   AT     : {(genome.count('A')+genome.count('T'))/L*100:.3f}%")

# ─────────────────────────────────────────────────────────────
# SECTION 2: BASIC STATISTICS
# ─────────────────────────────────────────────────────────────

print("\n📊 Computing basic genome statistics...")

counts = {b: genome.count(b) for b in "ACGTN"}
gc     = (counts["G"] + counts["C"]) / L * 100
at     = (counts["A"] + counts["T"]) / L * 100
gc_skew_global = (counts["G"] - counts["C"]) / (counts["G"] + counts["C"])
at_skew_global = (counts["A"] - counts["T"]) / (counts["A"] + counts["T"])

stats = {
    "Organism"         : "E. coli K-12 MG1655",
    "Accession"        : "NC_000913.3",
    "Genome length (bp)": f"{L:,}",
    "GC content (%)"   : f"{gc:.3f}",
    "AT content (%)"   : f"{at:.3f}",
    "G count"          : f"{counts['G']:,}",
    "C count"          : f"{counts['C']:,}",
    "A count"          : f"{counts['A']:,}",
    "T count"          : f"{counts['T']:,}",
    "N count"          : f"{counts['N']:,}",
    "GC skew global"   : f"{gc_skew_global:.6f}",
    "AT skew global"   : f"{at_skew_global:.6f}",
    "Purine (A+G) %"   : f"{(counts['A']+counts['G'])/L*100:.3f}",
    "Pyrimidine (C+T) %": f"{(counts['C']+counts['T'])/L*100:.3f}",
}
for k, v in stats.items():
    print(f"   {k:25s}: {v}")

# ─────────────────────────────────────────────────────────────
# SECTION 3: SLIDING WINDOW GC CONTENT
# ─────────────────────────────────────────────────────────────

print("\n📈 Computing GC sliding window (10 kb windows)...")

WIN  = 10_000
STEP = 5_000
gc_windows, gc_pos = [], []

for start in range(0, L - WIN, STEP):
    w = genome[start:start+WIN]
    gc_w = (w.count("G") + w.count("C")) / len(w) * 100
    gc_windows.append(gc_w)
    gc_pos.append((start + WIN//2) / 1e6)   # Mb

gc_windows = np.array(gc_windows)
gc_pos     = np.array(gc_pos)
print(f"   Windows: {len(gc_windows)} | "
      f"GC range: {gc_windows.min():.1f}%–{gc_windows.max():.1f}%")

# ─────────────────────────────────────────────────────────────
# SECTION 4: GC SKEW + CUMULATIVE GC SKEW
# ─────────────────────────────────────────────────────────────

print("🔄 Computing GC skew and cumulative GC skew...")

SKEW_WIN  = 10_000
SKEW_STEP = 5_000
skew_vals, skew_pos = [], []

for start in range(0, L - SKEW_WIN, SKEW_STEP):
    w = genome[start:start+SKEW_WIN]
    g, c = w.count("G"), w.count("C")
    skew = (g - c) / (g + c) if (g + c) > 0 else 0
    skew_vals.append(skew)
    skew_pos.append((start + SKEW_WIN//2) / 1e6)

skew_vals  = np.array(skew_vals)
skew_pos   = np.array(skew_pos)
cum_skew   = np.cumsum(skew_vals)

# Ori = minimum of cumulative skew, Ter = maximum
ori_idx = np.argmin(cum_skew)
ter_idx = np.argmax(cum_skew)
ori_mb  = skew_pos[ori_idx]
ter_mb  = skew_pos[ter_idx]
print(f"   Predicted replication origin (oriC): ~{ori_mb:.2f} Mb")
print(f"   Predicted terminus (ter)           : ~{ter_mb:.2f} Mb")
print(f"   Known oriC position                : ~3.93 Mb (literature)")

# ─────────────────────────────────────────────────────────────
# SECTION 5: K-MER ANALYSIS
# ─────────────────────────────────────────────────────────────

print("\n🔤 Computing k-mer frequencies (di & trinucleotides)...")

# Dinucleotides
di_obs   = Counter(genome[i:i+2] for i in range(L-1)
                   if "N" not in genome[i:i+2])
di_total = sum(di_obs.values())

# Expected frequencies (product of monomer freqs)
mono_freq = {b: genome.count(b)/L for b in "ACGT"}
di_exp    = {a+b: mono_freq[a]*mono_freq[b]*di_total
             for a in "ACGT" for b in "ACGT"}

# Observed / Expected ratio (over-/under-representation)
di_oe = {k: di_obs.get(k,0)/di_exp[k] for k in di_exp if di_exp[k]>0}
di_df = pd.DataFrame({"Dinucleotide":list(di_oe.keys()),
                       "Observed":   [di_obs.get(k,0) for k in di_oe],
                       "Expected":   [di_exp[k]        for k in di_oe],
                       "OE_ratio":   list(di_oe.values())
                      }).sort_values("OE_ratio", ascending=False)

print("   Top over-represented dinucleotides:")
for _, row in di_df.head(5).iterrows():
    print(f"   {row['Dinucleotide']}: O/E={row['OE_ratio']:.4f}")
print("   Top under-represented dinucleotides:")
for _, row in di_df.tail(5).iterrows():
    print(f"   {row['Dinucleotide']}: O/E={row['OE_ratio']:.4f}")

# ─────────────────────────────────────────────────────────────
# SECTION 6: CODON USAGE (RSCU)
# ─────────────────────────────────────────────────────────────

print("\n🧬 Computing codon usage (RSCU)...")

# Count all codons (both strands, reading frame 0)
codon_counts = Counter()
for i in range(0, L-2, 3):
    codon = genome[i:i+3]
    if "N" not in codon and len(codon)==3:
        codon_counts[codon] += 1

# Genetic code (standard)
GENETIC_CODE = {
    "TTT":"Phe","TTC":"Phe","TTA":"Leu","TTG":"Leu",
    "CTT":"Leu","CTC":"Leu","CTA":"Leu","CTG":"Leu",
    "ATT":"Ile","ATC":"Ile","ATA":"Ile","ATG":"Met",
    "GTT":"Val","GTC":"Val","GTA":"Val","GTG":"Val",
    "TCT":"Ser","TCC":"Ser","TCA":"Ser","TCG":"Ser",
    "CCT":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "ACT":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "GCT":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "TAT":"Tyr","TAC":"Tyr","TAA":"Stop","TAG":"Stop",
    "CAT":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "AAT":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
    "GAT":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "TGT":"Cys","TGC":"Cys","TGA":"Stop","TGG":"Trp",
    "CGT":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "AGT":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
    "GGT":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly",
}

# RSCU = observed / (total for AA / n synonymous codons)
from collections import defaultdict
aa_codons = defaultdict(list)
for codon, aa in GENETIC_CODE.items():
    if aa != "Stop":
        aa_codons[aa].append(codon)

rscu = {}
for aa, codons in aa_codons.items():
    total = sum(codon_counts.get(c,0) for c in codons)
    n     = len(codons)
    for c in codons:
        obs  = codon_counts.get(c, 0)
        rscu[c] = (obs / (total/n)) if total > 0 else 0

rscu_df = pd.DataFrame({"Codon":list(rscu.keys()),
                         "AA":   [GENETIC_CODE[c] for c in rscu],
                         "Count":[codon_counts.get(c,0) for c in rscu],
                         "RSCU": list(rscu.values())
                        }).sort_values("RSCU", ascending=False)

print("   Top preferred codons (RSCU > 1.5):")
for _, row in rscu_df[rscu_df["RSCU"]>1.5].head(8).iterrows():
    print(f"   {row['Codon']} ({row['AA']:3s}): RSCU={row['RSCU']:.3f}  n={row['Count']:,}")

# ─────────────────────────────────────────────────────────────
# SECTION 7: SHANNON ENTROPY COMPLEXITY
# ─────────────────────────────────────────────────────────────

print("\n🔢 Computing sequence complexity (Shannon entropy)...")

ENT_WIN  = 5_000
ENT_STEP = 2_500
ent_vals, ent_pos = [], []

for start in range(0, L - ENT_WIN, ENT_STEP):
    w    = genome[start:start+ENT_WIN]
    freq = np.array([w.count(b)/ENT_WIN for b in "ACGT"])
    freq = freq[freq > 0]
    H    = -np.sum(freq * np.log2(freq))
    ent_vals.append(H)
    ent_pos.append((start + ENT_WIN//2) / 1e6)

ent_vals = np.array(ent_vals)
ent_pos  = np.array(ent_pos)
print(f"   Mean entropy: {ent_vals.mean():.4f} bits (max=2.0 for uniform)")
print(f"   Min entropy : {ent_vals.min():.4f} (low-complexity region)")

# ─────────────────────────────────────────────────────────────
# SECTION 8: SAVE DATA FILES
# ─────────────────────────────────────────────────────────────

pd.DataFrame({"position_Mb":gc_pos,"GC_pct":gc_windows}).to_csv(
    "outputs/gc_content_windows.csv", index=False)
pd.DataFrame({"position_Mb":skew_pos,"GC_skew":skew_vals,
               "cumulative_GC_skew":cum_skew}).to_csv(
    "outputs/gc_skew.csv", index=False)
di_df.to_csv("outputs/dinucleotide_oe.csv",  index=False)
rscu_df.sort_values("Codon").to_csv("outputs/codon_usage_rscu.csv", index=False)
pd.DataFrame({"position_Mb":ent_pos,"shannon_entropy":ent_vals}).to_csv(
    "outputs/sequence_complexity.csv", index=False)
pd.DataFrame(list(stats.items()), columns=["Metric","Value"]).to_csv(
    "outputs/genome_stats.csv", index=False)
print("\n✅ All data files saved")

# ─────────────────────────────────────────────────────────────
# SECTION 9: DASHBOARD (9 panels)
# ─────────────────────────────────────────────────────────────

print("\n🎨 Generating dashboard...")

fig = plt.figure(figsize=(24, 20))
fig.suptitle(
    "Genome Sequence Analysis — Escherichia coli K-12 MG1655\n"
    "NC_000913.3 | 4,641,652 bp | Complete Genome | REAL DATA\n"
    "GC content · GC skew · k-mer · Codon usage · Sequence complexity",
    fontsize=15, fontweight="bold", y=0.99
)

# ── Plot 1: GC content sliding window ──
ax1 = fig.add_subplot(3, 3, 1)
ax1.fill_between(gc_pos, gc_windows, gc, alpha=0.3,
                 color=np.where(gc_windows >= gc, "#E74C3C", "#3498DB")[0])
ax1.fill_between(gc_pos, gc_windows, gc,
                 where=gc_windows >= gc, alpha=0.4, color="#E74C3C", label="GC > mean")
ax1.fill_between(gc_pos, gc_windows, gc,
                 where=gc_windows < gc,  alpha=0.4, color="#3498DB", label="GC < mean")
ax1.plot(gc_pos, gc_windows, lw=0.8, color="#2C3E50", alpha=0.6)
ax1.axhline(gc, color="black", lw=1.5, linestyle="--",
            label=f"Mean GC={gc:.2f}%")
ax1.set_xlabel("Genome position (Mb)")
ax1.set_ylabel("GC content (%)")
ax1.set_title("GC Content Sliding Window\n(10 kb windows, 5 kb step)",
              fontweight="bold", fontsize=10)
ax1.legend(fontsize=8)
ax1.set_xlim(0, L/1e6)

# ── Plot 2: GC skew ──
ax2 = fig.add_subplot(3, 3, 2)
ax2.fill_between(skew_pos, skew_vals, 0,
                 where=skew_vals >= 0, alpha=0.5, color="#E74C3C", label="G > C (leading)")
ax2.fill_between(skew_pos, skew_vals, 0,
                 where=skew_vals < 0,  alpha=0.5, color="#3498DB", label="C > G (lagging)")
ax2.axhline(0, color="black", lw=1)
ax2.set_xlabel("Genome position (Mb)")
ax2.set_ylabel("GC skew [(G−C)/(G+C)]")
ax2.set_title("GC Skew Analysis\n(identifies replication orientation)",
              fontweight="bold", fontsize=10)
ax2.legend(fontsize=8)
ax2.set_xlim(0, L/1e6)

# ── Plot 3: Cumulative GC skew (ori/ter detection) ──
ax3 = fig.add_subplot(3, 3, 3)
ax3.plot(skew_pos, cum_skew, lw=2, color="#9B59B6")
ax3.fill_between(skew_pos, cum_skew, alpha=0.15, color="#9B59B6")
ax3.axvline(ori_mb, color="#E74C3C", lw=2, linestyle="--",
            label=f"Predicted oriC: {ori_mb:.2f} Mb")
ax3.axvline(ter_mb, color="#3498DB", lw=2, linestyle="--",
            label=f"Predicted ter: {ter_mb:.2f} Mb")
ax3.axvline(3.93, color="#2ECC71", lw=1.5, linestyle=":",
            label="Known oriC: 3.93 Mb")
ax3.set_xlabel("Genome position (Mb)")
ax3.set_ylabel("Cumulative GC Skew")
ax3.set_title("Cumulative GC Skew\n(min=oriC, max=terminus)",
              fontweight="bold", fontsize=10)
ax3.legend(fontsize=8)
ax3.set_xlim(0, L/1e6)

# ── Plot 4: Nucleotide composition bar ──
ax4 = fig.add_subplot(3, 3, 4)
bases  = ["A","T","G","C"]
bcols  = ["#3498DB","#E74C3C","#2ECC71","#F39C12"]
bvals  = [counts[b]/L*100 for b in bases]
bars4  = ax4.bar(bases, bvals, color=bcols, edgecolor="black",
                 linewidth=0.6, alpha=0.87, width=0.5)
for bar, val in zip(bars4, bvals):
    ax4.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.1,
             f"{val:.2f}%", ha="center", fontsize=11, fontweight="bold")
ax4.set_ylabel("Frequency (%)")
ax4.set_title(f"Nucleotide Composition\nGC={gc:.2f}% | AT={at:.2f}%",
              fontweight="bold", fontsize=10)
ax4.set_ylim(0, 30)
ax4.axhline(25, color="gray", lw=1, linestyle="--", alpha=0.5,
            label="Expected 25% (uniform)")
ax4.legend(fontsize=8)

# ── Plot 5: Dinucleotide O/E heatmap ──
ax5 = fig.add_subplot(3, 3, 5)
oe_matrix = np.zeros((4,4))
bases4    = list("ACGT")
for i, b1 in enumerate(bases4):
    for j, b2 in enumerate(bases4):
        oe_matrix[i,j] = di_oe.get(b1+b2, 1.0)
sns.heatmap(oe_matrix, ax=ax5, cmap="RdBu_r", center=1.0,
            vmin=0.5, vmax=1.5, annot=True, fmt=".2f",
            xticklabels=bases4, yticklabels=bases4,
            linewidths=0.5, cbar_kws={"label":"O/E ratio","shrink":0.8},
            annot_kws={"size":10})
ax5.set_title("Dinucleotide O/E Ratios\n(>1=over-represented, <1=under-represented)",
              fontweight="bold", fontsize=10)
ax5.set_xlabel("3′ base"); ax5.set_ylabel("5′ base")

# ── Plot 6: RSCU heatmap (all sense codons) ──
ax6 = fig.add_subplot(3, 3, 6)
sense_rscu = rscu_df[rscu_df["AA"] != "Stop"].copy()
# Pivot: AA as rows, codon position as columns
aa_order = sorted(sense_rscu["AA"].unique())
codon_pivot = {}
for _, row in sense_rscu.iterrows():
    codon_pivot[(row["AA"], row["Codon"])] = row["RSCU"]
# Build matrix sorted by AA
aa_groups = sense_rscu.groupby("AA")["Codon"].apply(list).to_dict()
rows_rscu, row_labels = [], []
for aa in aa_order:
    codons = sorted(aa_groups[aa])
    for c in codons:
        rows_rscu.append(rscu_df[rscu_df["Codon"]==c]["RSCU"].values[0])
        row_labels.append(f"{c}\n({aa})")
rscu_arr = np.array(rows_rscu).reshape(-1, 1)
im = ax6.imshow(rscu_arr, cmap="RdYlGn", aspect="auto",
                vmin=0, vmax=3)
ax6.set_yticks(range(len(row_labels)))
ax6.set_yticklabels(row_labels, fontsize=5)
ax6.set_xticks([]); ax6.set_xlabel("")
plt.colorbar(im, ax=ax6, label="RSCU", shrink=0.8)
ax6.set_title("Codon Usage Bias (RSCU)\n(green=preferred, red=avoided)",
              fontweight="bold", fontsize=10)

# ── Plot 7: Shannon entropy ──
ax7 = fig.add_subplot(3, 3, 7)
ax7.plot(ent_pos, ent_vals, lw=1, color="#1ABC9C", alpha=0.7)
ax7.fill_between(ent_pos, ent_vals, ent_vals.min(),
                 alpha=0.2, color="#1ABC9C")
ax7.axhline(ent_vals.mean(), color="#E74C3C", lw=1.5,
            linestyle="--", label=f"Mean H={ent_vals.mean():.3f}")
ax7.axhline(2.0, color="gray", lw=1, linestyle=":",
            label="Max H=2.0 (uniform)")
# Highlight low-complexity regions
low_thresh = ent_vals.mean() - 2*ent_vals.std()
ax7.fill_between(ent_pos, ent_vals, low_thresh,
                 where=ent_vals<low_thresh, alpha=0.5,
                 color="#E74C3C", label="Low complexity")
ax7.set_xlabel("Genome position (Mb)")
ax7.set_ylabel("Shannon Entropy H (bits)")
ax7.set_title("Sequence Complexity (Shannon Entropy)\n(5 kb windows)",
              fontweight="bold", fontsize=10)
ax7.legend(fontsize=8)
ax7.set_xlim(0, L/1e6)

# ── Plot 8: Top 20 preferred / avoided codons ──
ax8 = fig.add_subplot(3, 3, 8)
top10  = rscu_df[rscu_df["AA"]!="Stop"].nlargest(10,"RSCU")
bot10  = rscu_df[rscu_df["AA"]!="Stop"].nsmallest(10,"RSCU")
combo  = pd.concat([top10, bot10]).reset_index(drop=True)
colors8= ["#2ECC71" if r>1 else "#E74C3C" for r in combo["RSCU"]]
labels8= [f"{row['Codon']}\n({row['AA']})" for _,row in combo.iterrows()]
ax8.barh(range(len(combo)), combo["RSCU"].values,
         color=colors8, edgecolor="black", linewidth=0.3, alpha=0.87)
ax8.set_yticks(range(len(combo)))
ax8.set_yticklabels(labels8, fontsize=7)
ax8.axvline(1.0, color="black", lw=1.5)
ax8.set_xlabel("RSCU")
ax8.set_title("Preferred (RSCU>1) vs Avoided (RSCU<1) Codons\n(Top/Bottom 10)",
              fontweight="bold", fontsize=10)

# ── Plot 9: Summary stats table ──
ax9 = fig.add_subplot(3, 3, 9)
ax9.axis("off")
rows9 = list(stats.items())
tbl9  = ax9.table(cellText=rows9,
                   colLabels=["Metric","Value"],
                   cellLoc="left", loc="center")
tbl9.auto_set_font_size(False); tbl9.set_fontsize(8.5); tbl9.scale(1.8, 1.85)
for j in range(2):
    tbl9[(0,j)].set_facecolor("#2C3E50")
    tbl9[(0,j)].set_text_props(color="white", fontweight="bold")
for i in range(1, len(rows9)+1, 2):
    for j in range(2):
        tbl9[(i,j)].set_facecolor("#F2F3F4")
ax9.set_title("E. coli K-12 MG1655 — Genome Stats",
              fontweight="bold", fontsize=11, pad=20)

plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig("outputs/genome_analysis_dashboard.png",
            dpi=150, bbox_inches="tight")
plt.close()
print("✅ Dashboard saved → outputs/genome_analysis_dashboard.png")

# ─────────────────────────────────────────────────────────────
# FINAL SUMMARY
# ─────────────────────────────────────────────────────────────

print("\n" + "="*60)
print("FINAL SUMMARY — E. coli K-12 GENOME ANALYSIS")
print("="*60)
print(f"\nGenome length  : {L:,} bp ({L/1e6:.3f} Mb)")
print(f"GC content     : {gc:.3f}%")
print(f"GC skew (global): {gc_skew_global:+.6f}")
print(f"Predicted oriC : {ori_mb:.2f} Mb  (known: 3.93 Mb)")
print(f"Predicted ter  : {ter_mb:.2f} Mb")
print(f"\nK-mer analysis:")
print(f"  Most over-represented  : {di_df.iloc[0]['Dinucleotide']} "
      f"(O/E={di_df.iloc[0]['OE_ratio']:.4f})")
print(f"  Most under-represented : {di_df.iloc[-1]['Dinucleotide']} "
      f"(O/E={di_df.iloc[-1]['OE_ratio']:.4f})")
print(f"\nCodon usage (top preferred codons):")
for _, row in rscu_df[rscu_df["AA"]!="Stop"].nlargest(5,"RSCU").iterrows():
    print(f"  {row['Codon']} ({row['AA']:3s}): RSCU={row['RSCU']:.3f}")
print(f"\nSequence complexity:")
print(f"  Mean entropy : {ent_vals.mean():.4f} bits")
print(f"  Low-complexity windows: "
      f"{(ent_vals < ent_vals.mean()-2*ent_vals.std()).sum()}")
print("\n✅ All outputs saved!")
