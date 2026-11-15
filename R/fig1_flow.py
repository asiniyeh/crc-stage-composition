"""Figure 1: study design and sample flow (numbers from results/flow_counts.csv)."""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

out = sys.argv[1]
fig, ax = plt.subplots(figsize=(11, 7.2))
ax.set_xlim(0, 110); ax.set_ylim(0, 72); ax.axis("off")

def box(x, y, w, h, text, bold_first=True, fc="white", ec="#444444", fs=9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=ec, lw=1))
    lines = text.split("\n")
    n = len(lines)
    for i, line in enumerate(lines):
        ax.text(x + w / 2, y + h - (i + 0.85) * h / (n + 0.3), line, ha="center", va="center", fontsize=fs,
                fontweight="bold" if (i == 0 and bold_first) else "normal")

def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", color="#444444", lw=1))

def excl(x, y, text):
    ax.text(x, y, text, fontsize=7.8, color="#555555", va="center", style="italic")

# column headers
for x, t in [(2, "Discovery: TCGA-COAD/READ (RNA-seq)"), (40, "Validation: GSE39582 (array)"), (76, "Validation: GSE17536 (array)")]:
    ax.text(x, 69.5, t, fontsize=10, fontweight="bold", color="#0B4F8A")

# TCGA column
box(2, 58, 30, 7, "701 RNA-seq samples\nopen-access STAR counts")
arrow(17, 58, 17, 53.6); excl(19, 55.8, "− 54 normal / recurrent / metastatic")
box(2, 46, 30, 7, "647 primary tumours")
arrow(17, 46, 17, 41.6); excl(19, 43.8, "− 23 repeat aliquots")
box(2, 34, 30, 7, "624 patients")
arrow(17, 34, 17, 29.6); excl(19, 31.8, "− 23 stage missing")
box(2, 22, 30, 7, "601 tumours with stage I–IV\nTable 1, Figure 2", fc="#EAF2FA", ec="#0B4F8A")
arrow(17, 22, 17, 17.6); excl(19, 19.8, "− 104 side / CMS / score missing")
box(2, 6, 30, 11, "497 tumours (267 early, 230 late)\nDESeq2 model A vs model B\npanel, GSEA, classifier training\nsensitivity analyses S1–S5", fc="#EAF2FA", ec="#0B4F8A")

# GSE39582
box(40, 58, 30, 7, "585 samples")
arrow(55, 58, 55, 47.6); excl(57, 52.8, "− 19 non-tumoral mucosa")
box(40, 40, 30, 7, "566 colon tumours")
arrow(55, 40, 55, 29.6); excl(57, 34.8, "− 4 stage 0")
box(40, 22, 30, 7, "562 tumours with stage I–IV\nTable 1, Figure 2", fc="#FDF1E6", ec="#B35C00")

# GSE17536
box(76, 58, 30, 7, "177 samples")
arrow(91, 58, 91, 29.6)
box(76, 22, 30, 7, "177 tumours with stage I–IV\nTable 1, Figure 2", fc="#FDF1E6", ec="#B35C00")

# validation box spanning both array cohorts
box(40, 6, 66, 11, "External validation (limma, same adjustment)\npanel replication and meta-analysis (Figure 5)\nfrozen classifier (Figure S1) · panel score vs RFS / OS (Table 3, Figure 6)",
    fc="#FDF1E6", ec="#B35C00")
arrow(55, 22, 55, 17.6); arrow(91, 22, 91, 17.6)
arrow(32, 11.5, 39.3, 11.5)
ax.text(35.6, 13.2, "frozen\npanel", fontsize=7.5, ha="center", color="#555555")

plt.savefig(out + ".png", dpi=300, bbox_inches="tight")
plt.savefig(out + ".pdf", bbox_inches="tight")
print("saved", out)
