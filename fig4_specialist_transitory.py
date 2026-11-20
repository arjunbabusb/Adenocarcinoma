#!/usr/bin/env python3
"""
Replacement for Figure 4: per cell type, what fraction of cells is Specialist / Generalist / Transitory / Unassigned.

Usage:
    python fig4_specialist_transitory.py collated_labels.csv --celltype-col Cell_type --label-col Arch_after_transitory
Input: the lab's collated label file (one row per cell), NOT the expression matrix.
Writes: fig4_proportions.png/.pdf, fig4_counts_table.csv, fig4_per_label_table.csv (appendix)
"""
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("labels")
ap.add_argument("--celltype-col", default="Cell_type")
ap.add_argument("--label-col", default="Arch_after_transitory")
a = ap.parse_args()

df = pd.read_csv(a.labels)
for c in (a.celltype_col, a.label_col):
    if c not in df.columns:
        raise SystemExit(f"column '{c}' not found. Columns: {list(df.columns)}")

def cls(x):
    x = str(x)
    if x.startswith("Archetype"):  return "Specialist"
    if x.startswith("Generalist"): return "Generalist"
    if x.startswith("Transitory"): return "Transitory"
    return "Unassigned"

df["class"] = df[a.label_col].map(cls)
order_cls = ["Specialist", "Generalist", "Transitory", "Unassigned"]
colors = ["#1f4e79", "#7f7f7f", "#e0a458", "#d9d9d9"]
tab = pd.crosstab(df[a.celltype_col], df["class"]).reindex(columns=order_cls, fill_value=0)
prop = tab.div(tab.sum(axis=1), axis=0) * 100
tab.assign(total=tab.sum(axis=1)).to_csv("fig4_counts_table.csv")
pd.crosstab(df[a.label_col], df[a.celltype_col]).to_csv("fig4_per_label_table.csv")

fig, ax = plt.subplots(figsize=(7.5, 4.8))
types = list(prop.index)
x = range(len(types))
bottom = pd.Series(0.0, index=prop.index)
for c, col in zip(order_cls, colors):
    ax.bar(x, prop[c].values, bottom=bottom.values, color=col, label=c, width=0.6)
    for i, ct in enumerate(types):
        n, p_ = tab.loc[ct, c], prop.loc[ct, c]
        if n == 0:
            continue
        mid = bottom[ct] + p_ / 2
        txt = f"{n} ({p_:.1f}%)"
        if p_ >= 5:   # fits inside the segment
            ax.text(i, mid, txt, ha="center", va="center", fontsize=8,
                    color="white" if c == "Specialist" else "black")
        else:         # small segment: label to the right with a short connector
            ax.annotate(txt, xy=(i + 0.3, mid), xytext=(i + 0.38, mid),
                        ha="left", va="center", fontsize=7.5,
                        arrowprops=dict(arrowstyle="-", lw=0.6, color="0.4"))
    bottom += prop[c]
ax.set_xticks(list(x))
ax.set_xticklabels([f"{t}\n(n = {int(tab.loc[t].sum()):,})" for t in types])
ax.set_xlim(-0.5, len(types) - 0.5 + 0.35)
ax.set_ylim(0, 100)
ax.set_ylabel("% of cells of that type")
ax.set_title("EBUS_06: Pareto label class by cell type")
ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
fig.tight_layout()
fig.savefig("fig4_proportions.png", dpi=300)
fig.savefig("fig4_proportions.pdf")
print(tab.assign(total=tab.sum(axis=1)))
