#!/usr/bin/env python3
"""
02_prepare_cell_table.py
Build a cells x circuit-node table for the EBUS_06 myeloid groups defined by the lab's Pareto labels.

Inputs (same folder):  Collated_EBUS_06.csv   (labels: cell_type, Arch_after_transitory)
                       EBUS_06.csv            (full matrix, genes x cells, log-normalised; unzip EBUS_06.zip first)
Output:                myeloid_node_table_EBUS06.csv   (one row per cell: group + node values)
                       node_detection_EBUS06.csv       (% of cells in each group with the gene detected)
Group mapping: Archetype4 -> Specialist, Generalist -> Generalist, Transitory_4_5 -> Transitory
(Archetype3, Transitory_1_4 and unlabeled NA cells are excluded, see the notes in the chat).
"""
import pandas as pd, numpy as np

NODES = {  # circuit node -> gene(s) in the matrix (mean of the log-normalised values if several)
    "STAT1": ["STAT1"], "STAT3": ["STAT3"], "STAT6": ["STAT6"], "SOCS1": ["SOCS1"], "SOCS3": ["SOCS3"],
    "NFKB": ["NFKB1", "RELA"], "AP1": ["JUN", "FOS"], "TNFA": ["TNF"], "IL10": ["IL10"], "TGFB": ["TGFB1"],
    "HIF1A": ["HIF1A"], "SMAD23": ["SMAD2", "SMAD3"], "ERK": ["MAPK1", "MAPK3"], "TLR4": ["TLR4"],
}
GROUPS = {"Archetype4": "Specialist", "Generalist": "Generalist", "Transitory_4_5": "Transitory"}
genes = sorted({g for v in NODES.values() for g in v})

col = pd.read_csv("Collated_EBUS_06.csv", index_col=0)
m = col[(col.cell_type == "Myeloid") & col.Arch_after_transitory.isin(GROUPS)].copy()
m["group"] = m.Arch_after_transitory.map(GROUPS)

rows = {}
for ch in pd.read_csv("EBUS_06.csv", index_col=0, chunksize=2000):   # genes x cells, read in chunks (low memory)
    for g in genes:
        if g in ch.index: rows[g] = ch.loc[g]
E = pd.DataFrame(rows).loc[m.index]
missing = [g for g in genes if g not in E.columns]
assert not missing, f"genes missing from matrix: {missing}"

det = (E > 0).groupby(m.group).mean().T * 100
det.round(1).to_csv("node_detection_EBUS06.csv")
tab = pd.DataFrame({n: E[g].mean(axis=1) for n, g in NODES.items()})
tab.insert(0, "group", m.group)
tab.to_csv("myeloid_node_table_EBUS06.csv")
print(m.group.value_counts().to_dict())
print((det.round(0)).to_string())
