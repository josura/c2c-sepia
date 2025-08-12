#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# -----------------------
# CONFIG
# -----------------------
# Fitted matrix exported by your previous script (set SAVE_CSV=True there)
FITTED_TSV = "/tmp/AT1_metabolites_iterationMatrix_fitted.tsv"

# metabolite files (scFEA outputs)
metabolites_1hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_1h_metabolites_module168_cell1646_20241014-125146.csv"
metabolites_6hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_6h_metabolites_module168_cell1037_20241014-125629.csv"
metabolites_7hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_7h_metabolites_module168_cell1722_20241014-123944.csv"
metabolites_10hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_10h_metabolites_module168_cell1240_20241014-171445.csv"

# RNA metadata (celltype per well)
rna_1h_metadataFile  = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-1h-metadata.tsv"
rna_6h_metadataFile  = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-6h-metadata.tsv"
rna_7h_metadataFile  = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-7h-metadata.tsv"
rna_10h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-10h-metadata.tsv"

# name mapping files
metabolites_names_csv = "/home/josura/Projects/ccc/fluxes/scFEA/data/cName_complete_mouse_c70_m168.csv"
namemap_csv = "/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/name_map_full_mouse_fixed.csv"

# figure output
SAVE_HTML = False
OUTPUT_DIR = "/tmp/at1_plots"  # created if saving

# -----------------------
# HELPERS
# -----------------------
def prep_metabolites_avg(met_file, meta_file):
    meta = pd.read_csv(meta_file, sep="\t", index_col=0)
    well_celltype = pd.DataFrame({"celltype": meta["cell_type"], "well": meta.index})
    mets = pd.read_csv(met_file, sep=",", index_col=0)
    mets["well"] = mets.index
    mets = mets.merge(well_celltype, on="well", how="left")
    mets.index = mets["well"]
    mets = mets.drop("well", axis=1)
    return mets.groupby(["celltype"]).mean()

def rename_metabolite_columns_to_original(met_1h_avg, met_6h_avg, met_7h_avg, met_10h_avg,
                                          metabolites_names_csv, namemap_csv):
    # Align columns to canonical cNames
    cname = pd.read_csv(metabolites_names_csv, sep=",")
    if len(met_1h_avg.columns) != len(cname.columns):
        raise ValueError("Metabolite columns count does not match cName file columns.")
    newcols = cname.values[0]
    for df in (met_1h_avg, met_6h_avg, met_7h_avg, met_10h_avg):
        df.columns = newcols

    # Map cName Query -> original Match
    namemap = pd.read_csv(namemap_csv, sep=",")
    original_cols = []
    for q in met_1h_avg.columns:
        m = namemap.loc[namemap["Query"] == q, "Match"]
        original_cols.append(m.values[0] if len(m) > 0 else q)
    for df in (met_1h_avg, met_6h_avg, met_7h_avg, met_10h_avg):
        df.columns = original_cols
    return met_1h_avg, met_6h_avg, met_7h_avg, met_10h_avg, namemap

def build_hmdb_to_original(namemap_df, hmdb_cols):
    nm = namemap_df.dropna(subset=['HMDB'])
    mapping = {}
    for hmdb in hmdb_cols:
        m = nm.loc[nm["HMDB"] == hmdb, "Match"]
        mapping[hmdb] = m.values[0] if len(m) > 0 else hmdb
    return mapping

def batched(iterable, n):
    for i in range(0, len(iterable), n):
        yield iterable[i:i+n]

# -----------------------
# MAIN
# -----------------------
if __name__ == "__main__":
    # 1) load fitted AT1 timeseries
    # Expect columns = HMDB IDs, 'time' column present
    fitted_df = pd.read_csv(FITTED_TSV, sep="\t")
    if "time" not in fitted_df.columns:
        raise ValueError("Fitted TSV must contain a 'time' column.")
    fitted_df = fitted_df.set_index("time")
    fitted_df.index = fitted_df.index.astype(float)
    fitted_df.sort_index(inplace=True)

    # 2) load real data (averaged by celltype)
    met_1h_avg  = prep_metabolites_avg(metabolites_1hFile,  rna_1h_metadataFile)
    met_6h_avg  = prep_metabolites_avg(metabolites_6hFile,  rna_6h_metadataFile)
    met_7h_avg  = prep_metabolites_avg(metabolites_7hFile,  rna_7h_metadataFile)
    met_10h_avg = prep_metabolites_avg(metabolites_10hFile, rna_10h_metadataFile)

    # 3) rename columns to original names for real dataframes
    met_1h_avg, met_6h_avg, met_7h_avg, met_10h_avg, namemap_df = rename_metabolite_columns_to_original(
        met_1h_avg, met_6h_avg, met_7h_avg, met_10h_avg, metabolites_names_csv, namemap_csv
    )

    # 4) build HMDB -> original name mapping for plot titles/lookup
    hmdb_to_original = build_hmdb_to_original(namemap_df, fitted_df.columns)

    # 5) pull AT1 row for each timepoint (real data)
    have_AT1 = all('AT1' in df.index for df in (met_1h_avg, met_6h_avg, met_7h_avg, met_10h_avg))
    if not have_AT1:
        print("Warning: AT1 not found in one or more real datasets; some markers will be missing.")

    # 6) plot in 3x3 subplots per page
    hmdb_cols = list(fitted_df.columns)

    if SAVE_HTML and not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    page = 1
    for cols9 in batched(hmdb_cols, 9):
        fig = make_subplots(rows=3, cols=3, subplot_titles=[
            hmdb_to_original.get(c, c) for c in cols9
        ])

        row = col = 1
        for hmdb in cols9:
            original_name = hmdb_to_original.get(hmdb, hmdb)

            # Skip virtuals
            if ("v-in" in original_name) or ("v-out" in original_name):
                # still add empty subplot title to keep grid consistent
                trace_x = fitted_df.index.values
                trace_y = fitted_df[hmdb].values
                fig.add_trace(go.Scatter(x=trace_x, y=trace_y, mode="lines",
                                         name=f"{original_name} (fitted)", showlegend=False),
                              row=row, col=col)
            else:
                # Fitted line
                fig.add_trace(
                    go.Scatter(
                        x=fitted_df.index.values,
                        y=fitted_df[hmdb].values,
                        mode="lines",
                        name=f"{original_name} (fitted)"
                    ),
                    row=row, col=col
                )

                # Real markers at 0,5,6,9 mapping to 1h,6h,7h,10h
                # We need to look up by original_name in the averaged real tables
                xs = []
                ys = []
                # 0 -> 1h
                if have_AT1 and original_name in met_1h_avg.columns:
                    xs.append(0.0)
                    ys.append(float(met_1h_avg.loc['AT1', original_name]))
                # 5 -> 6h
                if have_AT1 and original_name in met_6h_avg.columns:
                    xs.append(5.0)
                    ys.append(float(met_6h_avg.loc['AT1', original_name]))
                # 6 -> 7h
                if have_AT1 and original_name in met_7h_avg.columns:
                    xs.append(6.0)
                    ys.append(float(met_7h_avg.loc['AT1', original_name]))
                # 9 -> 10h
                if have_AT1 and original_name in met_10h_avg.columns:
                    xs.append(9.0)
                    ys.append(float(met_10h_avg.loc['AT1', original_name]))

                if len(xs) > 0:
                    fig.add_trace(
                        go.Scatter(
                            x=xs, y=ys, mode="markers",
                            name=f"{original_name} (real)"
                        ),
                        row=row, col=col
                    )

            # advance grid position
            if col == 3:
                row += 1
                col = 1
            else:
                col += 1

        fig.update_layout(
            height=900, width=1200,
            title_text="AT1 Metabolites — Fitted Timeseries vs Real Points (1h, 6h, 7h, 10h)"
        )

        if SAVE_HTML:
            out_path = os.path.join(OUTPUT_DIR, f"AT1_fitted_vs_real_page{page}.html")
            fig.write_html(out_path)
            print(f"Saved: {out_path}")
        else:
            fig.show()

        page += 1
