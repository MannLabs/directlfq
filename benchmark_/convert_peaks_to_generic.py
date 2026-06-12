"""Convert the PEAKS DIA feature export (HYE1-6_repeat/lfq.dia.features.csv) into
the directLFQ generic input format.

PEAKS is not one of directLFQ's auto-detected input types, so we reshape it into
the generic wide format: columns "protein", "ion", then one column per run holding
the ion intensity. Saving with the ".aq_reformat.tsv" suffix makes directLFQ treat
it as a pre-formatted generic table.
"""

import pandas as pd

INPUT_CSV = "HYE1-6_repeat/lfq.dia.features.csv"
OUTPUT_TSV = "HYE1-6_repeat/lfq.dia.features.csv.aq_reformat.tsv"

PROTEIN_SOURCE_COL = "Accession"
PEPTIDE_COL = "Peptide"
CHARGE_COL = "z"
INTENSITY_SUFFIX = " Normalized Area"

# Aggregate per-condition columns in the PEAKS export that also end in the
# intensity suffix but are not per-run sample columns.
AGGREGATE_PREFIX = "Group "


def _sample_columns(columns: list[str]) -> list[str]:
    return [
        c
        for c in columns
        if c.endswith(INTENSITY_SUFFIX) and not c.startswith(AGGREGATE_PREFIX)
    ]


def main() -> None:
    df = pd.read_csv(INPUT_CSV, encoding="latin1")

    sample_cols = _sample_columns(df.columns.tolist())
    print(f"rows in: {len(df)}  sample columns: {len(sample_cols)}")

    out = pd.DataFrame()
    out["protein"] = df[PROTEIN_SOURCE_COL]
    out["ion"] = df[PEPTIDE_COL].astype(str) + "_" + df[CHARGE_COL].astype(str)

    # Run-name headers strip the trailing " Normalized Area" suffix.
    for col in sample_cols:
        run_id = col[: -len(INTENSITY_SUFFIX)]
        out[run_id] = df[col]

    n_dupe = out["ion"].duplicated().sum()
    if n_dupe:
        print(f"WARNING: {n_dupe} duplicate ion ids (Peptide+z); directLFQ will drop these")

    out = out.dropna(subset=["protein"])
    out.to_csv(OUTPUT_TSV, sep="\t", index=False)
    print(f"wrote {len(out)} rows -> {OUTPUT_TSV}")


if __name__ == "__main__":
    main()