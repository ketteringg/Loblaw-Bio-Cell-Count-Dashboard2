"""
load_data.py

Initializes the SQLite database (cell_counts.db) using schema.sql and loads
all rows from cell-count.csv into it.

Run directly:
    python load_data.py

No CLI arguments or module-style execution (`python -m`) required or supported.
Paths resolve relative to this script's own location, so it can be invoked
from any working directory; cell-count.csv and schema.sql just need to sit
in the same directory as this script.
"""
import sqlite3
from pathlib import Path

import pandas as pd

# Paths anchored to this script's location, not the CWD (see module docstring).
ROOT = Path(__file__).parent
CSV_PATH = ROOT / "cell-count.csv"
SCHEMA_PATH = ROOT / "schema.sql"
DB_PATH = ROOT / "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

REQUIRED_COLUMNS = [
    "project", "subject", "condition", "age", "sex", "treatment",
    "response", "sample", "sample_type", "time_from_treatment_start",
    *POPULATIONS,
]
NON_NULLABLE_COLUMNS = [c for c in REQUIRED_COLUMNS if c != "response"]
EXPECTED_VALUES = {
    "sex": {"M", "F"},
    "treatment": {"miraclib", "phauximab", "none"},
    "sample_type": {"PBMC", "WB"},
    "response": {"yes", "no"},  # nulls checked separately
}


class DataValidationError(Exception):
    """Raised for data problems severe enough to stop the load entirely."""


def validate(df: pd.DataFrame) -> None:
    """
    Guards against a non-conforming CSV (missing columns, negative counts,
    nulls in required fields, duplicate IDs, non-numeric counts, response/
    treatment inconsistency). Every check below runs independently and adds
    its own message to `errors` -- a check is only skipped if it genuinely
    can't run (e.g. checking nulls in a column that doesn't exist), never
    because an earlier check already failed. This way a CSV with multiple
    simultaneous problems surfaces ALL of them in one pass, rather than
    reporting the first failure, getting fixed, rerun, and hitting the next
    one one at a time.

    Hard failures (in `errors`) raise DataValidationError, listing every
    failure found, before any DB writes happen. Softer issues (in
    `warnings`) are printed but don't stop the load, since they may reflect
    a legitimate new category rather than corrupt data.
    """
    errors = []
    warnings = []

    # --- structural: required columns present ---
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        errors.append(f"Missing required column(s): {missing_cols}")
    present_cols = set(df.columns)

    # --- duplicate sample IDs ---
    if "sample" in present_cols:
        dupe_samples = df["sample"][df["sample"].duplicated()].unique()
        if len(dupe_samples) > 0:
            errors.append(f"Duplicate sample id(s): {list(dupe_samples)[:5]}"
                           f"{' ...' if len(dupe_samples) > 5 else ''}")
    else:
        errors.append("Skipped duplicate-sample-id check: 'sample' column is missing.")

    # --- nulls in required (non-nullable) fields ---
    for col in NON_NULLABLE_COLUMNS:
        if col not in present_cols:
            continue  # already flagged by the missing-columns check above
        n_null = df[col].isna().sum()
        if n_null > 0:
            errors.append(f"Column '{col}' has {n_null} null value(s) but is required.")

    # --- cell count columns: numeric and non-negative ---
    for col in POPULATIONS:
        if col not in present_cols:
            continue  # already flagged by the missing-columns check above
        non_numeric = pd.to_numeric(df[col], errors="coerce").isna() & df[col].notna()
        if non_numeric.any():
            errors.append(f"Column '{col}' has {non_numeric.sum()} non-numeric value(s).")
            continue  # skip the negativity check on this column, it's not numeric
        negative = pd.to_numeric(df[col]) < 0
        if negative.any():
            errors.append(f"Column '{col}' has {negative.sum()} negative value(s).")

    # --- response only null when treatment == 'none' ---
    if {"response", "treatment"} <= present_cols:
        bad_null_response = df[df["response"].isna() & (df["treatment"] != "none")]
        if not bad_null_response.empty:
            errors.append(
                f"{len(bad_null_response)} row(s) have null response but a real "
                "treatment (expected null response only when treatment='none')."
            )
        bad_nonnull_response = df[df["response"].notna() & (df["treatment"] == "none")]
        if not bad_nonnull_response.empty:
            errors.append(
                f"{len(bad_nonnull_response)} row(s) have treatment='none' but a "
                "non-null response (untreated subjects shouldn't have a response)."
            )
    else:
        missing = {"response", "treatment"} - present_cols
        errors.append(f"Skipped response/treatment consistency check: missing {missing}.")

    # --- categorical values within expected sets (soft warning, not fatal) ---
    for col, allowed in EXPECTED_VALUES.items():
        if col not in present_cols:
            continue  # already flagged by the missing-columns check above
        actual = set(df[col].dropna().unique())
        unexpected = actual - allowed
        if unexpected:
            warnings.append(f"Column '{col}' has unexpected value(s) not in {allowed}: {unexpected}")

    # --- age sanity range (soft warning; not a hard business rule) ---
    if "age" in present_cols:
        implausible_age = df[(df["age"] < 0) | (df["age"] > 120)]
        if not implausible_age.empty:
            warnings.append(f"{len(implausible_age)} row(s) have an implausible age (<0 or >120).")

    for w in warnings:
        print(f"WARNING: {w}")

    if errors:
        raise DataValidationError(
            f"CSV failed validation with {len(errors)} issue(s), aborting load:\n  - "
            + "\n  - ".join(errors)
        )


def build_database() -> None:
    # Start from a clean slate every run. The dataset is small (10.5k rows),
    # so a full rebuild is fast and guarantees the .db never silently drifts
    # out of sync with cell-count.csv (e.g. if the CSV is swapped out).
    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    try:
        # Enforce foreign key constraints -- off by default in SQLite. With
        # this on, an insert referencing a nonexistent subject_id/sample_id
        # will raise immediately instead of silently corrupting the schema's
        # relational integrity.
        conn.execute("PRAGMA foreign_keys = ON")

        with open(SCHEMA_PATH) as f:
            conn.executescript(f.read())

        df = pd.read_csv(CSV_PATH)
        validate(df)  # raises DataValidationError and aborts before any writes if non-conforming

        # --- subjects: one row per subject, deduplicated ---
        # Verified against the source data: every subject has exactly one
        # project/condition/age/sex/treatment/response across all their
        # samples, so drop_duplicates on subject is safe here.
        subjects = df.drop_duplicates(subset="subject")[
            ["subject", "project", "condition", "age", "sex", "treatment", "response"]
        ].copy()
        subjects.columns = [
            "subject_id", "project", "condition", "age", "sex", "treatment", "response",
        ]
        # pandas reads a column with any NaN as float64; make sure age (which
        # has no NaNs) stays an int rather than silently becoming 64.0 etc.
        subjects["age"] = subjects["age"].astype(int)
        subjects.to_sql("subjects", conn, if_exists="append", index=False)

        # --- samples: one row per physical sample ---
        samples = df[["sample", "subject", "sample_type", "time_from_treatment_start"]].copy()
        samples.columns = ["sample_id", "subject_id", "sample_type", "time_from_treatment_start"]
        samples["time_from_treatment_start"] = samples["time_from_treatment_start"].astype(int)
        samples.to_sql("samples", conn, if_exists="append", index=False)

        # --- cell_counts: wide (5 population columns) -> long (1 row each) ---
        cell_counts = df.melt(
            id_vars=["sample"],
            value_vars=POPULATIONS,
            var_name="population",
            value_name="count",
        )
        cell_counts.columns = ["sample_id", "population", "count"]
        cell_counts["count"] = cell_counts["count"].astype(int)
        cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)

        conn.commit()

        # --- sanity check: confirm row counts match expectations ---
        n_subjects = conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
        n_samples = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        n_cell_counts = conn.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0]

        print(f"Database built at {DB_PATH}")
        print(f"  subjects:    {n_subjects}")
        print(f"  samples:     {n_samples}")
        print(f"  cell_counts: {n_cell_counts} (samples x {len(POPULATIONS)} populations)")

        assert n_samples == len(df), "sample count mismatch vs source CSV"
        assert n_cell_counts == len(df) * len(POPULATIONS), "cell_counts row count mismatch"

    finally:
        conn.close()


if __name__ == "__main__":
    build_database()
