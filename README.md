# Loblaw Bio Immune Cell Analysis

Analysis pipeline and interactive dashboard for Bob Loblaw's clinical trial
data, examining how immune cell populations relate to treatment response.

**Live dashboard:** https://loblaw-bio-cell-count-dashboard2-gyjnefkvjgzi5kkkv7wekn.streamlit.app/

Deployed on Streamlit Community Cloud. It builds `cell_counts.db`
automatically from `cell-count.csv` on first load (see "Self-initializing
database" below), so no setup step is needed to view it. For local
development, testing, or to reproduce the analysis yourself, see the
Quickstart below.

## Quickstart (GitHub Codespaces)

```bash
make setup       # installs dependencies (uses uv if available, falls back to pip)
make pipeline    # builds cell_counts.db, then generates every Part 2-4 output file/plot
make dashboard   # streamlit run app.py, launches the interactive dashboard
make test        # installs dev dependencies, builds cell_counts.db if missing, then runs pytest
```

Equivalent plain commands, if you would rather not use `make`:

```bash
pip install -r requirements.txt
python load_data.py
python generate_outputs.py
streamlit run app.py
```

`make pipeline` runs `load_data.py` (builds `cell_counts.db` from
`cell-count.csv`) followed by `generate_outputs.py`, which writes
`part2_frequency_table.csv`, `part3_stats_results.csv`,
`part3_boxplot_responders.png`, `part4_baseline_melanoma_samples.csv`, and
`part4_summary.txt`. See "Repository contents" below for what each one is.
You do not strictly have to run either first before the dashboard, though:
`app.py`
self-initializes the database on its own if it is missing (see
"Self-initializing database" below), which is what makes the Streamlit
Community Cloud deployment work without a separate setup step.

Tested with Python 3.12 and the exact package versions pinned in
`requirements.txt`.

## Self-initializing database

`app.py` checks for `cell_counts.db` on startup. If it's missing but
`cell-count.csv` is present, it builds the database automatically (calling
the same `load_data.build_database()` function `make pipeline` uses)
before rendering anything, showing a brief spinner while it does. If
neither file is present, it shows a clear error instead of crashing.

This matters specifically for cloud deployment: platforms like Streamlit
Community Cloud start `app.py` directly and have no built-in "run this
setup script first" step. Locally, this path essentially never fires,
since `make pipeline` / `python load_data.py` already builds the database
ahead of time.

## Tests

```bash
make test
```

runs the full test suite (`tests/test_analysis.py` and `tests/test_app.py`,
83 tests as of writing) via `pytest`. This also runs automatically on every
push and pull request via GitHub Actions (`.github/workflows/tests.yml`).
Check the "Actions" tab on the repo, or the checkmark next to any commit,
to see the result without running anything locally.

`make test` is self-sufficient: it installs dev dependencies and, if
`cell_counts.db` does not already exist, builds it first, via a real Make
file dependency on `cell-count.csv`/`schema.sql`/`load_data.py` (not just
an unconditional rebuild every time). Neither `make setup` nor
`make pipeline` needs to be run first. This is a second, independent
layer alongside `tests/conftest.py`'s own `ensure_db` fixture, which
already builds the database itself whenever pytest is run directly
(`.github/workflows/tests.yml`'s CI job does exactly this, calling
`pytest -v` directly rather than through `make test`).

`test_analysis.py` also includes tests that specifically guard against
committed output files going stale. `test_no_stale_prerename_output_
files_exist` checks that none of the old, pre-rename filenames
(`frequency_table.csv`, `stats_results.csv`, `boxplot_responders.png`,
`baseline_melanoma_samples.csv`) exist anywhere in the repo, full stop,
regardless of content. `test_committed_part3_stats_results_matches_
fresh_regeneration` and two neighbors instead check the CORRECTLY-named
committed files (`part3_stats_results.csv` etc.): each reads whatever is
actually committed right now and compares it against a fresh
recomputation from the current code, failing loudly if they disagree,
rather than trusting that whoever last changed the analysis also
remembered to regenerate and re-commit the output files by hand.

Both checks exist because a content-only check turned out not to be
enough. A stale, pre-melanoma-fix `stats_results.csv`, showing a
different cohort size and incorrectly flagging `b_cell` and `monocyte`
as significant, was found committed under its old filename, sitting
alongside the correct `part3_stats_results.csv`. Content-checking only
the correctly-named file caught nothing, because nothing was checking
whether the old-named duplicate existed at all: that old file survived,
undetected, across more than one submission before it was found by
manual inspection rather than by the test suite. Confirmed directly, for
both tests, that they actually catch what they are meant to: reintroducing
the exact stale file under its old name makes
`test_no_stale_prerename_output_files_exist` fail immediately, and
reintroducing it under the current, correct filename makes
`test_committed_part3_stats_results_matches_fresh_regeneration` fail
with a clear message naming the discrepancy. Both pass again once the
old file is removed and the current one is freshly regenerated.

`test_analysis.py` covers the query/analysis layer directly: known row
counts, the four-case stats handling (no samples, no response data, no
individuals in a group, small-n warning), the invariant that a population
filter doesn't change what percentages are computed against, specific
cross-checked numeric results (including the assignment's graded form
answer, as a regression check), the N-group pairwise comparison and paired
population-vs-population comparison functions, and `load_data.py`'s
validation logic (both the hard failures and the soft warnings).

`test_app.py` uses Streamlit's own `AppTest` to run the dashboard headlessly
and check for a handful of representative scenarios per tab (clean
load, the small-n warning path, the reset-filters button, each comparison
mode, mode switching). These are the same checks that were run manually
throughout development. More than one caught a real bug before it
shipped, now made permanent.

## Code structure and why

Three layers, each with one job:

- **`load_data.py`** owns getting the raw CSV into a validated, relational
  form. It's the only place that touches `cell-count.csv` directly or knows
  the wide-format-to-long-format conversion. Nothing downstream re-reads the
  CSV.
- **`analysis.py`** owns every query and statistical computation, and knows
  nothing about Streamlit or the dashboard. It's plain functions that take a
  connection or a DataFrame and return a DataFrame. That is what makes it
  possible to unit-test the actual statistics (`tests/test_analysis.py`)
  without spinning up the UI, and what `generate_outputs.py` calls directly
  to produce the graded Part 2-4 files without going through the dashboard
  at all.
- **`app.py`** owns presentation only: it calls `analysis.py` functions and
  renders the results. It's deliberately general-purpose (any cohort, any
  comparison, via filters) rather than having one fixed view per assignment
  part, because Parts 2-4 turn out to be special cases of "filter a cohort,
  then optionally compare groups within it," the thing the dashboard
  already does for anything.

That separation is also why there's a **`generate_outputs.py`**, distinct
from the dashboard. The assignment's outputs need to exist as concrete
files a grader can open without clicking through the UI to reconstruct the
exact required filter combination, so this script calls the same
`analysis.py` functions directly and writes them to disk. Both `app.py` and
`generate_outputs.py` are two different callers of the same underlying
logic, not two separately maintained copies of it: a change to
`analysis.py`'s statistics changes both consistently.

## Repository contents

| File | Purpose |
|---|---|
| `Makefile` | `make setup` / `make pipeline` / `make dashboard` / `make test`. `pipeline` builds the database and generates every required Part 2-4 output file/plot, in that order, with no manual steps in between |
| `.streamlit/config.toml` | Dashboard theme (colors, font) |
| `.github/workflows/tests.yml` | Runs the test suite on every push/PR |
| `.gitignore` | Editor/OS/venv ignores, plus `cell_counts.db` (see below) |
| `schema.sql` | Relational schema (3 tables) |
| `load_data.py` | Validates `cell-count.csv`, builds `cell_counts.db` |
| `analysis.py` | Query + analysis functions (framework agnostic, no Streamlit) |
| `generate_outputs.py` | Produces the required Part 2-4 files/plot by calling `analysis.py` directly (see "Code structure and why" above) |
| `app.py` | Streamlit dashboard: 5 tabs, single-cohort exploration plus 4 comparison modes |
| `tests/` | pytest suite (`test_analysis.py`, `test_app.py`, `conftest.py`) |
| `cell-count.csv` | Source data |
| `cell_counts.db` | Generated by `make pipeline`. Gitignored: it is an intermediate build artifact, not itself one of the assignment's requested outputs, and is quickly and deterministically rebuilt from `cell-count.csv` |
| `part2_frequency_table.csv` | Part 2's required output. Generated by `make pipeline`, and committed (not gitignored), since the assignment asks the submission to include "any input or output files generated." See "Part 2: frequency table" below for how `percentage` relates to relative frequency and to `total_count` |
| `part3_stats_results.csv`, `part3_boxplot_responders.png` | Part 3's required outputs. Generated by `make pipeline` and committed, same reasoning as above |
| `part4_baseline_melanoma_samples.csv`, `part4_summary.txt` | Part 4's required outputs. Generated by `make pipeline` and committed, same reasoning as above |
| `requirements.txt` | Pinned runtime dependencies |
| `requirements-dev.txt` | Runtime dependencies plus `pytest` |

Committing the part2/3/4 files reintroduces a staleness risk that
gitignoring them would have avoided: a committed file can in principle
drift out of sync with the code that produces it. This is a real risk,
not a hypothetical one, and is covered in "Tests" below.

## Data dictionary (`cell-count.csv`)

| Column | Type | Meaning | Valid values |
|---|---|---|---|
| `project` | text | Study/project identifier | `prj1`, `prj2`, `prj3` |
| `subject` | text | Unique subject/patient identifier | e.g. `sbj000` |
| `condition` | text | Subject's diagnosis | `melanoma`, `carcinoma`, `healthy` |
| `age` | integer | Subject's age in years | 50-79 in this dataset |
| `sex` | text | Subject's sex | `M`, `F` |
| `treatment` | text | Treatment arm | `miraclib`, `phauximab`, `none` (untreated, including all healthy subjects) |
| `response` | text | Treatment response | `yes`, `no`, or null (null exactly when `treatment='none'`, since response doesn't apply to untreated subjects) |
| `sample` | text | Unique sample identifier | e.g. `sample00000` |
| `sample_type` | text | Tissue/specimen type | `PBMC`, `WB` (whole blood) |
| `time_from_treatment_start` | integer | Days since treatment start at time of sample draw | `0`, `7`, `14` |
| `b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, `monocyte` | integer | Raw cell count for that population in this sample | non-negative integers |

Each subject has exactly one `project` / `condition` / `age` / `sex` /
`treatment` / `response` value across all their samples (verified against
the data), and exactly 3 samples (one per timepoint). See "Schema design"
below for how this gets normalized into `subjects` / `samples` /
`cell_counts` tables.

## Dashboard architecture: one page, 5 tabs, each independently filterable

Earlier versions of this dashboard used tabs differently. First a fixed
tab per Part 2/3/4 cohort, then two general-purpose tabs (Custom Explorer
and Cohort Comparison), then briefly a single page with a radio-button
mode selector instead of tabs at all. Each was replaced for a different
reason: the fixed tabs duplicated what a more general view could already
show; the two-tab version couldn't express every comparison someone might
want; and the radio-selector version had a real bug: switching away
from a view and back silently cleared whatever filters had been set,
because Streamlit only instantiates the currently-selected branch's
widgets under an `if/elif`, and drops `session_state` for any widget that
stops being instantiated. **Tabs render every tab's content on every
rerun** (just CSS-hiding the inactive ones), so widget state survives
switching tabs for free. Confirmed directly: selecting specific
populations in the By Population tab, then interacting with a completely
different tab, no longer clears that selection (`test_app.py`'s
`test_selections_persist_across_tabs`).

The dashboard is 5 tabs: **Default**, **Responder vs Non-responder**, **By
Population**, **By Date**, and **Custom**. Default is single-cohort
exploration (cohort summary, average cell counts, the frequency table,
and a per-population distribution boxplot, with no comparison). The other
four each compare 2 or more groups against each other, split a different
way:

- **Responder vs Non-responder**: always exactly 2 groups, the
  assignment's Part 3 axis.
- **By Population**: select 2+ cell populations (all 5 by default) and
  compare them directly against each other *within* the same cohort.
  This uses a **paired test** (Wilcoxon signed-rank), not the unpaired
  Mann-Whitney used everywhere else, because two populations'
  percentages from the same sample aren't independent samples. See
  "Statistical approach" below for why that distinction matters.
- **By Date**: select 2 or all 3 timepoints and compare them against
  each other, per population. Selecting 3 timepoints tests every pair
  (3 pairs), not just one.
- **Custom**: build 2 to 4 independent cohorts with any combination of
  filters and compare them all against each other.

Any number of groups beyond 2 is handled by testing every pairwise
combination, Bonferroni-corrected across all pairs and populations tested
together (see `compare_n_groups` in `analysis.py`). This correction gets
stricter fast as more groups are added: 3 groups is 3 pairs per
population, 4 groups is 6 pairs per population. That's the real
statistical cost of choosing pairwise tests over a single omnibus test
(e.g. Kruskal-Wallis), worth knowing before comparing many groups at
once.

Responder/Non-responder and Cohort A/B deliberately share the exact same
2 colors (blue / reddish purple, Okabe-Ito): both are "first group vs.
second group in a two-way split," so one consistent color pairing is used
throughout rather than two different ones.

The required, graded Part 2-4 answers can still be reproduced exactly:
Part 2 and Part 4's baseline summary from the Default tab's frequency
table and cohort summary with the right filters applied, and Part 3 from
the Responder vs Non-responder tab with the melanoma, miraclib, and PBMC
filters applied (documented in the assignment spec).

## Schema design

```
subjects (subject_id PK, project, condition, age, sex, treatment, response)
samples  (sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
cell_counts (sample_id FK, population, count, PRIMARY KEY(sample_id, population))
```

**Rationale:**

- **`subjects`** holds attributes that are fixed per person. Verified against
  the source data: every subject has exactly one `project` / `condition` /
  `age` / `sex` / `treatment` / `response` value across all their samples, so
  storing these once (rather than repeating them on every sample row) avoids
  redundancy and keeps updates consistent.
- **`samples`** is one row per physical sample draw, the unit that varies by
  timepoint (`time_from_treatment_start`) and `sample_type` (PBMC/WB).
- **`cell_counts`** is stored in long format (one row per population per
  sample) rather than as five wide columns. This makes "group by population"
  queries (Part 2's frequency table, Part 3's per-population stats) a plain
  SQL/pandas group-by instead of manually unpivoting five columns in code,
  and it's the standard normalized shape for repeated-measurement data like
  this.
- **`response` is nullable by design**, not an error state. It is `NULL`
  exactly when `treatment = 'none'`. Untreated (including healthy) subjects
  have no response outcome to record. `load_data.py` validates this
  relationship explicitly (in both directions) before writing any rows, so a
  future CSV that breaks this assumption fails loudly instead of silently
  loading bad data.
- **Foreign key constraints are enforced** (`PRAGMA foreign_keys = ON`), so
  the relational integrity in the schema is actually checked at insert time,
  not just declared as documentation.

**Scalability.** This schema and SQLite specifically hold up fine at the
current size (3,500 subjects, 10,500 samples, 52,500 cell-count rows) and
would still be comfortable at 10-100x that. Past that, three things would
actually need to change:

- **`project` becomes its own table.** Right now it's a plain string column
  on `subjects`. At hundreds of projects, project-level metadata (PI,
  institution, funding source, start date) needs somewhere to live that
  isn't duplicated across every subject row: a `projects` table
  (`project_id` PK, plus those attributes) with `subjects.project_id` as a
  foreign key, same normalization logic already applied to `subjects` and
  `samples`.
- **Indexes on the columns actually filtered on.** `subject_id`,
  `sample_type`, `time_from_treatment_start`, and `population` are the
  columns every query in this project group-bys or filters on. SQLite
  handles the current volume fine without extra indexes, but at real scale
  these would need explicit indexes (or, in Postgres, the equivalent).
  Otherwise every filtered query degrades to a full table scan as row
  counts grow.
- **SQLite itself becomes the limiting factor before the schema does.**
  SQLite is single-writer and file-based, which is exactly right for a
  self-contained, reproducible analysis pipeline like this one, but the
  wrong choice once multiple analysts need concurrent read/write access or
  the database needs to live somewhere other than "next to the code that
  built it." The schema itself is standard relational SQL and moves to
  Postgres (or similar) without any redesign. This is a deployment
  decision, not a schema one.

**For "various types of analytics"**, the shape that generalizes well is
already in place: `cell_counts` being long-format, keyed by `sample_id` +
`population`, means adding a 6th cell population is just new rows, not a
schema change (a wide, one-column-per-population table would need an
actual `ALTER TABLE` every time a new population is measured). The same
pattern extends to analytics beyond cell counts: a new data type per sample
(e.g. genomic variants, imaging features, longitudinal lab values) is a new
fact table keyed by `sample_id`, sitting alongside `cell_counts` rather than
inside it, a star-schema-style extension, not a rewrite of `subjects` or
`samples`.

## Data validation

`load_data.py` validates `cell-count.csv` before writing anything to the
database. Every check is independent (a missing column doesn't suppress
other findings; for example, a missing column and a negative count in an unrelated
column are both reported together), and **all failures found are listed in
one error message**, not just the first one encountered.

**Hard failures** (abort the load, nothing is written to the database):
- Missing required columns
- Negative or non-numeric cell counts
- Nulls in required (non-nullable) fields
- Duplicate sample IDs
- `response` inconsistent with `treatment` (non-null when `treatment='none'`,
  or null when it isn't)

**Soft warnings** (printed, load continues):
- Categorical values outside the expected set (e.g. an unexpected `sex` or
  `treatment` value). This is flagged rather than rejected, since it could
  reflect a legitimate new category rather than corrupt data
- Implausible age values (<0 or >120)

**What this does and does not cover if a different `cell-count.csv` is
used.** The checks above are structural and re-run automatically on
whatever CSV is actually present, so most non-conforming data is caught
without any code changes: a missing column, a negative or non-numeric
cell count, a duplicate sample ID, a null where one shouldn't be, or a
`response`/`treatment` inconsistency all still abort the load and report
every issue found, on any dataset with this same column structure, not
just the current one.

Two things are not automatically handled, and are worth being explicit
about rather than implying more automation exists than actually does:

- **A genuinely new cell population (a 6th population column, not just a
  new value in an existing categorical column) is silently dropped, not
  flagged.** `build_database()`'s wide-to-long conversion is
  `df.melt(..., value_vars=POPULATIONS, ...)`, and `POPULATIONS` is a
  fixed list of the 5 known population names. An extra column beyond
  those 5 is never read, never validated, never loaded, and never
  warned about; it simply doesn't exist as far as the database is
  concerned. Supporting a 6th population is a 2-file code change, not a
  data-driven one: `POPULATIONS` is independently defined in both
  `load_data.py` (kept dependency-light and self-contained on purpose,
  so it doesn't import from `analysis.py`) and `analysis.py` (the
  version everything else, `app.py` and `generate_outputs.py` included,
  imports from). Both would need the new population name added.
- **The test suite's specific expected numbers are pinned to the current
  dataset, not derived generically.** Tests like
  `test_database_has_expected_row_counts` (3,500 subjects / 10,500
  samples / 52,500 cell-count rows) and
  `test_baseline_melanoma_sample_count` (656) assert exact values,
  deliberately: they exist to catch a real computational bug reappearing
  (like the melanoma-filter regression described in "Statistical
  approach" above), which requires checking against a known-correct
  number, not just "some number came out." Swapping in a genuinely
  different `cell-count.csv` (different subject count, different cohort
  composition) would make these specific tests fail, correctly, because
  the ground truth actually changed; they'd need their expected values
  updated to match the new dataset, not just re-run. This is a
  deliberate regression-testing tradeoff, not an oversight: the
  structural checks in `load_data.py` above are what's meant to
  generalize to new data, not these value-pinned tests.

The output-side counterpart to this is `test_no_stale_prerename_output_
files_exist` and its neighbors ("Tests" above): once a new CSV is loaded
and the analysis genuinely changes, those tests are what catch a
committed output file that wasn't regenerated to match.

## Part 2: frequency table

`percentage` and "relative frequency" are the same quantity. The
assignment's own column definition says so directly: `percentage` is
"relative frequency in percentage," meaning the fraction each population
makes up of a sample's total cell count, expressed on a 0 to 100 scale
rather than a 0 to 1 fraction. `get_frequency_table()` computes it exactly
that way:

```
total_count = sum of count across all 5 populations, for that sample
percentage  = count / total_count * 100
```

Verified directly for a real sample: the 5 `percentage` values for any
one sample always sum to 100 (up to floating point precision), confirming
`percentage` genuinely is that population's relative share of the
sample's total, not some other transformed quantity. This is a plain
descriptive statistic. No transform (like the CLR used for the Part 3
significance test, see "Statistical approach" below) is applied here.

`total_count`, the sum across all 5 populations, is the assignment's own
explicit first step ("for each sample, calculate the total number of
cells by summing the counts across all five populations"), before
`percentage` is computed from it. `part2_frequency_table.csv` contains
this value directly, repeated once per population row since the table is
one row per sample-population pair; the 5 repeated values for any one
sample are always identical, since `total_count` is a per-sample
quantity, not a per-population one.

## Statistical approach (Part 3)

Responders vs. non-responders are compared per cell population using the
**Mann-Whitney U test** rather than a t-test. Cell-frequency data is
bounded and often skewed, so a non-parametric test that doesn't assume
normality is the safer default here.

**The test runs on a CLR-transformed value, not the raw percentage.** The
5 cell populations are compositional data: they sum to 100% of each
sample's total cell count, so an increase in one population mechanically
forces the others down. That "closure" constraint means the 5 populations
aren't independent, and testing directly on raw percentages (or raw
counts) treats them as if they were, which can manufacture or mask
apparent significance. The fix is the **centered log-ratio (CLR)
transform**: for each sample, `clr_i = ln(count_i) - mean(ln(counts))`
across that sample's own 5 populations, moving the data into ordinary,
unconstrained real space before testing.

This isn't a hypothetical concern for this dataset, it changes the
answer. Tested empirically, both ways, on the same melanoma, miraclib,
and PBMC cohort that Part 3 asks for (see "Part 3 asks for melanoma
patients specifically" below):

| Population | Raw-percentage test | CLR-based test |
|---|---|---|
| `cd4_t_cell` | not significant (p≈0.0133) | **significant** (p≈0.0025) |
| `b_cell` | not significant (p≈0.0557) | not significant (p≈0.1364) |
| `nk_cell` | not significant | not significant |
| `monocyte` | not significant | not significant |
| `cd8_t_cell` | not significant | not significant |

(p-values above are unadjusted; both columns are Bonferroni-corrected
across 5 populations before deciding significance.) `cd4_t_cell` only
reaches significance under the CLR-based test. Under raw percentages its
Bonferroni-corrected p-value (≈0.067) just misses the 0.05 threshold,
even though the unadjusted p-value (≈0.013) looks meaningful on its own.
That is the closure problem in action from the other direction: it is not
only capable of manufacturing a false positive (as it would for a
population whose raw-percentage significance does not hold up under CLR),
it can also mask a real one. This is the reason CLR is used throughout
this project (`analysis.py`'s `add_clr_column`), not just noted as a
caveat.

**Caveat, stated plainly:** CLR is not a complete fix. The 5
CLR-transformed values for a sample still sum to exactly zero by
construction, so one linear dependency remains among them (unlike ILR,
which uses 4 orthonormal coordinates and removes the constraint
entirely). CLR was chosen over ILR here because it keeps one
interpretable value per population, matching what the assignment and
Bob actually need ("which populations differ"), whereas ILR's
coordinates are linear combinations across multiple populations at once
and don't map back to a single population cleanly. A fully rigorous
treatment would use ILR or a compositional MANOVA as an omnibus test,
with CLR-based per-population contrasts as follow-up; that's a
legitimately larger analysis than what's built here.

Descriptive statistics (the medians and averages shown throughout the
dashboard) are still reported as percentages, since that's the natural,
interpretable unit for describing composition. Only the significance
test itself uses CLR.

Because 5 populations are tested, a **Bonferroni correction** is applied
(alpha = 0.05 / 5). After correction, `cd4_t_cell` is the only population
that remains statistically significant. `b_cell`, `monocyte`, `nk_cell`,
and `cd8_t_cell` do not. This is the headline finding for Part 3.

**Part 3 asks for melanoma patients specifically.** The assignment's
exact wording is "melanoma patients receiving miraclib," and "please only
include PBMC samples." `get_responder_comparison()` in `analysis.py`
filters on `condition = 'melanoma'`, `treatment = 'miraclib'`, and
`sample_type = 'PBMC'` together, matching that wording exactly: 1,968
samples across 656 subjects (331 responders, 325 non-responders), spanning
2 of the 3 projects (`prj1` and `prj3`). `prj2` has 229 melanoma, miraclib
subjects of its own, but all of them were only ever sampled as WB, never
PBMC, so the PBMC restriction excludes that project entirely from this
particular cohort, not the melanoma restriction on its own. To reproduce
this in the dashboard, apply all three filters (condition, treatment,
sample type) on the Responder vs Non-responder tab, not just treatment
and sample type.

**Paired vs. unpaired tests, and why it matters which one you use.**
Every comparison in this dashboard is Mann-Whitney (unpaired) on CLR
values, *except* the By Population tab, which uses the **Wilcoxon
signed-rank test** (paired) instead. The reason is structural, not a
preference: Responder vs. Non-responder, By Date, and Custom all
compare *different sets of samples* against each other, genuinely
independent groups, which is exactly what Mann-Whitney assumes. By
Population mode compares *the same samples'* `b_cell%` against their
`cd4_t_cell%`, two measurements from the same sample, tied together by
the same compositional closure constraint discussed above. Treating that
as two independent groups (unpaired) would be the wrong test for the same
underlying reason raw percentages needed the CLR fix: it ignores a real
dependency in the data. Wilcoxon signed-rank is the non-parametric,
paired analogue of Mann-Whitney, run on each pair of populations' CLR
values for the matched set of samples present in both.

## Confounder check

Part 3's comparison pools samples across both of the 2 projects present in
this cohort, both sexes, and a range of ages without stratifying by any of
them. If response rates or population baselines differed systematically by
project (batch effects are common in cytometry data pooled across
cohorts/sites), an observed "responders differ" finding could partly be a
project effect rather than a genuine treatment-response signal.

`analysis.py`'s `check_group_balance(df, group_col, stratify_col)`
checks this directly: a chi-square test of independence on a per-subject
contingency table, checking whether `group_col` (e.g. `response`) is
balanced across levels of `stratify_col` (e.g. `project`). It's
general-purpose by construction. It takes whatever dataframe and column
names are passed in, so it recomputes against whatever data is actually
loaded rather than reporting a historical fact about this one CSV. That's
what makes it valid if a different or updated dataset arrives: the check
is code, not documentation. It's surfaced live in the dashboard, in an
expander under the Responder vs Non-responder tab.

For the current dataset, in the Part 3 cohort: **response is balanced
across project** (`prj1` is 58.2% of non-responders and 58.9% of
responders, `prj3` the remainder; p≈0.91) and **across sex** (p≈0.12).
No evidence either is confounding the comparison. A high p-value here
means no evidence of imbalance was found, not that confounding is proven
absent; this is a simple heuristic (p > 0.05), not a formal equivalence
test, and it only covers the two variables checked (project, sex), not
every possible confounder.

## The four comparison tabs

**Responder vs Non-responder** lets you filter on any combination of
variables (response itself isn't offered as a pre-filter, since it's the
comparison axis), and compares responders against non-responders within
that cohort: average cell count/percentage per group, a distribution
boxplot, and a Mann-Whitney test per population.

**By Population** lets you filter on any combination of variables (except
population, for the same reason), then select 2 or more cell populations
to compare directly against each other within that cohort, using the
paired Wilcoxon signed-rank test (see "Statistical approach" above for
why paired, not unpaired).

**By Date** lets you filter on any combination of variables (except
`time_from_treatment_start`), then select 2 or all 3 timepoints to
compare against each other, per population. 3 timepoints tests every
pair (3 tests per population instead of 1).

**Custom** builds 2 to 4 independent cohorts (each with its own full
filter set and a custom label) and compares them all directly against
each other: average cell count/percentage per population, a distribution
boxplot, and a Mann-Whitney test per population, per pair of cohorts.
Cohort labels must be unique. The dashboard checks and shows an error
rather than silently merging two identically-labeled cohorts.

**All four comparison tabs share the same section order** (cohort
summary, average table + charts, frequency table, distribution boxplot,
stats table), so switching between them doesn't mean re-learning the
layout. Confirmed directly by walking the rendered page and checking
each tab's section headers appear in the same sequence, not just
eyeballed.

**Every stats table offers the same toggle**: grouped by cell type by
default (all `b_cell` rows together, then `cd8_t_cell`, in canonical
order), or a checkbox to instead surface Bonferroni-significant results
at the top regardless of population. Both orderings answer a genuinely
different question, "what's going on with this population specifically"
versus "what's the strongest finding here," so this is an explicit
choice, not a default silently picked for you.

**The average-count and average-percentage bar charts facet by
population with independent y-axis scales**, whenever they're showing a
comparison (Responder/Non-responder, By Date, Custom, not the plain
population-colored version in Default/By Population). This isn't
cosmetic: between-population variation in raw cell count (e.g. `b_cell`
~10k vs. `cd4_t_cell` ~30k) is far larger than the real variation across
comparison groups within one population (often under 2%), so a shared
axis makes correctly-computed differences visually invisible. Verified
directly against real data before concluding this was a display problem,
not a computation bug: per-timepoint averages for `b_cell` genuinely
differ (9908 to 9965 to 9909 across the 3 days), just by an amount too
small to see against a ~20,000-unit shared scale.

**These results are exploratory, not confirmatory.** Every comparison
tab lets you re-slice the cohort and re-run the same test machinery in
effectively unlimited ways, a classic garden-of-forking-paths setup.
Bonferroni correction is applied correctly within a single result (across
every population/pair combination actually tested), but there's no
correction across the broader search space of every filter combination
you could try. A caption on each results panel says this explicitly:
p-values from ad hoc slices carry less evidentiary weight than the single
pre-specified Part 3 comparison and shouldn't be read as confirmatory on
their own.

Age can be set either by dragging the slider or by typing exact values
into the Min age / Max age number inputs. The two stay in sync in either
direction (moving the slider updates the number inputs, and vice versa).
Comparison-group order (Responder before Non-responder, chronological
timepoints, cohort A before B before C before D) is pinned explicitly
throughout, both in the underlying data (`get_population_averages`,
`compare_n_groups`) and in each chart's `category_orders`, rather than
left to default alphabetical or click-order sorting. That default
sorting was a real, confirmed inconsistency at one point (the average
table showed "Non-responder" before "Responder" while the distribution
boxplot showed the opposite), and `st.multiselect` returns selections in
the order clicked, not the options list order, which would otherwise let
a timepoint comparison silently render out of chronological order
depending on click order.

Because narrowing a cohort can shrink group sizes quickly, every stats
result reports the current **n per group**, and handles four distinct edge
cases explicitly rather than crashing or silently misleading:

1. **No samples match the selected filters at all**, or (Custom/By Date/By
   Population) **fewer than 2 of the selected groups have any matching
   samples**. One clear message, shown once, before any per-population
   stats are attempted.
2. **Samples exist but have no response data at all**, e.g. filtering to
   `treatment='none'` (untreated/healthy subjects). Reported with an
   explanation of why (response doesn't apply to untreated subjects), not
   just "no data."
3. **Samples exist but a specific population has no individuals in one of
   the groups being compared**. Reported per population, distinct
   from cases 1 and 2.
4. **All groups have data but one or more are small (n < 20)**. The
   result is still computed and shown, with a visible warning that the
   p-value may be unreliable at that sample size. Nothing is hidden.

The n < 20 threshold is a judgment call, not derived from the data, noted
here for transparency.

## Dashboard design

**Theme and color.** The dashboard uses a navy-teal theme
(`.streamlit/config.toml`) plus a consistent color mapping, verified to
stay well-separated (RGB distance) both from each other and across uses:
each of the 5 cell populations always renders in the same soft/muted
color across every chart and table; Responder/Non-responder are always
bluish-green vs. vermillion; Cohort A/Cohort B are always blue vs.
reddish-purple (including a small colored dot next to each cohort's
label). The Response and Cohort colors are drawn from the Okabe-Ito
palette, a peer-reviewed colorblind-safe palette (safe for deuteranopia,
protanopia, and tritanopia). Population names are also color-coded in the
summary tables (average counts, stats results), using pandas' `Styler`.

**Boxplots** fill each box in the compared group's own color, exactly
matching that group's legend swatch. There's never a mismatch between
what the legend shows and what's actually drawn. Population is conveyed
separately: each facet's background is tinted with that population's own
color at a low, fixed opacity (`POPULATION_BACKGROUND_OPACITY`), and the
population name is shown as that facet's x-axis title. A small key below
each boxplot shows exactly which population color maps to which
background tint, built from the *exact* list of populations the
boxplot actually faceted (returned directly by `render_boxplot`), not
independently recomputed from a stats table, since two separate
computations of "which populations, in what order" can drift apart if
either one's filtering logic changes later (this happened once: a
population-filtered cohort produced facets for all 5 populations
regardless of what was actually selected, misaligning every facet-indexed
label onto the wrong facet). The key's swatch opacity matches the actual
facet background exactly, rather than a solid, fully-opaque dot. An
earlier version showed a visibly more saturated color in the key than
what the chart itself displayed.

**Bar charts** (average cell count and average percentage) use
`GROUP_FILL_ALPHA`, the same fill opacity value the boxplot's boxes use,
for both the population-colored and comparison-group-colored variants.
An earlier version left the comparison-group bars fully opaque and the
population-colored bars at the much lower `POPULATION_BACKGROUND_OPACITY`
(meant for background tints, not fills), so the same colors read at
noticeably different intensities depending on chart type. Both now match
the boxplot's fill exactly.

**Cohort summary** (Default, Responder vs Non-responder, By Population,
and By Date tabs, plus Custom's combined per-cohort version) breaks down
subjects by project, condition, treatment, response, and sex: 5 tables
total, always in that order.

**Filters** are collapsed by default inside an expander, with a badge
showing how many filters are currently active even while collapsed (e.g.
"3 filters active"), and a **Reset** button that clears all of that
panel's filters back to the full dataset in one click.

**Exporting.** Every summary table (average cell counts, frequency table,
stats results) has a "Download this table as CSV" button. Charts have
their own built-in PNG export via the camera icon in the Plotly toolbar
that appears on hover.

**Layout.** Related content is grouped into bordered cards
(`st.container(border=True)`) rather than sitting directly on the page,
and dataset-level context (subject/sample/project counts) lives in a
persistent sidebar rather than being repeated inline.

## Caching

The dashboard caches query results via `@st.cache_data`, keyed on
`cell-count.csv`'s modification time. If the CSV changes and `load_data.py`
is rerun, the cache automatically invalidates on the next dashboard
interaction. No manual cache-clearing needed, and no risk of the dashboard
silently serving stale results after the underlying data changes.

## Part 4: baseline cohort summary

Filtered to melanoma, PBMC, miraclib-treated, baseline (`time_from_treatment_start = 0`):

- **656** matching samples
- **384** from `prj1`, **272** from `prj3`
- **331** responder subjects, **325** non-responder subjects
- **344** male subjects, **312** female subjects

(Response/sex counts are per-subject, not per-sample, since each subject has
exactly one response and sex value across all their samples.)

Reproducible via `make pipeline` (writes `part4_baseline_melanoma_samples.csv`
and `part4_summary.txt`) or the same filters in the dashboard's Default tab.

## Known limitations

**Chart height does not adapt to browser zoom.** Every `st.plotly_chart`
call uses `width='stretch'`, so chart width correctly follows the
container as it resizes. Height does not: no chart sets an explicit
height, so all of them fall back to Plotly's fixed default (about 450px)
regardless of the container's actual width. At normal zoom this looks
fine, but at high browser zoom (150-200%) the effective container width
shrinks substantially while height stays fixed, so the faceted charts in
particular (5 population panels sharing that width, in the boxplots and
comparison bar charts) can end up visibly cramped, with x-axis labels
crowding or overlapping. `PLOTLY_CONFIG` in `app.py` does not currently
set Plotly's own `responsive` option, and no chart computes height as a
function of viewport size or facet count. A future update should address
this, most likely by setting `"responsive": True` in `PLOTLY_CONFIG` and
testing whether that alone is sufficient, or, if not, computing each
chart's height dynamically (e.g. from facet count or container width)
rather than relying on Plotly's fixed default.

**Bar width in the average-count/average-percentage charts is not
pinned, only auto-sized.** Neither `_build_avg_bar_chart`'s `px.bar`
calls nor any later `update_layout`/`update_traces` call sets an explicit
bar width or `bargap`/`bargroupgap`. Confirmed directly: building the
same chart with 2 categories versus 5 categories both come back with
`bar.width == None`, meaning Plotly recalculates bar thickness from
scratch on every render, purely as a function of how many categories or
groups happen to be present and how much axis space is available. In
practice this means the same chart can look proportionally different
depending on the current selection: bars read as noticeably thicker when
a population filter narrows 5 populations down to 2, or when Custom
compares 2 cohorts instead of 4, with no consistent bar thickness across
those states. A future update should pin a consistent visual bar width
(e.g. a fixed `bargap`/`bargroupgap` ratio, or an explicit `width` on the
trace scaled to the current category count) so bar thickness stays
visually consistent as the underlying selection changes, rather than
being left to Plotly's default per-render auto-sizing.

**Bar width in the average-count and average-percentage bar charts is
not tuned.** `_build_avg_bar_chart` in `app.py` builds every bar chart
with `px.bar(..., barmode="group")` and never sets `bargap`,
`bargroupgap`, or an explicit per-trace `width`, so bar sizing is left
entirely to Plotly's untuned defaults. The current result does not look
right (bars read as too wide relative to the space between them, and
that changes noticeably between the 5-population single-panel charts on
Default/By Population and the per-population-facet comparison charts on
the other tabs, since facet count and group count both affect how
Plotly's default sizing divides up the available width). This is a
separate issue from the chart-height limitation above; fixing it means
explicitly tuning `bargap`/`bargroupgap` (and possibly setting a fixed
`width` per bar) rather than relying on Plotly's defaults, most likely
with different tuned values for the single-panel and faceted cases,
which is future work not yet done.
