# Loblaw Bio Immune Cell Analysis

Analysis pipeline and interactive dashboard for Bob Loblaw's clinical trial
data, examining how immune cell populations relate to treatment response.

**Live dashboard:** https://loblaw-bio-cell-count-dashboard2-gyjnefkvjgzi5kkkv7wekn.streamlit.app/

Deployed on Streamlit Community Cloud. The app builds `cell_counts.db`
automatically from `cell-count.csv` on first load (see "Self-initializing
database"), so no setup step is needed to view it.

## Quickstart (GitHub Codespaces)

```bash
make setup       # installs dependencies (uses uv if available, falls back to pip)
make pipeline    # builds cell_counts.db, then generates every Part 2-4 output file/plot
make dashboard   # streamlit run app.py, launches the interactive dashboard
make test        # installs dev dependencies, builds cell_counts.db if missing, runs pytest
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
`part4_summary.txt`.

Tested with Python 3.12 and the exact package versions pinned in
`requirements.txt`.




## Part 2: frequency table

Reproducible via `make pipeline` (writes `part2_frequency_table.csv`) or 
the same filters in the dashboard's Default tab.


## Statistical approach (Part 3)

Reproducible via `make pipeline` (writes `part3_boxplot_responders.png` 
and `part3_stats_results.csv`) or the same filters in the dashboard's Default 
tab.


## Part 4: baseline cohort summary

Reproducible via `make pipeline` (writes `part4_baseline_melanoma_samples.csv` 
and `part4_summary.txt`) or the same filters in the dashboard's Default tab.





## Self-initializing database

`app.py` checks for `cell_counts.db` on startup. If it's missing but
`cell-count.csv` is present, it builds the database automatically (calling
the same `load_data.build_database()` function `make pipeline` uses)
before rendering anything. This matters for cloud deployment: platforms
like Streamlit Community Cloud start `app.py` directly and have no
built-in setup step. Locally, `make pipeline` builds the database ahead
of time and this path never fires.

## Tests

```bash
make test
```

runs the full pytest suite, which also runs on every push and pull
request via GitHub Actions (`.github/workflows/tests.yml`).
`tests/test_analysis.py` covers the query/analysis layer directly: known
row counts, the edge-case handling in the stats functions (no samples, no
response data, no individuals in a group, small-n warning), the invariant
that a population filter doesn't change what percentages are computed
against, specific cross-checked numeric results (including the
assignment's form-question answer), the N-group, paired-timepoint, and
paired population-vs-population comparisons, and `load_data.py`'s
validation logic. `tests/test_app.py` uses Streamlit's official AppTest to run the
dashboard headlessly and check representative scenarios per tab. `make
test` is self-sufficient: it installs dev dependencies and builds
`cell_counts.db` first if it's missing.

## Code structure and why

Three layers, each with one job:

- **`load_data.py`** owns getting the raw CSV into a validated, relational
  form. It's the only place that touches `cell-count.csv` directly or knows
  the wide-format-to-long-format conversion. Nothing downstream re-reads the
  CSV.
- **`analysis.py`** owns every query and statistical computation, and knows
  nothing about Streamlit or the dashboard. It's plain functions that take a
  connection or a DataFrame and return a DataFrame, which is what makes the
  statistics unit-testable (`tests/test_analysis.py`) without spinning up
  the UI.
- **`app.py`** owns presentation only: it calls `analysis.py` functions and
  renders the results. It's deliberately general-purpose (any cohort, any
  comparison, via filters) rather than having one fixed view per assignment
  part, because Parts 2-4 are special cases of "filter a cohort, then
  optionally compare groups within it."

That separation is also why there's a **`generate_outputs.py`**, distinct
from the dashboard: the assignment's outputs need to exist as concrete
files a grader can open without clicking through the UI, so this script
calls the same `analysis.py` functions directly and writes them to disk.
`app.py` and `generate_outputs.py` are two callers of the same underlying
logic, not two separately maintained copies of it.

## Repository contents

| File | Purpose |
|---|---|
| `Makefile` | `make setup` / `make pipeline` / `make dashboard` / `make test` |
| `.streamlit/config.toml` | Dashboard theme (colors, font) |
| `.github/workflows/tests.yml` | Runs the test suite on every push/PR |
| `schema.sql` | Relational schema (3 tables) |
| `load_data.py` | Validates `cell-count.csv`, builds `cell_counts.db` |
| `analysis.py` | Query + analysis functions (framework agnostic, no Streamlit) |
| `generate_outputs.py` | Produces the required Part 2-4 files/plot by calling `analysis.py` directly |
| `app.py` | Streamlit dashboard: 5 tabs, single-cohort exploration plus 4 comparison modes |
| `tests/` | pytest suite (`test_analysis.py`, `test_app.py`, `conftest.py`) |
| `cell-count.csv` | Source data |
| `cell_counts.db` | Generated by `make pipeline`; gitignored as an intermediate build artifact, quickly rebuilt from the CSV |
| `part2_frequency_table.csv` | Part 2's required output, generated by `make pipeline` and committed |
| `part3_stats_results.csv`, `part3_boxplot_responders.png` | Part 3's required outputs, generated and committed |
| `part4_baseline_melanoma_samples.csv`, `part4_summary.txt` | Part 4's required outputs, generated and committed |
| `requirements.txt` | Pinned runtime dependencies |
| `requirements-dev.txt` | Runtime dependencies plus `pytest` |

## Data dictionary (`cell-count.csv`)

| Column | Type | Meaning | Valid values |
|---|---|---|---|
| `project` | text | Study/project identifier | `prj1`, `prj2`, `prj3` |
| `subject` | text | Unique subject/patient identifier | e.g. `sbj000` |
| `condition` | text | Subject's diagnosis | `melanoma`, `carcinoma`, `healthy` |
| `age` | integer | Subject's age in years | 50-79 in this dataset |
| `sex` | text | Subject's sex | `M`, `F` |
| `treatment` | text | Treatment arm | `miraclib`, `phauximab`, `none` (untreated, including all healthy subjects) |
| `response` | text | Treatment response | `yes`, `no`, or null (null exactly when `treatment='none'`) |
| `sample` | text | Unique sample identifier | e.g. `sample00000` |
| `sample_type` | text | Tissue/specimen type | `PBMC`, `WB` (whole blood) |
| `time_from_treatment_start` | integer | Days since treatment start at sample draw | `0`, `7`, `14` |
| `b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, `monocyte` | integer | Raw cell count for that population in this sample | non-negative integers |

Each subject has exactly one `project` / `condition` / `age` / `sex` /
`treatment` / `response` value across all their samples, and exactly 3
samples (one per timepoint).

## Schema design

```
subjects (subject_id PK, project, condition, age, sex, treatment, response)
samples  (sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
cell_counts (sample_id FK, population, count, PRIMARY KEY(sample_id, population))
```

**Rationale:**

- **`subjects`** holds attributes that are fixed per person. Every subject
  has exactly one `project` / `condition` / `age` / `sex` / `treatment` /
  `response` value across all their samples, so storing these once (rather
  than repeating them on every sample row) avoids redundancy and keeps
  updates consistent.
- **`samples`** is one row per physical sample draw, the unit that varies by
  timepoint (`time_from_treatment_start`) and `sample_type` (PBMC/WB).
- **`cell_counts`** is stored in long format (one row per population per
  sample) rather than as five wide columns. This makes "group by population"
  queries (Part 2's frequency table, Part 3's per-population stats) a plain
  SQL/pandas group-by instead of manually unpivoting five columns in code,
  and it's the standard normalized shape for repeated-measurement data.
- **`response` is nullable by design**, not an error state. It is `NULL`
  exactly when `treatment = 'none'`: untreated (including healthy) subjects
  have no response outcome to record. `load_data.py` validates this
  relationship in both directions before writing any rows.
- **Foreign key constraints are enforced** (`PRAGMA foreign_keys = ON`), so
  the relational integrity in the schema is actually checked at insert time,
  not just declared as documentation.

**Scalability.** This schema and SQLite hold up fine at the current size
(3,500 subjects, 10,500 samples, 52,500 cell-count rows) and would still
be comfortable at 10-100x that. Past that, three things would need to
change:

- **`project` becomes its own table.** At hundreds of projects,
  project-level metadata (PI, institution, funding source, start date)
  needs somewhere to live that isn't duplicated across every subject row:
  a `projects` table with `subjects.project_id` as a foreign key, the same
  normalization logic already applied to `subjects` and `samples`.
- **Indexes on the columns actually filtered on.** `subject_id`,
  `sample_type`, `time_from_treatment_start`, and `population` are the
  columns every query in this project group-bys or filters on. SQLite
  handles the current volume without extra indexes, but at real scale
  these need explicit indexes, or every filtered query degrades to a full
  table scan.
- **SQLite itself becomes the limiting factor before the schema does.**
  SQLite is single-writer and file-based, which is exactly right for a
  self-contained, reproducible analysis pipeline, but the wrong choice
  once multiple analysts need concurrent access. The schema is standard
  relational SQL and moves to Postgres without redesign: a deployment
  decision, not a schema one.

**For "various types of analytics"**, the shape that generalizes is
already in place: `cell_counts` being long-format, keyed by `sample_id` +
`population`, means adding a 6th cell population is just new rows, not a
schema change. The same pattern extends beyond cell counts: a new data
type per sample (genomic variants, imaging features, longitudinal labs)
is a new fact table keyed by `sample_id`, sitting alongside `cell_counts`
rather than inside it, a star-schema-style extension rather than a
rewrite of `subjects` or `samples`.

## Data validation

`load_data.py` validates `cell-count.csv` before writing anything to the
database. Every check is independent (a missing column doesn't suppress
other findings), and all failures found are listed in one error message,
not just the first one encountered.

**Hard failures** (abort the load, nothing is written):
- Missing required columns
- Negative or non-numeric cell counts
- Nulls in required (non-nullable) fields
- Duplicate sample IDs
- `response` inconsistent with `treatment` (non-null when `treatment='none'`,
  or null when it isn't)

**Soft warnings** (printed, load continues):
- Categorical values outside the expected set, flagged rather than
  rejected since they could reflect a legitimate new category
- Implausible age values (<0 or >120)

These checks are structural and re-run on whatever CSV is present, so a
non-conforming replacement dataset fails loudly without code changes. Two
caveats worth stating explicitly: a genuinely new cell population (a 6th
count column) is silently ignored rather than flagged, because the
wide-to-long conversion melts a fixed `POPULATIONS` list (supporting a
6th population means adding its name to that list in both `load_data.py`
and `analysis.py`); and tests that pin exact values (row counts, cohort
sizes) are deliberately tied to the current dataset as regression checks,
so a genuinely different CSV would make them fail correctly until their
expected values are updated.

## Part 2: frequency table

Reminder: reproducible via `make pipeline` (writes `part2_frequency_table.csv`) or 
the same filters in the dashboard's Default tab.

`percentage` and "relative frequency" are the same quantity, per the
assignment's own column definition ("relative frequency in percentage").
`get_frequency_table()` computes it exactly as specified:

```
total_count = sum of count across all 5 populations, for that sample
percentage  = count / total_count * 100
```

The 5 `percentage` values for any one sample sum to 100. `total_count` is
a per-sample quantity, repeated once per population row since the table
is one row per sample-population pair. This is a plain descriptive
statistic; no transform (like the CLR used for Part 3's significance
test) is applied here.

## Statistical approach (Part 3)

Reminder: reproducible via `make pipeline` (writes `part3_boxplot_responders.png` 
and `part3_stats_results.csv`) or the same filters in the dashboard's Default 
tab.

Responders vs. non-responders are compared per cell population using the
**Mann-Whitney U test** rather than a t-test. Cell-frequency data is
bounded and often skewed, so a non-parametric test that doesn't assume
normality is the safer default.

**The test runs on a CLR-transformed value, not the raw percentage.** The
5 cell populations are compositional data: they sum to 100% of each
sample's total cell count, so an increase in one population mechanically
forces the others down. That "closure" constraint means the 5 populations
aren't independent, and testing directly on raw percentages treats them
as if they were, which can manufacture or mask apparent significance. The
fix is the **centered log-ratio (CLR) transform**: for each sample,
`clr_i = ln(count_i) - mean(ln(counts))` across that sample's own 5
populations, moving the data into ordinary, unconstrained real space
before testing.

This changes the answer for this dataset. Tested both ways on the same
melanoma, miraclib, PBMC cohort Part 3 asks for:

| Population | Raw-percentage test | CLR-based test |
|---|---|---|
| `cd4_t_cell` | not significant (p≈0.0133) | **significant** (p≈0.0025) |
| `b_cell` | not significant (p≈0.0557) | not significant (p≈0.1364) |
| `nk_cell` | not significant | not significant |
| `monocyte` | not significant | not significant |
| `cd8_t_cell` | not significant | not significant |

(p-values above are unadjusted; both columns are Bonferroni-corrected
across 5 populations before deciding significance.) `cd4_t_cell` only
reaches significance under the CLR-based test: under raw percentages its
Bonferroni-corrected p-value (≈0.067) just misses the 0.05 threshold. The
closure constraint can mask a real signal, not just manufacture a false
one, which is why CLR is used throughout this project
(`analysis.py`'s `add_clr_column`) rather than just noted as a caveat.

**Caveat:** CLR is not a complete fix. The 5 CLR-transformed values for a
sample still sum to exactly zero by construction, so one linear
dependency remains among them (unlike ILR, which uses 4 orthonormal
coordinates and removes the constraint entirely). CLR was chosen over ILR
because it keeps one interpretable value per population, matching what
the assignment and Bob actually need ("which populations differ"),
whereas ILR's coordinates are linear combinations across multiple
populations and don't map back to a single population cleanly.
Descriptive statistics (medians and averages shown in the dashboard) are
still reported as percentages, the natural unit for composition; only the
significance test itself uses CLR.

Because 5 populations are tested, a **Bonferroni correction** is applied
(alpha = 0.05 / 5). After correction, `cd4_t_cell` is the only population
that remains statistically significant. This is the headline finding for
Part 3.

**Part 3 asks for melanoma patients specifically.** The assignment's
wording is "melanoma patients receiving miraclib," and "please only
include PBMC samples." `get_responder_comparison()` filters on
`condition = 'melanoma'`, `treatment = 'miraclib'`, and
`sample_type = 'PBMC'` together: 1,968 samples across 656 subjects (331
responders, 325 non-responders), spanning `prj1` and `prj3`. `prj2` has
229 melanoma miraclib subjects of its own, but all of them were only ever
sampled as WB, never PBMC, so the PBMC restriction excludes that project
from this cohort. To reproduce this in the dashboard, apply all three
filters (condition, treatment, sample type) on the Responder vs
Non-responder tab.

**Paired vs. unpaired tests.** Which test a tab uses follows from
whether its comparison groups are independent:

- **Responder vs Non-responder** and **Custom** compare different sets
  of subjects and use the unpaired Mann-Whitney U test on CLR values.
- **By Population** compares measurements taken from the same samples
  (one sample's `b_cell%` against its own `cd4_t_cell%`, tied together
  by the closure constraint), so it uses the paired **Wilcoxon
  signed-rank** test on the matched samples.
- **By Date** also uses the paired **Wilcoxon signed-rank**, pairing on
  subject (`compare_n_groups_paired` in `analysis.py`). Every subject
  in this dataset is sampled exactly once at each timepoint with a
  constant sample type, so any cohort's "Day 0" and "Day 14" groups are
  the same subjects measured repeatedly -- repeated measures, not
  independent groups. The test therefore runs on each subject's own
  within-subject change rather than pooling the two days as if they
  were unrelated crowds.

Like the CLR choice above, the paired choice changes the answer for
this dataset. Among melanoma miraclib PBMC **responders**, `cd4_t_cell`
rises during treatment (a median within-subject change of roughly +1.2
percentage points from day 0 to day 14, while non-responders stay
flat). Tested unpaired, that trend does not survive Bonferroni
correction across the tab's 15 tests (p=0.0103, corrected ~0.15);
tested paired, it does (p=0.0032, corrected ~0.048). A regression test
pins this exact divergence
(`test_paired_by_date_finds_the_responder_cd4_trend_unpaired_misses`).
Note: an unpaired test on paired data costs sensitivity. It can miss 
real within-subject changes but does not manufacture false positives
which is also why this distinction never touched Part 3, where responders
and non-responders are different people and unpaired Mann-Whitney is 
the right tool.

## Confounder check

Part 3's comparison pools samples across 2 projects, both sexes, and a
range of ages without stratifying by any of them. If response rates or
population baselines differed systematically by project (batch effects
are common in cytometry data pooled across sites), an observed
"responders differ" finding could partly be a project effect rather than
a genuine treatment-response signal.

`analysis.py`'s `check_group_balance(df, group_col, stratify_col)` checks
this directly: a chi-square test of independence on a per-subject
contingency table. It's general-purpose by construction, recomputing
against whatever data is loaded rather than reporting a fixed historical
finding, and it's surfaced live in an expander under the Responder vs
Non-responder tab.

For the current dataset, in the Part 3 cohort, response is balanced
across project (p≈0.91) and across sex (p≈0.12): no evidence either is
confounding the comparison. A high p-value means no evidence of imbalance
was found, not that confounding is proven absent; this is a simple
heuristic (p > 0.05), not a formal equivalence test, and it only covers
the two variables checked.

## Dashboard: one page, 5 tabs, each independently filterable

The dashboard uses `st.tabs`, which renders every tab's content on every
rerun and just CSS-hides the inactive ones, so each tab's filter
selections survive switching tabs (an `if/elif` view switcher would drop
widget state for whichever views aren't currently rendered).

The 5 tabs are **Default**, **Responder vs Non-responder**, **By
Population**, **By Date**, and **Custom**. Default is single-cohort
exploration (cohort summary, average cell counts, the frequency table,
and a per-population distribution boxplot). The other four each compare 2
or more groups, split a different way:

- **Responder vs Non-responder**: always exactly 2 groups, the
  assignment's Part 3 axis. Response isn't offered as a pre-filter here,
  since it's the comparison axis itself.
- **By Population**: select 2+ cell populations and compare them directly
  against each other within the same cohort, using the paired Wilcoxon
  signed-rank test (see "Statistical approach" for why paired).
- **By Date**: select 2 or all 3 timepoints and compare them against each
  other, per population, using the paired Wilcoxon signed-rank test:
  every subject is sampled at every timepoint, so timepoint groups are
  the same subjects measured repeatedly (see "Statistical approach").
  Selecting 3 timepoints tests every pair.
- **Custom**: build 2 to 4 independent cohorts with any combination of
  filters and compare them all against each other. Cohort labels must be
  unique; the dashboard shows an error rather than silently merging
  identically-labeled cohorts.

Any number of groups beyond 2 is handled by testing every pairwise
combination, Bonferroni-corrected across all pairs and populations
tested together (`compare_n_groups` in `analysis.py` for independent
groups, `compare_n_groups_paired` for By Date's repeated measures). The
correction gets stricter fast as groups are added (3 groups is 3 pairs
per population, 4 groups is 6), which is the real statistical cost of
pairwise tests over a single omnibus test like Kruskal-Wallis.

All four comparison tabs share the same section order (cohort summary,
average table + charts, frequency table, distribution boxplot, stats
table). Every stats table offers the same toggle: grouped by cell type in
canonical order by default, or a checkbox to surface
Bonferroni-significant results at the top. Comparison-group order
(Responder before Non-responder, chronological timepoints, cohort A
before B) is pinned explicitly throughout rather than left to
alphabetical or click-order sorting.

The comparison bar charts facet by population on a shared y-axis, so
between-population scale differences read directly from bar heights;
differences between groups *within* a population are usually small and
are easier to read from the averages table, the stats table, and the
distribution boxplot.

**These results are exploratory, not confirmatory.** The comparison tabs
let you re-slice the cohort and re-run the same test machinery in
effectively unlimited ways, a classic garden-of-forking-paths setup.
Bonferroni correction is applied within a single result, but there's no
correction across the broader search space of every filter combination
you could try. A caption on each results panel says this explicitly:
p-values from ad hoc slices carry less evidentiary weight than the single
pre-specified Part 3 comparison.

Because narrowing a cohort can shrink group sizes quickly, every stats
result reports the current n per group and handles four edge cases
explicitly rather than crashing or silently misleading:

1. No samples match the selected filters at all (or fewer than 2 selected
   groups have any samples): one clear message, shown once.
2. Samples exist but have no response data at all, e.g. filtering to
   `treatment='none'`: reported with an explanation of why.
3. A specific population has no individuals in one of the groups being
   compared: reported per population.
4. All groups have data but one or more are small (n < 20): the result is
   still computed and shown, with a visible warning that the p-value may
   be unreliable. The n < 20 threshold is a judgment call, noted here for
   transparency.

## Dashboard design

The dashboard uses a navy-teal theme (`.streamlit/config.toml`) with a
consistent color mapping: each of the 5 cell populations always renders
in the same color across every chart and table, and
Responder/Non-responder and Cohort A/B pairs are drawn from the
Okabe-Ito colorblind-safe palette. Boxplot boxes fill in the compared
group's own color, exactly matching the legend; population is conveyed
separately by tinting each facet's background with that population's
color at low opacity, with a small key below each boxplot built from the
exact list of populations the figure actually faceted.

Filters are collapsed inside an expander with a badge showing how many
are active, plus a Reset button. Age can be set by slider or by typing
exact values into synced Min/Max number inputs. Every summary table has a
"Download this table as CSV" button, and charts have Plotly's built-in
PNG export. Related content is grouped into bordered cards, and
dataset-level context (subject/sample/project counts) lives in a
persistent sidebar.

## Caching

The dashboard caches query results via `@st.cache_data`, keyed on
`cell-count.csv`'s modification time. If the CSV changes and
`load_data.py` is rerun, the cache automatically invalidates on the next
interaction, so the dashboard never silently serves stale results after
the underlying data changes.

## Part 4: baseline cohort summary

Reminder: Reproducible via `make pipeline` (writes `part4_baseline_melanoma_samples.csv` 
and `part4_summary.txt`) or the same filters in the dashboard's Default tab.

Filtered to melanoma, PBMC, miraclib-treated, baseline
(`time_from_treatment_start = 0`):

- **656** matching samples
- **384** from `prj1`, **272** from `prj3`
- **331** responder subjects, **325** non-responder subjects
- **344** male subjects, **312** female subjects

(Response/sex counts are per-subject, not per-sample, since each subject
has exactly one response and sex value.)



## Known limitations

**Chart height does not adapt to browser zoom.** Charts stretch to the
container's width, but no chart sets an explicit height, so all fall back
to Plotly's fixed default (about 450px). At high browser zoom (150-200%)
the faceted charts in particular can end up cramped, with x-axis labels
crowding. A future update should set Plotly's `responsive` option and, if
needed, compute chart height from facet count.

**Bar width in the average-count/percentage charts is not tuned.**
`_build_avg_bar_chart` never sets `bargap`, `bargroupgap`, or an explicit
bar width, so bar thickness is left to Plotly's per-render auto-sizing
and varies with how many categories or groups are currently selected. A
future update should pin a consistent ratio, likely with different tuned
values for the single-panel and faceted cases.
