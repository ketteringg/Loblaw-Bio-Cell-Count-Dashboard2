"""
app.py

Interactive dashboard for the Loblaw Bio cell-count analysis.

Five tabs: Default (single-cohort exploration -- cohort summary, average
cell counts, frequency table, population distribution), and four
comparison tabs each splitting the data a different way -- Responder vs
Non-responder, By Population, By Date, and Custom (2-4 independently
filtered cohorts). Every comparison tab shares the same section order
(cohort summary, average table + charts, frequency table, distribution
boxplot, stats table) and the same stats-table controls (grouped by cell
type by default, with a toggle to surface significant results instead).
There are no fixed tabs for the assignment's Part 2/3/4 cohorts
specifically -- those are special cases of what these tabs already
produce by setting the equivalent filters, so a dedicated tab for each
would just be showing the same underlying data twice. See README.md for
how to reproduce those specific answers.

Run with:
    streamlit run app.py
"""
from pathlib import Path

import pandas as pd
import plotly.express as px
import sqlite3
import streamlit as st

from analysis import (
    POPULATIONS,
    SMALL_N_THRESHOLD,
    format_pvalue,
    check_group_balance,
    run_stats_test_safe,
    get_full_dataset,
    filter_dataset,
    get_cohort_summary,
    get_filtered_frequency_table,
    get_population_averages,
    compare_n_groups,
    compare_populations_paired,
)

ROOT = Path(__file__).parent
DB_PATH = ROOT / "cell_counts.db"
CSV_PATH = ROOT / "cell-count.csv"

st.set_page_config(page_title="Loblaw Bio", layout="wide")

# ---------- visual constants ----------
# Cohort (A/B) colors use the Okabe-Ito palette, a peer-reviewed
# colorblind-safe palette (safe for deuteranopia, protanopia, and
# tritanopia); Responder/Non-responder deliberately reuse the exact same
# 2 colors (see below). Population colors are a separate, soft/pastel
# palette -- verified programmatically (RGB Euclidean distance) to stay
# well-separated both from these 2 "meaning" colors and from each other,
# so no color is ever ambiguous between "which population" and "which
# group" at a glance, even though they're deliberately gentler than the
# meaning colors.
POP_COLORS = {
    "b_cell": "#7DB760",       # soft green
    "cd8_t_cell": "#6067B7",   # soft blue-violet
    "cd4_t_cell": "#BCB96B",   # soft gold
    "nk_cell": "#60B0B7",      # soft teal
    "monocyte": "#B85951",     # soft clay
}
COHORT_COLOR_SEQUENCE = ["#0072B2", "#CC79A7"]  # blue / reddish purple (Okabe-Ito)
# Responder/Non-responder deliberately reuse Cohort A/B's exact colors --
# both are "first group vs second group in a two-way split," so using the
# same color pairing throughout keeps that visual language consistent
# across every comparison mode rather than introducing a second pair.
RESPONSE_COLORS = {"Responder": COHORT_COLOR_SEQUENCE[0], "Non-responder": COHORT_COLOR_SEQUENCE[1]}

# General-purpose color sequence for comparisons with an arbitrary number
# of groups (By Date, Custom with 3+ cohorts): the full Okabe-Ito
# colorblind-safe palette (minus black, which we deliberately avoid using
# for a data color -- see population/gridline color choices elsewhere),
# in an order chosen so the first 2 entries stay visually distinct from
# both RESPONSE_COLORS and COHORT_COLOR_SEQUENCE above (not that they'd
# ever appear in the same chart, but for consistency).
N_GROUP_COLOR_SEQUENCE = [
    "#0072B2",  # blue
    "#CC79A7",  # reddish purple
    "#009E73",  # bluish green
    "#D55E00",  # vermillion
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
]

# Boxplot fill opacity (fill only -- outline/median stay fully opaque, see
# render_boxplot). Moderate rather than near-solid: 0.85 read as too heavy/
# opaque against the population background tint underneath -- this value
# still reads clearly as "filled," just without dominating the chart.
GROUP_FILL_ALPHA = 0.55

# Population background tint opacity (facet backgrounds in render_boxplot,
# population-colored table cells, population-colored bar fills, and the
# population key swatches -- all deliberately share this one constant so
# they stay visually consistent with each other).
#
# Moderate rather than very light: 0.15 read as too faint, especially as a
# bar-chart fill (bars need to be clearly visible as the primary data
# there, not just a background wash). Kept well short of full opacity,
# though, because of documented history: at full opacity, the strongly
# saturated, differently-hued backgrounds behind each boxplot facet
# created a real optical effect where identical gridlines appeared to
# have different thickness depending on which background they crossed
# (confirmed the gridlines themselves were byte-for-byte identical across
# every facet -- same width, color, range, tick spacing -- so the
# apparent inconsistency was purely the background contrast). That issue
# was specifically observed near full opacity, not at this more moderate
# level, but the margin is deliberate rather than assumed safe.
POPULATION_BACKGROUND_OPACITY = 0.3

CUSTOM_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'Inter', sans-serif;
}

h1, h2, h3 {
    letter-spacing: -0.01em;
}

section[data-testid="stSidebar"] {
    border-right: 1px solid #E3E8EA;
}

.accent-bar {
    height: 4px;
    width: 100%;
    background: linear-gradient(90deg, #1F6F78 0%, #4FA8A0 40%, #E0A458 100%);
    border-radius: 2px;
    margin: 4px 0 18px 0;
}

.filter-badge {
    display: inline-block;
    padding: 2px 10px;
    border-radius: 12px;
    font-size: 12px;
    font-weight: 500;
}
.filter-badge.active {
    background: #E4F2F1;
    color: #1F6F78;
}
.filter-badge.inactive {
    background: #F1F2F4;
    color: #6B7280;
}

.cohort-dot {
    display: inline-block;
    width: 10px;
    height: 10px;
    border-radius: 50%;
    margin-right: 6px;
}
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)


def hex_to_rgba(hex_color: str, alpha: float) -> str:
    """Converts a hex color to an rgba() string with the given opacity."""
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r}, {g}, {b}, {alpha})"


def style_population_column(df: pd.DataFrame, column: str = "population") -> "pd.io.formats.style.Styler":
    """Colors the given column's cells using POP_COLORS, so population
    names read consistently with the same colors used in every chart.
    Background opacity matches POPULATION_BACKGROUND_OPACITY, the same
    value used for the boxplot's facet backgrounds and the average-count
    bar charts, so population color opacity is consistent everywhere it
    appears. Text is dark, not white: the population palette is
    intentionally soft/light, and white text on these backgrounds falls
    well short of readable contrast (verified: ~2:1 at best, vs. a 4.5:1
    minimum for body text). Dark text gives 9:1+ on every population
    color, and stays well above that even at reduced background opacity
    (opacity can only lighten the effective background further)."""
    def _colorize(col):
        if col.name != column:
            return ["" for _ in col]
        return [
            f"background-color: {hex_to_rgba(POP_COLORS.get(v, '#F1F2F4'), POPULATION_BACKGROUND_OPACITY)}; "
            f"color: #1A2733; font-weight: 600;"
            for v in col
        ]
    return df.style.apply(_colorize)


def download_csv_button(df: pd.DataFrame, filename: str, label: str, key: str):
    """Small wrapper around st.download_button for a dataframe as CSV."""
    st.download_button(
        label, data=df.to_csv(index=False).encode("utf-8"),
        file_name=filename, mime="text/csv", key=key,
    )


PLOTLY_CONFIG = {
    "displaylogo": False,
    "toImageButtonOptions": {"format": "png", "scale": 2},
}


# ---------- connection + caching ----------
# Cache key includes the CSV's mtime, so cached results self-invalidate the
# moment cell-count.csv changes and load_data.py is rerun. No manual
# cache-clearing needed, and no risk of silently serving stale results.
def _csv_mtime() -> float:
    return CSV_PATH.stat().st_mtime if CSV_PATH.exists() else 0.0


@st.cache_resource
def get_connection(_mtime: float) -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


@st.cache_data
def cached_full_dataset(_mtime: float) -> pd.DataFrame:
    return get_full_dataset(get_connection(_mtime))


# ---------- shared rendering helpers ----------

def darken_hex(hex_color: str, factor: float = 0.35) -> str:
    """Returns a darker shade of hex_color (multiplies each RGB channel
    toward 0 by `factor`). Used for box outlines/median lines: Plotly's
    Box trace has exactly one `line` property that controls BOTH the
    outline and the median line (verified directly -- there's no separate
    median styling), so setting it identical to fillcolor makes the
    median line invisible against the fill. A darker shade keeps both
    visible without introducing a third color."""
    def to_rgb(h):
        h = h.lstrip("#")
        return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))
    r, g, b = to_rgb(hex_color)
    return f"#{round(r * (1 - factor)):02X}{round(g * (1 - factor)):02X}{round(b * (1 - factor)):02X}"


def render_boxplot(
    comparison_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    x_col: str = "response_label",
    group_order: list[str] | None = None,
    group_colors: list[str] | None = None,
):
    """Boxplot faceted by population, supporting 2 or more groups (e.g.
    Responder/Non-responder, or Cohort A/B/C, or Day 0/Day 7/Day 14). Box
    fill is each group's own signature color -- exactly matching the
    legend swatch, so there's no mismatch between what the legend shows
    and what's actually drawn. The population name
    is shown as each facet's x-axis title at the very bottom; the group
    names (Cohort A / Cohort B etc.) are shown just above that, in each
    group's own color -- replacing Plotly's native x-axis tick labels,
    which don't support per-tick coloring. The top facet-title area is
    reserved for the p-value only. Box outlines (and median lines --
    Plotly's Box trace shares one `line` property for both) use a darker,
    fully opaque shade of the same group color; the fill itself is drawn
    with slightly less than full opacity (an rgba() fillcolor, while
    line/marker stay solid hex) so only the fill -- not the outline or
    median line -- is affected. Each facet's background is tinted with
    that population's own color (same hex value used everywhere else),
    at a moderate opacity (POPULATION_BACKGROUND_OPACITY), using
    layer="below" -- which Plotly specifically defines as below gridlines
    (not just below traces), so gridlines and the box traces themselves
    both still render on top regardless of the background's opacity.

    Returns (fig, present_populations): present_populations is the exact,
    already-computed list of populations this figure actually faceted,
    in the exact order they were graphed. Callers should pass this
    directly to render_population_key rather than recomputing their own
    guess at the same list -- two independent computations of "which
    populations, in what order" can drift apart if either one's filtering
    logic changes later, even though they happen to agree today.
    """
    groups = group_order or sorted(comparison_df[x_col].dropna().unique())
    default_colors = ["#6B7280", "#9CA3AF", "#D1D5DB", "#4B5563", "#111827", "#78716C"]
    colors = group_colors or default_colors[:len(groups)]

    # Computed BEFORE building the figure, and used directly in
    # category_orders below -- this must be the actual list of facets
    # Plotly creates. Passing the full POPULATIONS constant instead (all
    # 5, regardless of what's actually in comparison_df) was a real bug:
    # Plotly Express creates one facet per category_orders entry even
    # when a category has zero matching rows, so a population-filtered
    # cohort (e.g. just b_cell + cd4_t_cell selected) produced 5 facets,
    # 3 of them empty with Plotly's raw unprocessed "population=X" title
    # text (never reached by the p-value-replacement loop below, which
    # only iterates present populations). Worse, the later loops that
    # position backgrounds/labels/titles by index (present_populations[i]
    # -> facet slot i) silently landed on the WRONG facet whenever an
    # empty facet was interspersed before a present one in canonical
    # order -- e.g. cd4_t_cell's own label rendering on cd8_t_cell's
    # empty facet, confirmed directly against a real filtered cohort.
    present_populations = [p for p in POPULATIONS if p in comparison_df["population"].unique()]

    fig = px.box(
        comparison_df, x=x_col, y="percentage", facet_col="population",
        facet_col_wrap=5, points="outliers",
        category_orders={"population": present_populations, x_col: groups},
        labels={"percentage": "% of total cells", x_col: ""},
        color=x_col,  # forces one trace per (facet, group); overridden below
    )

    # Traces are ordered facet-major, len(groups) per facet, in the order
    # given by category_orders -- verified empirically against this
    # Plotly version. Fill uses an rgba() string with GROUP_FILL_ALPHA so
    # only the fill is translucent; line/marker stay fully opaque solid
    # hex, since Box traces have exactly one blanket `opacity` that would
    # otherwise fade everything (fill, outline, median line) together.
    for i, trace in enumerate(fig.data):
        group_idx = groups.index(trace.name)
        base_color = colors[group_idx]
        outline = darken_hex(base_color)
        trace.fillcolor = hex_to_rgba(base_color, GROUP_FILL_ALPHA)
        trace.line.color = outline
        trace.marker.color = outline
        trace.showlegend = False

    # Legend: one swatch per group, in the exact (fully opaque) colors
    # used for the outline/marker -- no population swatches here, since
    # population is now conveyed by the x-axis titles, not by any color.
    for idx, group in enumerate(groups):
        fig.add_scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=10, color=colors[idx]),
            name=group, showlegend=True,
        )

    # Share one y-axis scale across every population facet, so box
    # heights are directly comparable at a glance. Gridlines still render
    # on every facet (for visual alignment), but only the leftmost facet
    # shows the numeric tick labels -- repeating identical labels on
    # every facet once the scale is shared is redundant.
    #
    # The range is set explicitly (not left to autorange) with headroom
    # above the actual data max, so the group-name text labels added
    # below have dedicated empty space to sit in rather than overlapping
    # the topmost whiskers/outlier points.
    data_max = comparison_df["percentage"].max()
    data_min = comparison_df["percentage"].min()
    top_padding = (data_max - data_min) * 0.18
    y_range = [max(0, data_min - (data_max - data_min) * 0.05), data_max + top_padding]

    # Gridline color is set explicitly (not left at Plotly's default pale
    # gray) because the default has poor contrast against the population
    # background tints -- verified via WCAG contrast ratio against all 5
    # population colors at full opacity (the worst case): the default
    # gray scores as low as 1.76:1 (well under the 3:1 minimum for visible
    # graphical elements), while pure black is the only candidate tested
    # that clears 3:1 against every population color (worst case 4.12:1,
    # on cd8_t_cell's blue-violet background). Now that the background
    # opacity is well under full opacity (POPULATION_BACKGROUND_OPACITY =
    # 0.3), contrast is even better than that worst-case check, so black
    # remains a safe choice. zeroline is explicitly disabled: the y-axis range's lower
    # bound lands exactly on 0 (see y_range above), so Plotly's separate
    # zero-line was stacking directly on top of the regular gridline at
    # that same position, rendering as a visibly thicker line there than
    # at any other tick -- confirmed by checking the computed range, not
    # just guessed. Tick labels are shown on every facet (not just the
    # first), even though the scale is shared, per explicit request.
    fig.update_yaxes(
        matches="y", range=y_range,
        gridcolor="#000000", gridwidth=0.5,
        zeroline=False, showticklabels=True,
    )

    # Native x-axis tick labels (which would show "Cohort A" / "Cohort B"
    # etc. under each box, all in one uniform color) are hidden -- Plotly
    # doesn't support coloring individual tick labels differently, so the
    # group names are drawn as custom annotations instead, colored to each
    # group's own color, positioned just below the plot area.
    fig.update_xaxes(showticklabels=False)

    fig.update_layout(legend_title_text="", margin=dict(b=90))

    # Background tint per facet, in that facet's population color -- same
    # hue/hex value used everywhere else, at reduced opacity (see
    # POPULATION_BACKGROUND_OPACITY above for why).
    for i, pop in enumerate(present_populations):
        xaxis_suffix = "" if i == 0 else str(i + 1)
        fig.add_shape(
            type="rect", xref=f"x{xaxis_suffix} domain", yref=f"y{xaxis_suffix} domain",
            x0=0, x1=1, y0=0, y1=1,
            fillcolor=POP_COLORS[pop], opacity=POPULATION_BACKGROUND_OPACITY, line_width=0,
            layer="below",
        )

    # Group-name labels, color-coded to each group's own color, positioned
    # just below the plot area (replacing the hidden native tick labels).
    # The population name (x-axis title, added below) uses `standoff` to
    # sit further down still, so the two never collide.
    for i, pop in enumerate(present_populations):
        xaxis_suffix = "" if i == 0 else str(i + 1)
        for idx, group in enumerate(groups):
            fig.add_annotation(
                x=group, y=0,
                xref=f"x{xaxis_suffix}", yref=f"y{xaxis_suffix} domain",
                text=group, showarrow=False,
                font=dict(size=10, color=colors[idx]),
                yanchor="top", yshift=-6,
            )

    # Population name moves to each facet's x-axis title, below the
    # color-coded group labels added above -- not a color, and not at the
    # top. The top facet-title area is reserved for the p-value only.
    # standoff pushes the title further from the axis than Plotly's
    # default, leaving room for the group labels in between.
    for i, pop in enumerate(present_populations):
        xaxis_suffix = "" if i == 0 else str(i + 1)
        fig.layout[f"xaxis{xaxis_suffix}"].title = dict(
            text=pop, font=dict(size=11, color="#374151"), standoff=28,
        )

    for i, pop in enumerate(present_populations):
        p_rows = stats_df.loc[stats_df["population"] == pop, "p_value"]
        p_text = f"p={p_rows.values[0]:.4f}" if not p_rows.empty and not pd.isna(p_rows.values[0]) else ""
        fig.for_each_annotation(
            lambda a, pop=pop, p_text=p_text: a.update(text=p_text)
            if a.text.split("=")[-1] == pop else a
        )
    return fig, present_populations


def render_single_population_boxplot(comparison_df: pd.DataFrame):
    """Boxplot of percentage distribution per population, colored directly
    by each population's own color -- used when NOT splitting by
    responder/non-responder (no group to compare against, so there's
    nothing for a background tint to distinguish from). No faceting
    either: with only one box per population, the x-axis category names
    already identify each population, so a 5-facet layout would just add
    empty structure around a single box per panel. Fill uses
    GROUP_FILL_ALPHA and a darkened outline/median, matching every other
    boxplot in the dashboard for visual consistency."""
    present_populations = [p for p in POPULATIONS if p in comparison_df["population"].unique()]
    fig = px.box(
        comparison_df, x="population", y="percentage", color="population",
        color_discrete_map=POP_COLORS, points="outliers",
        category_orders={"population": present_populations},
        labels={"percentage": "% of total cells", "population": ""},
    )
    for trace in fig.data:
        pop_color = POP_COLORS.get(trace.name)
        if pop_color:
            outline = darken_hex(pop_color)
            trace.fillcolor = hex_to_rgba(pop_color, GROUP_FILL_ALPHA)
            trace.line.color = outline
            trace.marker.color = outline
    fig.update_layout(showlegend=False)
    fig.update_yaxes(gridcolor="#000000", gridwidth=0.5, zeroline=False)
    return fig


def render_population_key(populations: list[str]):
    """Renders a small colored-swatch key for the population background
    tints, as plain Streamlit HTML entirely OUTSIDE the Plotly figure --
    not part of its (auto-positioned) legend. This guarantees it can
    never overlap the plot area, unlike an in-figure legend, which can
    crowd or overlap the chart on narrow viewports or with many entries.

    `populations` should be the exact list render_boxplot/
    render_n_group_boxplot returned as their second value (the
    present_populations they actually faceted, already in graphed order)
    -- not independently recomputed from a stats table. Passing
    results["population"].tolist() used to work by coincidence (both
    computations happened to apply the same canonical-order filter), but
    that's fragile: two independent computations of "which populations,
    in what order" can silently drift apart if either one's logic changes
    later. This function still deduplicates and canonical-orders its
    input defensively, but the caller contract is now "pass what was
    actually graphed."

    Swatch color uses POPULATION_BACKGROUND_OPACITY, matching the actual
    facet background opacity exactly -- rather than a solid, fully opaque
    dot, which would show a noticeably more saturated color than what the
    chart itself displays."""
    seen = [p for p in POPULATIONS if p in populations]
    swatches = "".join(
        f"<span style='display:inline-flex; align-items:center; margin-right:14px;'>"
        f"<span style='display:inline-block; width:10px; height:10px; "
        f"border-radius:50%; background:{hex_to_rgba(POP_COLORS[pop], POPULATION_BACKGROUND_OPACITY)}; "
        f"margin-right:5px;'></span>"
        f"<span style='font-size:12px; color:#4B5563;'>{pop}</span></span>"
        for pop in seen if pop in POP_COLORS
    )
    st.markdown(
        f"<div style='margin-top:6px; margin-bottom:2px;'>"
        f"<span style='font-size:12px; color:#6B7280; margin-right:10px;'>"
        f"Facet background:</span>{swatches}</div>",
        unsafe_allow_html=True,
    )


def _build_avg_bar_chart(avg_table, y_col, y_label, title, color_col, color_map, category_orders, is_comparison, yaxis_tickformat=None):
    """Builds one average-count or average-percentage bar chart. When
    is_comparison is True (color_col is a comparison group, not
    population), facets by population instead of putting population on
    the shared x-axis, deliberately WITHOUT matches='y' (unlike the
    boxplot), so each population gets its own independently-scaled
    y-axis. This matters concretely: between-population variation (e.g.
    b_cell ~10k cells vs cd4_t_cell ~30k) is far larger than
    within-population variation across comparison groups (often under
    1-2%), so a shared y-axis makes real, correctly-computed differences
    visually indistinguishable. Verified directly against real data
    (per-day averages genuinely differ, just by a small amount) before
    concluding this was a display problem, not a computation bug.

    Both branches use marker_opacity=GROUP_FILL_ALPHA, the same value
    the boxplot fills use (see render_boxplot), so the same population
    or group colors read at the same visual intensity whether they show
    up in a bar or a box. Plotly's Bar trace has a single marker.opacity
    that fades the fill and the outline together (unlike Box, which has
    a separate fillcolor property); the matching-color outline trick
    here still applies at that same faded opacity rather than staying
    fully solid, which is a minor difference from the boxplot's
    darkened, fully opaque outline, but keeps this change scoped to
    opacity specifically rather than also changing the outline's color."""
    if is_comparison:
        fig = px.bar(
            avg_table, x=color_col, y=y_col, color=color_col,
            color_discrete_map=color_map, facet_col="population", facet_col_wrap=5,
            facet_col_spacing=0.045,
            title=title, labels={y_col: y_label}, category_orders=category_orders,
        )
        fig.update_yaxes(matches=None, gridcolor="#000000", gridwidth=0.5)
        if yaxis_tickformat:
            fig.update_yaxes(tickformat=yaxis_tickformat)
        # showticklabels=False only hides the per-bar tick labels
        # (Responder/Non-responder under each bar, redundant with the
        # legend). It does NOT touch the axis title, which Plotly
        # auto-populates from the x= column name independently for each
        # of the 5 facets, so without title_text="" every facet showed
        # its own "response_label"/"cohort"/"timepoint" title and, packed
        # into a narrow shared row, they visually ran together into
        # unreadable repeated text. Confirmed directly against a real
        # screenshot: this was a second, separate bug from the population
        # facet titles above the chart, not the same issue restated.
        fig.update_xaxes(showticklabels=False, title_text="")
        fig.update_traces(marker_line_width=1.5, marker_opacity=GROUP_FILL_ALPHA)
        for trace in fig.data:
            trace.marker.line.color = trace.marker.color
        # Facet titles (population names) are staggered onto two
        # alternating heights rather than relying on font size and
        # spacing alone. With 5 titles packed into roughly half the page
        # width (this chart sits next to its sibling), every adjacent
        # pair competes for the same narrow horizontal band, and a font
        # small enough to reliably fit the longest names ("cd8_t_cell",
        # "cd4_t_cell", 10 characters) either way stops being legible.
        # Staggering means an adjacent pair no longer needs to fit
        # side-by-side in that band at all: each title only needs to
        # clear the *next* one over (two positions away), which has a
        # full facet width of horizontal room between them. Confirmed
        # directly that shrinking font size and adding facet_col_spacing
        # alone were not enough (an earlier attempt still overlapped on
        # the longest adjacent pair), so this is a structurally different
        # fix, not another size tweak.
        for i, a in enumerate(fig.layout.annotations):
            a.text = a.text.split("=")[-1]
            a.font = dict(size=11)
            a.y = 1.09 if i % 2 == 0 else 1.02
        fig.update_layout(margin=dict(t=90))
    else:
        fig = px.bar(
            avg_table, x="population", y=y_col, color=color_col,
            color_discrete_map=color_map, barmode="group",
            title=title, labels={y_col: y_label}, category_orders=category_orders,
        )
        fig.update_yaxes(gridcolor="#000000", gridwidth=0.5)
        if yaxis_tickformat:
            fig.update_yaxes(tickformat=yaxis_tickformat)
        fig.update_layout(showlegend=False)
        fig.update_traces(marker_opacity=GROUP_FILL_ALPHA)
    return fig


def render_avg_charts(avg_table: pd.DataFrame, color_col: str, color_map: dict, key_prefix: str, group_order: list | None = None):
    """Average cell count / average percentage bar charts, colored
    consistently by whatever column is being compared (population, or
    response, or cohort, or timepoint). All bars, whether colored by
    population or by a comparison group, use marker_opacity=GROUP_FILL_ALPHA:
    the same fill opacity the boxplot uses, so the same color reads at
    the same visual intensity in both chart types (no hatch pattern: an
    earlier version used pattern_shape as a colorblind safety measure,
    but Plotly's default pattern assigns the *first* category an
    empty/solid pattern and the *second* a diagonal-hatch fill, which
    reads as "not fully filled in" rather than as a deliberate texture,
    removed in favor of just keeping the Okabe-Ito color choices
    themselves colorblind-safe, which they already are).

    group_order pins the legend/bar order for color_col explicitly. This
    matters: without it, Plotly Express sorts a string column's values
    alphabetically by default, which silently put "Non-responder" before
    "Responder" here (while the distribution boxplot elsewhere used an
    explicit Responder-first order), a real inconsistency, not a
    stylistic choice, caught by comparing the two side by side."""
    # Uses only the populations actually present in avg_table, not the
    # full POPULATIONS constant. The same bug fixed earlier for
    # render_boxplot's facets applies identically here: Plotly Express
    # creates one x-axis category per category_orders entry even with
    # zero matching rows, so a population-filtered cohort (e.g. just
    # b_cell + cd4_t_cell selected) reserved x-axis space for all 5
    # populations, leaving empty gaps where cd8_t_cell/nk_cell/monocyte
    # would be -- confirmed directly against a real filtered chart.
    present_populations = [p for p in POPULATIONS if p in avg_table["population"].unique()]
    category_orders = {"population": present_populations}
    if group_order:
        category_orders[color_col] = group_order
    is_comparison = color_col != "population"

    chart_col1, chart_col2 = st.columns(2)
    with chart_col1:
        fig_count = _build_avg_bar_chart(
            avg_table, "avg_count", "avg. cell count", "Average number of cells",
            color_col, color_map, category_orders, is_comparison, yaxis_tickformat=",",
        )
        st.plotly_chart(fig_count, width='stretch', key=f"{key_prefix}_avg_count_chart")
    with chart_col2:
        fig_pct = _build_avg_bar_chart(
            avg_table, "avg_percentage", "avg. % of total cells", "Average relative frequency",
            color_col, color_map, category_orders, is_comparison,
        )
        st.plotly_chart(fig_pct, width='stretch', key=f"{key_prefix}_avg_pct_chart")


def render_cohort_summary_block(df: pd.DataFrame):
    """Cohort summary card: sample/subject counts + breakdowns by
    project/condition/treatment/response/sex. Shared across every
    single-cohort tab (Default, Responder vs Non-responder, By
    Population, By Date) so this section always appears in the same
    place with the same shape -- Custom mode uses
    combined_cohort_breakdown_table instead, since it has multiple
    independent cohorts rather than one."""
    summary = get_cohort_summary(df)
    with st.container(border=True):
        c1, c2 = st.columns(2)
        c1.metric("Samples in cohort", f"{summary['n_samples']:,}")
        c2.metric("Subjects in cohort", f"{summary['n_subjects']:,}")

        b1, b2, b3, b4, b5 = st.columns(5)
        with b1:
            st.write("**By project**")
            st.dataframe(summary["samples_per_project"], hide_index=True, width='stretch')
        with b2:
            st.write("**By condition**")
            st.dataframe(summary["subjects_by_condition"], hide_index=True, width='stretch')
        with b3:
            st.write("**By treatment**")
            st.dataframe(summary["subjects_by_treatment"], hide_index=True, width='stretch')
        with b4:
            st.write("**By response**")
            st.dataframe(summary["subjects_by_response"], hide_index=True, width='stretch')
        with b5:
            st.write("**By sex**")
            st.dataframe(summary["subjects_by_sex"], hide_index=True, width='stretch')


def combined_cohort_breakdown_table(cohort_dfs_by_label: dict, breakdown_key: str) -> pd.DataFrame:
    """Combines get_cohort_summary's breakdown table across 2+ cohorts
    into ONE table with a 'cohort' column, instead of one small table per
    cohort per dimension -- used by Custom mode's summary section so it
    stays at 4 breakdown tables total (same as every other tab), not
    4 x (number of cohorts)."""
    frames = []
    for label, df in cohort_dfs_by_label.items():
        summary = get_cohort_summary(df)
        piece = summary[breakdown_key].copy()
        piece.insert(0, "cohort", label)
        frames.append(piece)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def render_avg_table_and_charts(
    avg_table: pd.DataFrame, key_prefix: str, color_col: str, color_map: dict,
    group_order: list | None = None, download_filename: str = "average_cell_counts.csv",
):
    """Average cell count/percentage table (+ download) and bar charts,
    given an already-built avg_table. Shared rendering so this section
    looks and behaves identically across every tab -- the table itself
    (not just the charts) is always shown and always downloadable, which
    an earlier version of By Date and Custom mode omitted."""
    if avg_table.empty:
        st.info("No data available for this cohort.")
        return
    display_cols = {
        "avg_count": "avg. number of cells",
        "avg_percentage": "avg. % of total cells",
        "n_samples": "n samples",
    }
    avg_display = avg_table.rename(columns=display_cols)
    st.dataframe(
        style_population_column(avg_display), width='stretch', hide_index=True,
        column_config={"avg. number of cells": st.column_config.NumberColumn(format="localized")},
    )
    download_csv_button(avg_table, download_filename, "Download this table as CSV", key=f"{key_prefix}_avg_download")
    st.write("")
    render_avg_charts(avg_table, color_col, color_map, key_prefix=key_prefix, group_order=group_order)


def render_frequency_table_block(freq_df: pd.DataFrame, key_prefix: str, download_filename: str = "frequency_table.csv"):
    """Sample x population frequency table (+ download). Shared rendering
    so this section appears identically across every tab."""
    # Not population-colored: this table can have thousands of rows (one
    # per sample x population), which exceeds pandas Styler's cell-count
    # limit and, more importantly, wouldn't be very readable as row-by-row
    # color banding anyway. Coloring is reserved for the short, scannable
    # summary tables instead.
    st.dataframe(
        freq_df, width='stretch', height=300,
        column_config={
            "count": st.column_config.NumberColumn(format="localized"),
            "total_count": st.column_config.NumberColumn(format="localized"),
        },
    )
    download_csv_button(freq_df, download_filename, "Download this table as CSV", key=f"{key_prefix}_freq_download")


def render_stats_messages(status: str, results: pd.DataFrame | None):
    """Shared rendering of the four run_stats_test_safe cases."""
    if status == "no_samples":
        st.info("No samples match the selected filters.")
        return
    if status == "no_response_data":
        st.info(
            "No responder/non-responder data available for the selected "
            "filters. The selected subjects were not treated, so response "
            "doesn't apply to them."
        )
        return
    for _, row in results.iterrows():
        if row["status"] == "no_individuals":
            st.warning(f"**{row['population']}**: No individuals meet these criteria.")
        elif row["small_n_warning"]:
            st.warning(
                f"**{row['population']}**: small sample size "
                f"(n={row['n_responders']} responders, n={row['n_non_responders']} "
                f"non-responders). p-value may be unreliable (threshold: n<{SMALL_N_THRESHOLD})."
            )


def render_comparison_stats_table(
    results: pd.DataFrame, key_prefix: str, download_filename: str,
    group_col: str | None = "population",
):
    """Renders a comparison stats table (Mann-Whitney/Wilcoxon results)
    with a toggle between two orderings: grouped by group_col (default --
    e.g. all b_cell rows together, then all cd8_t_cell rows, in canonical
    population order), or with the most statistically significant results
    surfaced at the top regardless of group. Both are genuinely useful for
    different purposes -- grouped-by-type answers "what's going on with
    this population specifically," surfaced-by-significance answers
    "what's the strongest finding here" -- so this is an explicit choice
    rather than one silently picked for the user. Shared across every
    comparison tab so the toggle behaves identically everywhere."""
    sort_significant = st.checkbox(
        "Bring significant results to the top", key=f"{key_prefix}_sig_sort",
        help="Off (default): grouped by cell type. On: sorted by "
             "statistical significance, most significant first, "
             "regardless of cell type.",
    )
    display_results = results.drop(columns=["status", "small_n_warning"]).copy()

    if sort_significant:
        display_results = display_results.sort_values(
            ["significant_bonferroni", "p_value_bonferroni"], ascending=[False, True],
        )
    elif group_col and group_col in display_results.columns:
        sort_key = display_results[group_col]
        if group_col == "population":
            sort_key = pd.Categorical(sort_key, categories=POPULATIONS, ordered=True)
        display_results = (
            display_results.assign(_sort_key=sort_key)
            .sort_values(["_sort_key", "p_value_bonferroni"])
            .drop(columns=["_sort_key"])
        )
    else:
        display_results = display_results.sort_values("p_value_bonferroni")

    display_results = display_results.reset_index(drop=True)
    display_results["p_value"] = display_results["p_value"].apply(format_pvalue)
    display_results["p_value_bonferroni"] = display_results["p_value_bonferroni"].apply(format_pvalue)

    if group_col == "population" and "population" in display_results.columns:
        st.dataframe(style_population_column(display_results), width='stretch', hide_index=True)
    else:
        st.dataframe(display_results, width='stretch', hide_index=True)
    download_csv_button(display_results, download_filename, "Download this table as CSV", key=f"{key_prefix}_download")


def render_n_group_messages(status: str, results: pd.DataFrame | None):
    """Shared rendering of compare_n_groups' / compare_populations_paired's
    status cases -- generalizes render_stats_messages/the old
    render_cohort_comparison_messages to an arbitrary number of groups."""
    if status == "insufficient_groups":
        st.info(
            "Fewer than 2 of the selected groups have any matching samples. "
            "Select at least 2 groups with data, or broaden your filters."
        )
        return
    has_population = results is not None and "population" in results.columns
    for _, row in results.iterrows():
        pop_prefix = f"**{row['population']}**, " if has_population else ""
        if row["status"] == "no_individuals":
            st.warning(
                f"{pop_prefix}{row['group_a']} vs {row['group_b']}: not present "
                f"in both groups (n={row['n_a']} in {row['group_a']}, "
                f"n={row['n_b']} in {row['group_b']})."
            )
        elif row["small_n_warning"]:
            st.warning(
                f"{pop_prefix}{row['group_a']} vs {row['group_b']}: small sample "
                f"size (n={row['n_a']}, n={row['n_b']}). p-value may be unreliable "
                f"(threshold: n<{SMALL_N_THRESHOLD})."
            )


def render_n_group_boxplot(group_dfs: dict, stats_df: pd.DataFrame, colors: list, x_label: str = "group"):
    """Boxplot for 2 or more arbitrary groups, faceted by population --
    generalizes the old render_cohort_comparison_boxplot (which was fixed
    at exactly 2 cohorts) to any number of groups via render_boxplot's own
    N-group support. Returns (fig, present_populations), same as
    render_boxplot -- this just forwards it."""
    frames = []
    for label, df in group_dfs.items():
        d = df.copy()
        d[x_label] = label
        frames.append(d)
    combined = pd.concat(frames, ignore_index=True)
    labels = list(group_dfs.keys())
    return render_boxplot(
        combined, stats_df, x_col=x_label, group_order=labels,
        group_colors=colors[:len(labels)],
    )


def render_population_comparison_boxplot(df: pd.DataFrame, populations: list):
    """Boxplot comparing 2+ selected populations directly against each
    other, for the paired population-vs-population comparison (see
    compare_populations_paired). Visually identical to
    render_single_population_boxplot, restricted to just the selected
    populations -- pairwise p-values are shown in the stats table below
    rather than as on-chart significance brackets, since there's no
    faceting here to anchor per-pair annotations to."""
    subset = df[df["population"].isin(populations)]
    return render_single_population_boxplot(subset)


def age_range_widgets(age_lo: int, age_hi: int, key_prefix: str) -> tuple[int, int] | None:
    """
    Renders BOTH a slider and two number inputs (min/max) for age, kept in
    sync via session_state. Dragging the slider updates the number inputs
    and vice versa. Returns the (min, max) tuple, or None if it still
    equals the full range (i.e. no restriction), matching filter_dataset's
    "None means no filter" contract.
    """
    range_key = f"{key_prefix}_age_range"
    min_key = f"{key_prefix}_age_min_input"
    max_key = f"{key_prefix}_age_max_input"

    if range_key not in st.session_state:
        st.session_state[range_key] = (age_lo, age_hi)
    if min_key not in st.session_state:
        st.session_state[min_key] = st.session_state[range_key][0]
    if max_key not in st.session_state:
        st.session_state[max_key] = st.session_state[range_key][1]

    def _sync_from_slider():
        lo, hi = st.session_state[range_key]
        st.session_state[min_key] = lo
        st.session_state[max_key] = hi

    def _sync_from_inputs():
        lo, hi = st.session_state[min_key], st.session_state[max_key]
        if lo > hi:
            lo, hi = hi, lo
        st.session_state[range_key] = (lo, hi)

    st.slider("Age range", age_lo, age_hi, key=range_key, on_change=_sync_from_slider)
    c_min, c_max = st.columns(2)
    c_min.number_input(
        "Min age", min_value=age_lo, max_value=age_hi,
        key=min_key, on_change=_sync_from_inputs,
    )
    c_max.number_input(
        "Max age", min_value=age_lo, max_value=age_hi,
        key=max_key, on_change=_sync_from_inputs,
    )

    current = st.session_state[range_key]
    return current if current != (age_lo, age_hi) else None


def cohort_filter_widgets(full: pd.DataFrame, key_prefix: str, exclude: set | None = None) -> dict:
    """Renders the standard set of filter widgets, with keys prefixed so
    multiple independent instances can coexist on the same page without
    Streamlit's duplicate-widget-ID error. `exclude` skips specific
    fields entirely (not filtered on, not rendered) -- used when a field
    is the comparison axis itself for a given mode (e.g. population for
    "By Population" mode, time_from_treatment_start for "By Date"), so
    the same variable can't be both a pre-filter and the thing being
    compared at once."""
    exclude = exclude or set()
    condition_opts = sorted(full["condition"].unique())
    treatment_opts = sorted(full["treatment"].unique())
    sample_type_opts = sorted(full["sample_type"].unique())
    sex_opts = sorted(full["sex"].unique())
    project_opts = sorted(full["project"].unique())
    time_opts = sorted(full["time_from_treatment_start"].unique())
    age_lo, age_hi = int(full["age"].min()), int(full["age"].max())

    filters = {}
    c1, c2 = st.columns(2)
    if "condition" not in exclude:
        filters["condition"] = c1.multiselect("Condition", condition_opts, key=f"{key_prefix}_condition") or None
    if "treatment" not in exclude:
        filters["treatment"] = c2.multiselect("Treatment", treatment_opts, key=f"{key_prefix}_treatment") or None
    c3, c4 = st.columns(2)
    if "sample_type" not in exclude:
        filters["sample_type"] = c3.multiselect("Sample type", sample_type_opts, key=f"{key_prefix}_sample_type") or None
    if "sex" not in exclude:
        filters["sex"] = c4.multiselect("Sex", sex_opts, key=f"{key_prefix}_sex") or None
    c5, c6 = st.columns(2)
    if "project" not in exclude:
        filters["project"] = c5.multiselect("Project", project_opts, key=f"{key_prefix}_project") or None
    if "time_from_treatment_start" not in exclude:
        filters["time_from_treatment_start"] = c6.multiselect("Time from treatment start", time_opts, key=f"{key_prefix}_time") or None
    if "response" not in exclude:
        filters["response"] = st.multiselect("Response", ["yes", "no"], key=f"{key_prefix}_response") or None
    if "population" not in exclude:
        filters["population"] = st.multiselect(
            "Cell type (population)", POPULATIONS, key=f"{key_prefix}_population",
            help="Leave empty to include all 5. Selecting specific populations "
                 "only narrows which rows are shown/tested. Percentages still "
                 "reflect each sample's full cell count, not just the "
                 "selected populations.",
        ) or None
    if "age_range" not in exclude:
        filters["age_range"] = age_range_widgets(age_lo, age_hi, key_prefix)

    return filters


FILTER_KEY_SUFFIXES = [
    "condition", "treatment", "sample_type", "sex", "project", "response", "time", "population",
]


def reset_filter_state(key_prefix: str, age_lo: int, age_hi: int):
    """Clears every filter widget under key_prefix back to its default
    (empty selection / full age range) by writing directly to
    session_state, then triggers a rerun so the widgets pick up the reset
    values on their next render."""
    for suffix in FILTER_KEY_SUFFIXES:
        key = f"{key_prefix}_{suffix}"
        if key in st.session_state:
            st.session_state[key] = []
    if f"{key_prefix}_age_range" in st.session_state:
        st.session_state[f"{key_prefix}_age_range"] = (age_lo, age_hi)
        st.session_state[f"{key_prefix}_age_min_input"] = age_lo
        st.session_state[f"{key_prefix}_age_max_input"] = age_hi


def filters_with_expander(full: pd.DataFrame, key_prefix: str, label: str = "Filters", exclude: set | None = None) -> dict:
    """Wraps cohort_filter_widgets in a collapsible expander, shows a
    colored badge summarizing how many filters are active even while
    collapsed, and offers a one-click reset back to the full dataset."""
    age_lo, age_hi = int(full["age"].min()), int(full["age"].max())

    with st.expander(label, expanded=False):
        filters = cohort_filter_widgets(full, key_prefix, exclude=exclude)

    active_n = sum(1 for v in filters.values() if v)
    badge_col, reset_col = st.columns([4, 1])
    with badge_col:
        if active_n:
            st.markdown(
                f"<span class='filter-badge active'>{active_n} filter"
                f"{'s' if active_n != 1 else ''} active</span>",
                unsafe_allow_html=True,
            )
        else:
            st.markdown(
                "<span class='filter-badge inactive'>No filters. Showing the full dataset</span>",
                unsafe_allow_html=True,
            )
    with reset_col:
        if active_n:
            st.button(
                "Reset", key=f"{key_prefix}_reset_btn",
                on_click=reset_filter_state, args=(key_prefix, age_lo, age_hi),
            )
    return filters


def _db_has_data() -> bool:
    """Checks not just that cell_counts.db exists, but that it actually
    contains rows. A file that exists but is empty or only has an empty
    schema (e.g. from a build that was interrupted between
    sqlite3.connect creating the file and the inserts finishing, such as
    a container restart, resource limit, or crash mid-build) would
    otherwise be silently treated as "already built" forever afterward,
    since load_data.py's build_database() creates the .db file the
    instant it opens a connection, well before any data is written.
    Confirmed directly: this exact failure mode produced a live
    "0 subjects, 0 samples, 0 projects" app with a ValueError crashing
    on int(full["age"].min()) (min of an empty column is NaN), because
    the old check only asked "does the file exist," not "does it have
    data," and a stale empty file always answered yes to the former."""
    if not DB_PATH.exists():
        return False
    try:
        conn = sqlite3.connect(DB_PATH)
        try:
            n = conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
            return n > 0
        finally:
            conn.close()
    except sqlite3.Error:
        # Missing table, corrupt file, etc. Anything that fails this
        # check should be treated the same as "needs rebuilding," not
        # as an unrelated crash.
        return False


if not _db_has_data():
    if CSV_PATH.exists():
        # Self-initialize: this matters for cloud deployment (e.g.
        # Streamlit Community Cloud) where there's no separate "run
        # load_data.py first" step available. Locally, `make pipeline` /
        # `python load_data.py` still works fine and this branch just
        # won't fire since the .db will already exist and have data.
        with st.spinner("First run: building the database from cell-count.csv..."):
            from load_data import build_database
            build_database()
    else:
        st.error(
            f"Neither {DB_PATH.name} nor {CSV_PATH.name} was found. "
            "Make sure cell-count.csv is included in the deployment."
        )
        st.stop()

mtime = _csv_mtime()
full = cached_full_dataset(mtime)

# ---------- sidebar ----------
st.sidebar.markdown(
    "<div style='padding-bottom: 8px;'>"
    "<span style='font-size: 22px; font-weight: 700; color: #1F6F78;'>LOBLAW BIO</span><br/>"
    "<span style='font-size: 12px; color: #6B7280;'>Clinical Trial Analytics</span>"
    "</div>",
    unsafe_allow_html=True,
)
st.sidebar.caption(
    "Immune cell population data from Loblaw Bio's clinical trial, "
    "covering melanoma, carcinoma, and healthy cohorts."
)
st.sidebar.divider()
st.sidebar.markdown("**Dataset**")
st.sidebar.write(f"{full['subject_id'].nunique():,} subjects")
st.sidebar.write(f"{full['sample'].nunique():,} samples")
st.sidebar.write(f"{full['project'].nunique()} projects")
st.sidebar.write(f"{len(POPULATIONS)} cell populations tracked")

# ---------- header ----------
st.title("Loblaw Bio")
st.caption("Clinical trial cell population analysis")
st.markdown("<div class='accent-bar'></div>", unsafe_allow_html=True)

tab_default, tab_resp, tab_pop, tab_date, tab_custom = st.tabs(
    ["Default", "Responder vs Non-responder", "By Population", "By Date", "Custom"]
)
# Using st.tabs() (not a radio + if/elif) is deliberate: Streamlit renders
# every tab's content on every rerun (just CSS-hides the inactive ones),
# so widget state naturally survives switching tabs. A radio-driven
# if/elif only instantiates the currently-selected branch's widgets each
# run -- Streamlit clears session_state for any widget that disappears
# from the script that way, so filter selections were being silently
# reset every time you switched views. Confirmed directly (not assumed):
# selecting specific populations in By Population mode, switching to
# Default, then switching back reset the selection under the old radio
# structure; verified this doesn't happen with tabs.


# ---------- Default: single-cohort exploration (no comparison) ----------
with tab_default:
    st.subheader("Build your own cohort")
    st.caption(
        "Filter by any combination of variables in the dataset. Cohort "
        "counts, average cell counts, and the frequency table all "
        "recompute live. Switch to one of the other tabs above to "
        "compare groups against each other."
    )

    filters = filters_with_expander(full, "default")
    filtered = filter_dataset(full, **filters)

    st.divider()

    if filtered.empty:
        st.info("No samples match the selected filters.")
    else:
        render_cohort_summary_block(filtered)

        st.write("")
        with st.container(border=True):
            st.write("**Average cell count and relative frequency for this cohort**")
            avg_table = get_population_averages(filtered)
            render_avg_table_and_charts(avg_table, "default", "population", POP_COLORS)

        st.write("")
        with st.container(border=True):
            st.write("**Frequency table for this cohort**")
            freq_filtered = get_filtered_frequency_table(filtered)
            render_frequency_table_block(freq_filtered, "default")

        st.write("")
        with st.container(border=True):
            st.write("**Population distribution for this cohort**")
            st.caption(
                "Percentage distribution of each cell population, colored "
                "by population. Not a statistical comparison -- switch to "
                "one of the other tabs above to test 2 or more groups "
                "against each other."
            )
            fig = render_single_population_boxplot(filtered)
            st.plotly_chart(fig, width='stretch', key="default_boxplot", config=PLOTLY_CONFIG)


# ---------- Responder vs Non-responder ----------
with tab_resp:
    st.subheader("Responder vs non-responder comparison")
    st.caption(
        "Build a cohort using any combination of filters, then compare "
        "responders against non-responders within it. Response isn't "
        "offered as a pre-filter here -- it's the comparison axis for "
        "this view."
    )

    filters = filters_with_expander(full, "resp_mode", exclude={"response"})
    filtered = filter_dataset(full, **filters)

    st.divider()

    if filtered.empty:
        st.info("No samples match the selected filters.")
    else:
        render_cohort_summary_block(filtered)

        st.write("")
        st.caption(
            "Exploratory result: this test ran against whatever cohort you "
            "built above, not a pre-specified comparison. This tool lets you "
            "filter and re-test in effectively unlimited ways; p-values from "
            "ad hoc slices carry less evidentiary weight than a single "
            "pre-registered comparison (like the fixed miraclib+PBMC "
            "comparison documented in the README) and shouldn't be read as "
            "confirmatory on their own."
        )
        status, results = run_stats_test_safe(filtered)
        render_stats_messages(status, results)

        if status == "ok":
            with st.expander("Confounder check: is response balanced across project / sex?", expanded=False):
                st.caption(
                    "A high p-value here means no evidence of imbalance was "
                    "found, not that confounding is proven absent -- this is a "
                    "simple heuristic (p > 0.05), not a formal equivalence test."
                )
                for stratify_col, label in [("project", "Project"), ("sex", "Sex")]:
                    balance = check_group_balance(filtered, group_col="response", stratify_col=stratify_col)
                    if balance["p_value"] is None:
                        st.write(f"**{label}**: only one level present in this cohort, check doesn't apply.")
                    else:
                        icon = "\u2705" if balance["balanced"] else "\u26a0\ufe0f"
                        verdict = "no evidence of imbalance" if balance["balanced"] else "possible imbalance, interpret with extra caution"
                        st.write(f"{icon} **{label}**: p={balance['p_value']:.4f} ({verdict})")
                        st.dataframe(balance["contingency_table"], width='stretch')

            st.write("")
            with st.container(border=True):
                st.write("**Average cell count and relative frequency**")
                avg_table = get_population_averages(filtered, split_by_response=True)
                render_avg_table_and_charts(
                    avg_table, "resp_mode", "response_label", RESPONSE_COLORS,
                    group_order=["Responder", "Non-responder"],
                    download_filename="responder_averages.csv",
                )

            st.write("")
            with st.container(border=True):
                st.write("**Frequency table for this cohort**")
                freq_filtered = get_filtered_frequency_table(filtered)
                render_frequency_table_block(freq_filtered, "resp_mode")

            st.write("")
            with st.container(border=True):
                st.write("**Distribution comparison**")
                fig, graphed_populations = render_boxplot(
                    filtered, results, group_order=["Responder", "Non-responder"],
                    group_colors=[RESPONSE_COLORS["Responder"], RESPONSE_COLORS["Non-responder"]],
                )
                st.plotly_chart(fig, width='stretch', key="resp_boxplot", config=PLOTLY_CONFIG)
                render_population_key(graphed_populations)

            st.write("")
            with st.container(border=True):
                st.write("**Mann-Whitney U results** (Bonferroni-corrected across populations tested)")
                render_comparison_stats_table(results, "resp_stats", "responder_comparison_stats.csv")


# ---------- By Population: compare populations directly against each other ----------
with tab_pop:
    st.subheader("Compare cell populations directly")
    st.caption(
        "Build a cohort using any combination of filters, then compare 2 "
        "or more cell populations' relative frequencies directly against "
        "each other within it. This uses a paired test (Wilcoxon "
        "signed-rank), not the unpaired Mann-Whitney used in the other "
        "comparison views -- two populations' percentages from the same "
        "sample aren't independent of each other (see README, "
        "'Statistical approach', for why)."
    )

    filters = filters_with_expander(full, "pop_mode", exclude={"population"})
    filtered = filter_dataset(full, **filters)

    st.divider()

    if filtered.empty:
        st.info("No samples match the selected filters.")
    else:
        render_cohort_summary_block(filtered)

        st.write("")
        with st.container(border=True):
            st.write("**Average cell count and relative frequency**")
            avg_table = get_population_averages(filtered)
            render_avg_table_and_charts(avg_table, "pop_mode", "population", POP_COLORS)

        st.write("")
        with st.container(border=True):
            st.write("**Frequency table for this cohort**")
            freq_filtered = get_filtered_frequency_table(filtered)
            render_frequency_table_block(freq_filtered, "pop_mode")

        st.write("")
        selected_pops_raw = st.multiselect(
            "Populations to compare", POPULATIONS, default=POPULATIONS,
            key="pop_mode_selected_populations",
        )
        # Sorted to the canonical population order, not the order clicked --
        # st.multiselect returns selections in click order, so without this
        # a user who clicks nk_cell before b_cell would see that order
        # reflected in the comparison, purely as an artifact of click order.
        selected_pops = [p for p in POPULATIONS if p in selected_pops_raw]
        if len(selected_pops) < 2:
            st.info("Select at least 2 populations to compare.")
        else:
            st.caption(
                "Exploratory result: this cohort was built from whatever "
                "filters you chose above, not a pre-specified comparison. "
                "p-values from ad hoc slices carry less evidentiary weight "
                "than a single pre-registered comparison and shouldn't be "
                "read as confirmatory on their own."
            )
            status, results = compare_populations_paired(filtered, selected_pops)
            render_n_group_messages(status, results)

            if status == "ok":
                with st.container(border=True):
                    st.write("**Population distributions**")
                    fig = render_population_comparison_boxplot(filtered, selected_pops)
                    st.plotly_chart(fig, width='stretch', key="pop_mode_boxplot", config=PLOTLY_CONFIG)

                st.write("")
                with st.container(border=True):
                    st.write("**Wilcoxon signed-rank results** (Bonferroni-corrected across pairs tested)")
                    render_comparison_stats_table(
                        results, "pop_mode_stats", "population_comparison_stats.csv", group_col="group_a",
                    )


# ---------- By Date: compare timepoints directly, any number at once ----------
with tab_date:
    st.subheader("Compare timepoints directly")
    st.caption(
        "Build a cohort using any combination of filters, then compare 2 "
        "or more timepoints against each other, per population. Selecting "
        "3+ timepoints tests every pair, Bonferroni-corrected across all "
        "pairs and populations tested together -- this correction gets "
        "stricter fast as more timepoints are selected at once (see "
        "compare_n_groups in analysis.py)."
    )

    filters = filters_with_expander(full, "date_mode", exclude={"time_from_treatment_start"})
    filtered = filter_dataset(full, **filters)

    st.divider()

    if filtered.empty:
        st.info("No samples match the selected filters.")
    else:
        render_cohort_summary_block(filtered)

        st.write("")
        time_opts = sorted(full["time_from_treatment_start"].unique())
        selected_times_raw = st.multiselect(
            "Timepoints to compare", time_opts, default=time_opts,
            key="date_mode_selected_times",
        )
        # Sorted chronologically, not in click order -- st.multiselect
        # returns selections in the order they were clicked, so without
        # this a user who clicks Day 14 before Day 0 would see that
        # non-chronological order reflected in the comparison below.
        selected_times = [t for t in time_opts if t in selected_times_raw]
        if len(selected_times) < 2:
            st.info("Select at least 2 timepoints to compare.")
        else:
            st.caption(
                "Exploratory result: these timepoints were selected from "
                "whatever cohort you built above, not a pre-specified "
                "comparison. p-values from ad hoc slices carry less "
                "evidentiary weight than a single pre-registered comparison "
                "and shouldn't be read as confirmatory on their own."
            )
            group_dfs = {
                f"Day {t}": filter_dataset(filtered, time_from_treatment_start=[t])
                for t in selected_times
            }
            date_color_map = {
                label: N_GROUP_COLOR_SEQUENCE[i % len(N_GROUP_COLOR_SEQUENCE)]
                for i, label in enumerate(group_dfs.keys())
            }
            status, results = compare_n_groups(group_dfs)
            render_n_group_messages(status, results)

            if status == "ok":
                with st.container(border=True):
                    st.write("**Average cell count and relative frequency**")
                    avg_frames = []
                    for label, df in group_dfs.items():
                        avg = get_population_averages(df)
                        if not avg.empty:
                            avg_frames.append(avg.assign(timepoint=label))
                    avg_combined = pd.concat(avg_frames, ignore_index=True) if avg_frames else pd.DataFrame()
                    render_avg_table_and_charts(
                        avg_combined, "date_mode", "timepoint", date_color_map,
                        group_order=list(group_dfs.keys()),
                        download_filename="date_comparison_averages.csv",
                    )

                st.write("")
                with st.container(border=True):
                    st.write("**Frequency table**")
                    freq_frames = []
                    for label, df in group_dfs.items():
                        freq = get_filtered_frequency_table(df)
                        if not freq.empty:
                            freq_frames.append(freq.assign(timepoint=label))
                    freq_combined = pd.concat(freq_frames, ignore_index=True) if freq_frames else pd.DataFrame()
                    render_frequency_table_block(freq_combined, "date_mode", download_filename="date_comparison_frequency.csv")

                st.write("")
                with st.container(border=True):
                    st.write("**Distribution comparison**")
                    fig, graphed_populations = render_n_group_boxplot(
                        group_dfs, results, list(date_color_map.values()), x_label="timepoint",
                    )
                    st.plotly_chart(fig, width='stretch', key="date_mode_boxplot", config=PLOTLY_CONFIG)
                    render_population_key(graphed_populations)

                st.write("")
                with st.container(border=True):
                    st.write("**Mann-Whitney U results** (Bonferroni-corrected across all pairs and populations tested)")
                    render_comparison_stats_table(results, "date_mode_stats", "date_comparison_stats.csv")


# ---------- Custom: build 2-4 independent cohorts, any filters ----------
with tab_custom:
    st.subheader("Compare custom cohorts")
    st.caption(
        "Build 2 to 4 independent cohorts using any combination of "
        "filters, and compare their cell population averages and "
        "distributions directly against each other."
    )

    with st.expander("Add more cohorts (up to 4 total)", expanded=False):
        ec1, ec2 = st.columns(2)
        enable_c = ec1.checkbox("Enable Cohort C", key="custom_enable_c")
        enable_d = ec2.checkbox("Enable Cohort D", key="custom_enable_d")

    active_slots = ["A", "B"] + (["C"] if enable_c else []) + (["D"] if enable_d else [])

    cohort_labels = {}
    cohort_dfs_by_label = {}
    color_by_label = {}

    slot_cols = st.columns(len(active_slots))
    for i, slot in enumerate(active_slots):
        with slot_cols[i]:
            color = N_GROUP_COLOR_SEQUENCE[i % len(N_GROUP_COLOR_SEQUENCE)]
            st.markdown(
                f"<span class='cohort-dot' style='background:{color};'></span>"
                f"<strong>Cohort {slot}</strong>",
                unsafe_allow_html=True,
            )
            label = st.text_input(f"Label for Cohort {slot}", value=f"Cohort {slot}", key=f"custom_label_{slot}")
            slot_filters = filters_with_expander(full, f"custom_{slot}", label=f"Filters (Cohort {slot})")
            cohort_labels[slot] = label
            cohort_dfs_by_label[label] = filter_dataset(full, **slot_filters)
            color_by_label[label] = color

    st.divider()

    labels_used = [cohort_labels[s] for s in active_slots]
    if len(set(labels_used)) != len(labels_used):
        st.error("Cohort labels must be unique. Please give each cohort a different label above.")
    else:
        with st.container(border=True):
            for slot in active_slots:
                label = cohort_labels[slot]
                df = cohort_dfs_by_label[label]
                m1, m2 = st.columns(2)
                m1.metric(f"{label}: samples", f"{df['sample'].nunique():,}")
                m2.metric(f"{label}: subjects", f"{df['subject_id'].nunique():,}")

            b1, b2, b3, b4, b5 = st.columns(5)
            with b1:
                st.write("**By project**")
                st.dataframe(combined_cohort_breakdown_table(cohort_dfs_by_label, "samples_per_project"), hide_index=True, width='stretch')
            with b2:
                st.write("**By condition**")
                st.dataframe(combined_cohort_breakdown_table(cohort_dfs_by_label, "subjects_by_condition"), hide_index=True, width='stretch')
            with b3:
                st.write("**By treatment**")
                st.dataframe(combined_cohort_breakdown_table(cohort_dfs_by_label, "subjects_by_treatment"), hide_index=True, width='stretch')
            with b4:
                st.write("**By response**")
                st.dataframe(combined_cohort_breakdown_table(cohort_dfs_by_label, "subjects_by_response"), hide_index=True, width='stretch')
            with b5:
                st.write("**By sex**")
                st.dataframe(combined_cohort_breakdown_table(cohort_dfs_by_label, "subjects_by_sex"), hide_index=True, width='stretch')

        st.write("")
        st.caption(
            "Exploratory result: these cohorts were built from whatever "
            "filters you chose above, not a pre-specified comparison. "
            "p-values from ad hoc slices carry less evidentiary weight "
            "than a single pre-registered comparison and shouldn't be read "
            "as confirmatory on their own."
        )
        status, results = compare_n_groups(cohort_dfs_by_label)
        render_n_group_messages(status, results)

        if status == "ok":
            group_order = [cohort_labels[s] for s in active_slots]
            colors_in_order = [color_by_label[label] for label in group_order]

            st.write("")
            with st.container(border=True):
                st.write("**Average cell count and relative frequency**")
                avg_frames = []
                for label in group_order:
                    avg = get_population_averages(cohort_dfs_by_label[label])
                    if not avg.empty:
                        avg_frames.append(avg.assign(cohort=label))
                avg_combined = pd.concat(avg_frames, ignore_index=True) if avg_frames else pd.DataFrame()
                render_avg_table_and_charts(
                    avg_combined, "custom", "cohort", color_by_label,
                    group_order=group_order,
                    download_filename="custom_comparison_averages.csv",
                )

            st.write("")
            with st.container(border=True):
                st.write("**Frequency table**")
                freq_frames = []
                for label in group_order:
                    freq = get_filtered_frequency_table(cohort_dfs_by_label[label])
                    if not freq.empty:
                        freq_frames.append(freq.assign(cohort=label))
                freq_combined = pd.concat(freq_frames, ignore_index=True) if freq_frames else pd.DataFrame()
                render_frequency_table_block(freq_combined, "custom", download_filename="custom_comparison_frequency.csv")

            st.write("")
            with st.container(border=True):
                st.write("**Distribution comparison**")
                fig, graphed_populations = render_n_group_boxplot(cohort_dfs_by_label, results, colors_in_order, x_label="cohort")
                st.plotly_chart(fig, width='stretch', key="custom_boxplot", config=PLOTLY_CONFIG)
                render_population_key(graphed_populations)

            st.write("")
            with st.container(border=True):
                st.write("**Mann-Whitney U results** (Bonferroni-corrected across all pairs and populations tested)")
                render_comparison_stats_table(results, "custom_stats", "custom_comparison_stats.csv")
