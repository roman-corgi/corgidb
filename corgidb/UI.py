import base64
from io import BytesIO
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from IPython.display import display, HTML
import ipywidgets as widgets

from roman_pointing.roman_observability import (
    get_target_coords,
    compute_roman_angles,
    compute_keepout,
)
from Reference_Star_Selection_Tool import (
    load_catalog,
    select_ref_star,
    get_observable_windows,
    build_skycoord,
    SUN_MIN,
    SUN_MAX,
    MAX_PITCH_DIFF,
    REF_GRADES,
    SORT_MODES,
    SORT_MODE_LABELS,
)


SORT_MODE_TO_COLUMN = {
    'valid_days':     'valid_days',
    'closest_mag':    'mag',
    'brightest':      'mag',
    'faintest':       'mag',
    'closest_pitch':  'pitch',
    'farthest_pitch': 'pitch',
}

STYLE_ACTIVE_TH = (
    "padding:7px 10px;border:1px solid #ddd;"
    "background:#e65100;color:white;font-weight:bold"
)
STYLE_NORMAL_TH = (
    "padding:7px 10px;border:1px solid #ddd;"
    "background:#1565c0;color:white"
)
STYLE_ACTIVE_TD = (
    "padding:6px 10px;border:1px solid #ddd;"
    "background:#fff3e0;font-weight:bold"
)
STYLE_NORMAL_TD = "padding:6px 10px;border:1px solid #ddd"

GRADE_COLORS = {"A": "#2e7d32", "B": "#f57c00", "C": "#c62828"}


def fig_to_html(fig):
    """Encode a Matplotlib figure as an inline HTML ``<img>`` tag.

    Args:
        fig (matplotlib.figure.Figure): The figure to encode.

    Returns:
        str: An HTML ``<img>`` tag containing the figure as a base64-encoded PNG.
    """
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return f'<img src="data:image/png;base64,{b64}" style="max-width:100%"/>'


def html_panel(content, title=""):
    """Wrap content in a styled card panel.

    Args:
        content (str): Inner HTML content.
        title (str, optional): Panel heading text. Defaults to no heading.

    Returns:
        str: HTML string for the card panel.
    """
    header = f"<h3 style='margin:0 0 8px 0;color:#333'>{title}</h3>" if title else ""
    return (
        "<div style='background:#f9f9f9;border:1px solid #ddd;border-radius:6px;"
        "padding:16px;margin-bottom:12px;font-family:monospace;font-size:13px'>"
        f"{header}{content}</div>"
    )


def html_badge(text, color):
    """Render a small inline badge with a solid background.

    Args:
        text (str): Badge label (may contain HTML entities).
        color (str): CSS background color string.

    Returns:
        str: HTML ``<span>`` badge.
    """
    return (
        f"<span style='background:{color};color:white;border-radius:4px;"
        f"padding:2px 7px;font-size:12px;font-weight:bold'>{text}</span>"
    )


def grade_color(grade):
    """Return the display color for a grade letter.

    Args:
        grade (str): Grade letter, one of 'A', 'B', or 'C'.

    Returns:
        str: CSS color string.
    """
    return GRADE_COLORS.get(grade, "#555")


def mag_header(effective_sort, sci_mag_available):
    """Return the appropriate magnitude column header label.

    Args:
        effective_sort (str): The resolved sort mode in use.
        sci_mag_available (bool): Whether the science target magnitude is known.

    Returns:
        str: Either ``'Mag Diff'`` or ``'Magnitude'``.
    """
    return "Mag Diff" if (effective_sort == 'closest_mag' and sci_mag_available) else "Magnitude"


def mag_cell(ref, effective_sort, sci_mag_available):
    """Format the magnitude value for a single reference star table cell.

    Args:
        ref (dict): Reference star result dictionary.
        effective_sort (str): The resolved sort mode in use.
        sci_mag_available (bool): Whether the science target magnitude is known.

    Returns:
        str: HTML-formatted magnitude or magnitude difference string.
    """
    if effective_sort == 'closest_mag' and sci_mag_available and ref.get('mag_diff') is not None:
        return f"&Delta;{ref['mag_diff']:.2f}"
    if ref.get('mag') is not None:
        return f"{ref['mag']:.2f}"
    return "&mdash;"


def collect_best_refs(result, ref_coords_map):
    """Collect the best reference star and its coordinates for each window.

    Args:
        result (dict): Return value from :func:`select_ref_star`.
        ref_coords_map (dict): Mapping of star name to
            ``astropy.coordinates.SkyCoord``.

    Returns:
        list of tuple: Each tuple contains:
            - window_index (int)
            - star_name (str)
            - coord (astropy.coordinates.SkyCoord)
            - win_start (astropy.time.Time or None)
            - win_end (astropy.time.Time or None)
    """
    out = []
    for i, win in enumerate(result.get("observable_windows", [])):
        best = win.get("best_ref")
        if best is None:
            continue
        coord = ref_coords_map.get(best["reference_star"])
        if coord is None:
            continue
        out.append((
            i,
            best["reference_star"],
            coord,
            win.get("_win_start"),
            win.get("_win_end"),
        ))
    return out


def plot_solar_angle(ts, sci_solar, result, ref_coords_map):
    """Plot solar angle over time for the science target and best references.

    Args:
        ts (astropy.time.Time): Time array for the full analysis period.
        sci_solar (astropy.units.Quantity): Science target solar angles.
        result (dict): Return value from :func:`select_ref_star`.
        ref_coords_map (dict): Mapping of star name to
            ``astropy.coordinates.SkyCoord``.

    Returns:
        str: HTML ``<img>`` tag containing the rendered figure.
    """
    fig, ax = plt.subplots(figsize=(11, 5))
    dates = [datetime.fromisoformat(t.iso) for t in ts]
    ax.plot(dates, sci_solar.to(u.deg).value, lw=2,
            label="Science target", color="#1565c0")

    colors = plt.cm.tab10.colors
    for idx, (i, name, coord, ws, we) in enumerate(collect_best_refs(result, ref_coords_map)):
        if ws is None or we is None:
            continue
        try:
            dur = we.mjd - ws.mjd
            _, ref_sun, _, _ = compute_roman_angles(coord, ws.isot, dur, time_step=1.0)
            wdates = [
                datetime.fromisoformat((ws + d * u.day).iso)
                for d in np.arange(0, dur, 1.0)
            ]
            ax.plot(wdates, ref_sun.to(u.deg).value, lw=1.5, linestyle="--",
                    color=colors[(idx + 1) % len(colors)],
                    label=f"Ref: {name} (Win {i + 1})")
        except Exception:
            pass

    xlim = ax.get_xlim()
    ax.axhline(SUN_MIN, color="k", linestyle="--", lw=1)
    ax.axhline(SUN_MAX, color="k", linestyle="--", lw=1)
    ax.fill_between(xlim, [SUN_MIN] * 2, [0] * 2,
                    hatch="/", color="none", edgecolor="k", alpha=0.25)
    ax.fill_between(xlim, [SUN_MAX] * 2, [180] * 2,
                    hatch="\\", color="none", edgecolor="k", alpha=0.25)
    ax.set_xlabel("Date")
    ax.set_ylabel("Solar Angle (deg)")
    ax.set_title("Solar Angle vs Time\n(dashed = best reference star per window)")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
    fig.subplots_adjust(right=0.75)
    return fig_to_html(fig)


def plot_pitch_angle(ts, sci_pitch, result, ref_coords_map):
    """Plot pitch angle over time for the science target and best references.

    Args:
        ts (astropy.time.Time): Time array for the full analysis period.
        sci_pitch (astropy.units.Quantity): Science target pitch angles.
        result (dict): Return value from :func:`select_ref_star`.
        ref_coords_map (dict): Mapping of star name to
            ``astropy.coordinates.SkyCoord``.

    Returns:
        str: HTML ``<img>`` tag containing the rendered figure.
    """
    fig, ax = plt.subplots(figsize=(11, 5))
    dates = [datetime.fromisoformat(t.iso) for t in ts]
    ax.plot(dates, sci_pitch.to(u.deg).value, lw=2,
            label="Science target", color="#1565c0")

    colors = plt.cm.tab10.colors
    for idx, (i, name, coord, ws, we) in enumerate(collect_best_refs(result, ref_coords_map)):
        if ws is None or we is None:
            continue
        try:
            dur = we.mjd - ws.mjd
            _, _, _, ref_p = compute_roman_angles(coord, ws.isot, dur, time_step=1.0)
            wdates = [
                datetime.fromisoformat((ws + d * u.day).iso)
                for d in np.arange(0, dur, 1.0)
            ]
            ax.plot(wdates, ref_p.to(u.deg).value, lw=1.5, linestyle="--",
                    color=colors[(idx + 1) % len(colors)],
                    label=f"Ref: {name} (Win {i + 1})")
        except Exception:
            pass

    ax.set_xlabel("Date")
    ax.set_ylabel("Pitch Angle (deg)")
    ax.set_title("Pitch Angle vs Time\n(dashed = best reference star per window)")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
    fig.subplots_adjust(right=0.75)
    return fig_to_html(fig)


def plot_keepout_map(ts, keepout_arr, sci_name):
    """Plot a colour-coded Roman keepout map for the science target.

    Args:
        ts (astropy.time.Time): Time array for the full analysis period.
        keepout_arr (numpy.ndarray): Boolean array; ``True`` = observable.
        sci_name (str): Science target name used in the plot title.

    Returns:
        str: HTML ``<img>`` tag containing the rendered figure.
    """
    dates = [datetime.fromisoformat(t.iso) for t in ts]
    date_nums = mdates.date2num(dates)
    fig, ax = plt.subplots(figsize=(12, 2.5))
    cmap = matplotlib.colors.ListedColormap(["black", "green"])
    extended = np.append(
        date_nums,
        date_nums[-1] + (date_nums[-1] - date_nums[-2])
    )
    ax.pcolormesh(
        extended, [0, 1],
        keepout_arr.reshape(1, -1).astype(int),
        cmap=cmap, shading="flat"
    )
    transition_indices = np.where(np.diff(keepout_arr.astype(int)) != 0)[0] + 1
    tick_indices = np.unique(
        np.concatenate(([0], transition_indices, [len(ts) - 1]))
    )
    ax.set_xticks(mdates.date2num([dates[i] for i in tick_indices]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right", fontsize=9)
    ax.set_yticks([0.5])
    ax.set_yticklabels([sci_name])
    ax.set_ylim(0, 1)
    cbar = plt.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap),
        ticks=[0.25, 0.75], ax=ax
    )
    cbar.ax.set_yticklabels(["Unavailable", "Available"])
    ax.set_title(
        f"Roman Keepout Map — {sci_name}\n"
        f"{ts[0].iso[:10]} to {ts[-1].iso[:10]}"
    )
    fig.tight_layout()
    return fig_to_html(fig)


def plot_pitch_diff(ts, sci_pitch, result, ref_coords_map):
    """Plot pitch angle difference between science target and best references.

    Args:
        ts (astropy.time.Time): Time array for the full analysis period.
        sci_pitch (astropy.units.Quantity): Science target pitch angles.
        result (dict): Return value from :func:`select_ref_star`.
        ref_coords_map (dict): Mapping of star name to
            ``astropy.coordinates.SkyCoord``.

    Returns:
        str: HTML ``<img>`` tag containing the rendered figure.
    """
    fig, ax = plt.subplots(figsize=(11, 5))
    sci_pitch_vals = sci_pitch.to(u.deg).value
    ax.axhline(MAX_PITCH_DIFF, color="red", linestyle="--", lw=1.5,
               label=f"Max allowed ({MAX_PITCH_DIFF} deg)")

    colors = plt.cm.tab10.colors
    for idx, (i, name, coord, ws, we) in enumerate(collect_best_refs(result, ref_coords_map)):
        if ws is None or we is None:
            continue
        try:
            dur = we.mjd - ws.mjd
            _, _, _, ref_p = compute_roman_angles(coord, ws.isot, dur, time_step=1.0)
            ref_pitch_vals = ref_p.to(u.deg).value
            start_idx = int(ws.mjd - ts[0].mjd)
            sci_slice = sci_pitch_vals[start_idx:start_idx + len(ref_pitch_vals)]
            min_len = min(len(sci_slice), len(ref_pitch_vals))
            diff = np.abs(sci_slice[:min_len] - ref_pitch_vals[:min_len])
            wdates = [
                datetime.fromisoformat((ws + d * u.day).iso)
                for d in np.arange(0, min_len)
            ]
            ax.plot(wdates, diff, lw=1.5,
                    color=colors[(idx + 1) % len(colors)],
                    label=f"Ref: {name} (Win {i + 1})")
        except Exception:
            pass

    ax.set_xlabel("Date")
    ax.set_ylabel("Pitch Angle Difference (deg)")
    ax.set_title("Pitch Angle Difference: Science Target vs Best Reference per Window")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
    fig.subplots_adjust(right=0.75)
    return fig_to_html(fig)


def availability_html(result):
    """Build an HTML availability calendar for all observable windows.

    Args:
        result (dict): Return value from :func:`select_ref_star`.

    Returns:
        str: HTML string showing a per-window, per-date availability grid.
    """
    wins = result.get("observable_windows", [])
    if not wins:
        return "<p style='color:#c62828'>No observable windows found.</p>"

    sections = []
    for i, win in enumerate(wins):
        valid_refs = win.get("valid_refs", [])
        if not valid_refs:
            sections.append(
                f"<h4 style='color:#c62828;margin:12px 0 4px'>Window {i + 1}: "
                f"{win['start']} to {win['end']} — No valid reference stars</h4>"
            )
            continue

        sections.append(
            f"<h4 style='margin:16px 0 6px 0;color:#1565c0'>"
            f"Window {i + 1}: {win['start']} to {win['end']} "
            f"({win['duration_days']:.1f} days)</h4>"
        )

        star_names = [r["reference_star"] for r in valid_refs]
        grades = {r["reference_star"]: r.get("grade", "?") for r in valid_refs}
        n_days_map = {r["reference_star"]: r.get("n_valid_days", 0) for r in valid_refs}
        all_dates = sorted({d for r in valid_refs for d in r.get("valid_dates", [])})
        avail_lookup = {
            r["reference_star"]: set(r.get("valid_dates", []))
            for r in valid_refs
        }

        th_stars = "".join(
            "<th style='padding:5px 8px;border:1px solid #ddd;white-space:nowrap;"
            "background:#1565c0;color:white;font-size:11px'>"
            f"{name}<br><span style='font-weight:normal;font-size:10px'>"
            f"Grade {grades[name]} &middot; {n_days_map[name]}d</span></th>"
            for name in star_names
        )
        header_row = (
            "<tr><th style='padding:5px 8px;border:1px solid #ddd;"
            "background:#1565c0;color:white'>Date</th>"
            f"{th_stars}</tr>"
        )

        body_rows = ""
        for date in all_dates:
            cells = "".join(
                (
                    "<td style='padding:4px 8px;border:1px solid #ddd;"
                    "text-align:center;background:#e8f5e9'>&#10003;</td>"
                    if date in avail_lookup.get(name, set()) else
                    "<td style='padding:4px 8px;border:1px solid #ddd;"
                    "text-align:center;color:#ccc'>&mdash;</td>"
                )
                for name in star_names
            )
            body_rows += (
                f"<tr><td style='padding:4px 8px;border:1px solid #ddd;"
                f"font-family:monospace;white-space:nowrap'>{date}</td>{cells}</tr>"
            )

        sections.append(
            "<div style='overflow-x:auto;margin-bottom:20px'>"
            "<table style='border-collapse:collapse;font-size:12px'>"
            f"<thead>{header_row}</thead><tbody>{body_rows}</tbody>"
            "</table></div>"
        )

    return "".join(sections)


def results_html(result, band_label):
    """Build the main reference star results table as an HTML string.

    The column corresponding to the active sort metric receives an orange
    header and orange cell background. No duplicate columns are added.

    Columns: Window, Start, End, Duration, Reference Star, Grade,
    Magnitude/Mag Diff, Valid Days, Min Pitch Diff.

    Args:
        result (dict): Return value from :func:`select_ref_star`.
        band_label (str): ``'V'`` or ``'I'``, used in the summary header.

    Returns:
        str: HTML string for the full results summary and table.
    """
    sci = result["science_target"]
    vis = result["visibility_pct"]
    sort_method = result["sort_method"]
    wins = result["observable_windows"]
    active_grades = result.get("allowed_grades", REF_GRADES)
    sort_mode = result.get("sort_mode", "valid_days")
    effective = result.get("effective_sort", sort_mode)
    sci_mag = result.get("sci_mag")
    sci_mag_available = sci_mag is not None

    mag_hdr = mag_header(effective, sci_mag_available)
    sci_mag_str = f"{sci_mag:.2f}" if sci_mag_available else "N/A"
    active_col = SORT_MODE_TO_COLUMN.get(effective, '')

    def th(label, col_id):
        style = STYLE_ACTIVE_TH if col_id == active_col else STYLE_NORMAL_TH
        suffix = " &#9660;" if col_id == active_col else ""
        return f"<th style='{style}'>{label}{suffix}</th>"

    def td(value, col_id, extra_style=""):
        style = STYLE_ACTIVE_TD if col_id == active_col else STYLE_NORMAL_TD
        if extra_style:
            style += ";" + extra_style
        return f"<td style='{style}'>{value}</td>"

    header = (
        f"<div style='margin-bottom:14px;line-height:1.8'>"
        f"<b>Science target:</b> {sci} &nbsp;|&nbsp;"
        f"<b>Band:</b> {band_label} &nbsp;|&nbsp;"
        f"<b>Contrast:</b> {result.get('contrast', '&mdash;')} &nbsp;|&nbsp;"
        f"<b>Visibility:</b> {vis:.1f}%<br>"
        f"<b>Allowed grades:</b> {', '.join(active_grades)} &nbsp;|&nbsp;"
        f"<b>Sort mode:</b> {sort_mode} &nbsp;&mdash;&nbsp;"
        f"{SORT_MODE_LABELS.get(sort_mode, '')}<br>"
        f"<b>Science {band_label}-mag:</b> {sci_mag_str} &nbsp;|&nbsp;"
        f"<b>Grade column:</b> {result.get('grade_column', '&mdash;')}<br>"
        f"<b>Sort method:</b> {sort_method}"
        f"</div>"
    )

    table_header = (
        "<tr>"
        + f"<th style='{STYLE_NORMAL_TH}'>Window</th>"
        + f"<th style='{STYLE_NORMAL_TH}'>Start</th>"
        + f"<th style='{STYLE_NORMAL_TH}'>End</th>"
        + f"<th style='{STYLE_NORMAL_TH}'>Duration</th>"
        + f"<th style='{STYLE_NORMAL_TH}'>Reference Star</th>"
        + f"<th style='{STYLE_NORMAL_TH}'>Grade</th>"
        + th(mag_hdr, 'mag')
        + th("Valid Days", 'valid_days')
        + th("Min Pitch Diff", 'pitch')
        + "</tr>"
    )

    rows = ""
    for i, win in enumerate(wins):
        valid_refs = win.get("valid_refs", [])
        n_valid = len(valid_refs)

        if not valid_refs:
            rows += (
                f"<tr>"
                f"<td style='{STYLE_NORMAL_TD}'><b>Window {i + 1}</b></td>"
                f"<td style='{STYLE_NORMAL_TD}'>{win['start']}</td>"
                f"<td style='{STYLE_NORMAL_TD}'>{win['end']}</td>"
                f"<td style='{STYLE_NORMAL_TD}'>{win['duration_days']:.1f} d</td>"
                f"<td colspan='5' style='{STYLE_NORMAL_TD};"
                f"color:#c62828;font-style:italic'>"
                f"No suitable reference stars found</td>"
                f"</tr>"
            )
            continue

        for j, ref in enumerate(valid_refs):
            ref_name = ref["reference_star"]
            grade = ref.get("grade", "—")
            n_days = ref.get("n_valid_days", 0)
            min_pdiff = ref.get("min_pitch_diff", None)

            grade_badge = html_badge(grade, grade_color(grade))
            best_badge = html_badge("&#9733; best", "#1565c0") if j == 0 else ""
            mag_str = mag_cell(ref, effective, sci_mag_available)
            pdiff_str = (
                f"{min_pdiff:.4f}&deg;"
                if min_pdiff is not None and min_pdiff < 999
                else "&mdash;"
            )

            if j == 0:
                win_cell = (
                    f"<td rowspan='{n_valid}' style='{STYLE_NORMAL_TD};"
                    f"vertical-align:top'><b>Window {i + 1}</b></td>"
                    f"<td rowspan='{n_valid}' style='{STYLE_NORMAL_TD};"
                    f"vertical-align:top'>{win['start']}</td>"
                    f"<td rowspan='{n_valid}' style='{STYLE_NORMAL_TD};"
                    f"vertical-align:top'>{win['end']}</td>"
                    f"<td rowspan='{n_valid}' style='{STYLE_NORMAL_TD};"
                    f"vertical-align:top'>{win['duration_days']:.1f} d</td>"
                )
            else:
                win_cell = ""

            rows += (
                f"<tr>{win_cell}"
                f"<td style='{STYLE_NORMAL_TD};font-weight:bold;color:#1b5e20'>"
                f"{ref_name} {best_badge}</td>"
                f"<td style='{STYLE_NORMAL_TD};text-align:center'>{grade_badge}</td>"
                + td(mag_str, 'mag', 'text-align:right')
                + td(n_days, 'valid_days', 'text-align:center')
                + td(pdiff_str, 'pitch')
                + "</tr>"
            )

    table = (
        "<div style='overflow-x:auto'>"
        "<table style='border-collapse:collapse;width:100%;font-size:13px'>"
        f"<thead>{table_header}</thead><tbody>{rows}</tbody>"
        "</table></div>"
    )
    return header + table


class ReferenceStarPickerUI:
    """Jupyter ipywidgets UI for the ReferenceStarPicker module.

    Provides a form-based interface for configuring and running
    :func:`select_ref_star`, then displays tabbed output panels for the
    results table, angle plots, keepout map, and availability calendar.

    Args:
        engine (sqlalchemy.engine.base.Engine, optional): Connected plandb
            engine. Required at runtime; may be passed later. Defaults to None.
    """

    def __init__(self, engine=None):
        self.engine = engine
        self.catalog = None

        self.w_target = widgets.Text(
            value="47 Uma",
            description="Science Target:",
            style={"description_width": "140px"},
            layout=widgets.Layout(width="380px"),
            placeholder="SIMBAD-resolvable name",
        )
        self.w_band = widgets.Dropdown(
            options=[
                ("Band 1 — V band (NFB)",  1),
                ("Band 3 — I band (spec)", 3),
                ("Band 4 — I band (wide)", 4),
            ],
            value=1,
            description="Band:",
            style={"description_width": "140px"},
            layout=widgets.Layout(width="320px"),
        )
        self.w_contrast = widgets.Dropdown(
            options=[("High contrast", "high"), ("Medium contrast", "med")],
            value="high",
            description="Contrast:",
            style={"description_width": "140px"},
            layout=widgets.Layout(width="280px"),
        )
        self.w_start = widgets.Text(
            value="2027-01-01T00:00:00",
            description="Analysis Start:",
            style={"description_width": "140px"},
            layout=widgets.Layout(width="380px"),
            placeholder="YYYY-MM-DDTHH:MM:SS",
        )
        self.w_days = widgets.BoundedFloatText(
            value=365, min=1, max=3650,
            description="Duration (days):",
            style={"description_width": "140px"},
            layout=widgets.Layout(width="280px"),
        )
        self.w_timestep = widgets.BoundedFloatText(
            value=1.0, min=0.1, max=30.0, step=0.1,
            description="Time step (days):",
            style={"description_width": "140px"},
            layout=widgets.Layout(width="280px"),
        )
        self.w_grades = widgets.SelectMultiple(
            options=[
                ("A — best PSF quality",       "A"),
                ("B — acceptable PSF quality", "B"),
                ("C — marginal PSF quality",   "C"),
            ],
            value=("A", "B", "C"),
            description="Allowed grades:",
            style={"description_width": "140px"},
            layout=widgets.Layout(width="340px", height="82px"),
            rows=3,
        )
        self.w_sort_mode = widgets.Dropdown(
            options=[
                ("Most valid days first (default)",       "valid_days"),
                ("Closest magnitude to science target",   "closest_mag"),
                ("Brightest first (ascending magnitude)", "brightest"),
                ("Faintest first (descending magnitude)", "faintest"),
                ("Smallest pitch angle difference first", "closest_pitch"),
                ("Largest pitch angle difference first",  "farthest_pitch"),
            ],
            value="valid_days",
            description="Sort by:",
            style={"description_width": "140px"},
            layout=widgets.Layout(width="440px"),
        )
        self.w_run = widgets.Button(
            description="Find Reference Stars",
            button_style="primary",
            icon="search",
            layout=widgets.Layout(width="220px", height="36px"),
        )
        self.w_status = widgets.HTML(value="")

        self.out_results = widgets.Output()
        self.out_solar = widgets.Output()
        self.out_pitch = widgets.Output()
        self.out_pitchdiff = widgets.Output()
        self.out_keepout = widgets.Output()
        self.out_avail = widgets.Output()

        self.tabs = widgets.Tab(children=[
            self.out_results, self.out_solar, self.out_pitch,
            self.out_pitchdiff, self.out_keepout, self.out_avail,
        ])
        for i, title in enumerate([
            "Results", "Solar Angle", "Pitch Angle",
            "Pitch Difference", "Keepout Map", "Availability Calendar",
        ]):
            self.tabs.set_title(i, title)

        self.w_run.on_click(self.on_run)

    def display(self):
        """Render the full UI widget in the current Jupyter cell output."""
        title = widgets.HTML(
            "<h2 style='color:#1565c0;margin-bottom:4px'>"
            "&#128301; Reference Star Picker</h2>"
            "<p style='color:#555;margin-top:0'>"
            "Roman Space Telescope — Coronagraph Reference Star Selection</p>"
        )
        grades_hint = widgets.HTML(
            "<span style='font-size:11px;color:#777'>"
            "Ctrl/Cmd+click to select multiple grades. "
            "Order defines ranking priority."
            "</span>"
        )
        controls = widgets.VBox(
            [
                widgets.HBox([self.w_target,   self.w_band]),
                widgets.HBox([self.w_contrast, self.w_days]),
                widgets.HBox([self.w_start,    self.w_timestep]),
                widgets.HBox([self.w_sort_mode]),
                widgets.HBox(
                    [self.w_grades, grades_hint],
                    layout=widgets.Layout(align_items="center", gap="10px"),
                ),
                widgets.HBox(
                    [self.w_run, self.w_status],
                    layout=widgets.Layout(align_items="center", gap="16px"),
                ),
            ],
            layout=widgets.Layout(
                padding="12px", border="1px solid #ddd", margin_bottom="12px"
            ),
        )
        display(widgets.VBox([title, controls, self.tabs]))

    def set_status(self, msg, color="#555"):
        """Update the status label beneath the run button.

        Args:
            msg (str): Status message (may contain HTML).
            color (str, optional): CSS color for the message text.
                Defaults to ``'#555'``.
        """
        self.w_status.value = (
            f"<span style='color:{color};font-size:13px'>{msg}</span>"
        )

    def on_run(self, btn):
        """Handle a click on the 'Find Reference Stars' button.

        Loads the catalog if not already cached, runs the full reference star
        search, and populates all output tabs with results.

        Args:
            btn: The button widget that triggered this callback (unused).
        """
        self.w_run.disabled = True
        self.set_status("&#9203; Loading catalog...", "#f57c00")
        for out in [
            self.out_results, self.out_solar, self.out_pitch,
            self.out_pitchdiff, self.out_keepout, self.out_avail,
        ]:
            out.clear_output()

        try:
            if self.catalog is None:
                if self.engine is None:
                    self.set_status(
                        "&#10060; No database engine provided. "
                        "Pass engine= to ReferenceStarPickerUI().",
                        "#c62828",
                    )
                    return
                self.catalog = load_catalog(self.engine)

            sci_name = self.w_target.value.strip()
            band = self.w_band.value
            contrast = self.w_contrast.value
            start = self.w_start.value.strip()
            days = float(self.w_days.value)
            time_step = float(self.w_timestep.value)
            allowed_grades = list(self.w_grades.value)
            sort_mode = self.w_sort_mode.value
            band_label = "V" if band == 1 else "I"

            if not allowed_grades:
                self.set_status(
                    "&#10060; Select at least one grade (A, B, or C).", "#c62828"
                )
                return

            self.set_status(
                f"&#9203; Querying SIMBAD for '{sci_name}'...", "#f57c00"
            )
            coords = get_target_coords([sci_name])
            if sci_name not in coords:
                self.set_status(
                    f"&#10060; '{sci_name}' not found in SIMBAD.", "#c62828"
                )
                return
            sci_coord = coords[sci_name]

            self.set_status("&#9203; Computing Roman visibility...", "#f57c00")
            ts, keepout, solar_angles = compute_keepout(
                {sci_name: sci_coord}, start, days, time_step
            )
            sci_keepout = keepout[sci_name]
            sci_solar = solar_angles[sci_name]
            _, _, _, sci_pitch = compute_roman_angles(
                sci_coord, start, days, time_step
            )

            self.set_status("&#9203; Searching for reference stars...", "#f57c00")
            result = select_ref_star(
                sci_name, start, days,
                band=band, contrast=contrast,
                catalog=self.catalog, engine=self.engine,
                time_step=time_step,
                allowed_grades=allowed_grades,
                sort_mode=sort_mode,
            )

            if "error" in result:
                self.set_status(f"&#10060; {result['error']}", "#c62828")
                return

            wins_raw = get_observable_windows(ts, sci_keepout)
            for i, win in enumerate(result["observable_windows"]):
                win["_win_start"] = wins_raw[i][0] if i < len(wins_raw) else None
                win["_win_end"] = wins_raw[i][1] if i < len(wins_raw) else None

            ref_coords_map = {}
            for _, ref_row in self.catalog.iterrows():
                name = ref_row.get("main_id", "")
                if isinstance(name, str) and name.strip():
                    try:
                        ref_coords_map[name] = build_skycoord(ref_row)
                    except Exception:
                        pass

            self.set_status("&#9203; Rendering results...", "#f57c00")

            with self.out_results:
                display(HTML(html_panel(
                    results_html(result, band_label),
                    title="Reference Star Selection Results",
                )))

            with self.out_solar:
                display(HTML(plot_solar_angle(ts, sci_solar, result, ref_coords_map)))

            with self.out_pitch:
                display(HTML(plot_pitch_angle(ts, sci_pitch, result, ref_coords_map)))

            with self.out_pitchdiff:
                display(HTML(plot_pitch_diff(ts, sci_pitch, result, ref_coords_map)))

            with self.out_keepout:
                display(HTML(plot_keepout_map(ts, sci_keepout, sci_name)))

            with self.out_avail:
                display(HTML(html_panel(
                    availability_html(result),
                    title="Reference Star Availability — Valid Observation Dates per Window",
                )))

            vis = result["visibility_pct"]
            n_wins = len(result["observable_windows"])
            sci_mag = result.get("sci_mag")
            mag_note = (
                f" | sci {band_label}-mag: {sci_mag:.2f}"
                if sci_mag is not None else ""
            )
            self.set_status(
                f"&#10003; Done — {n_wins} window(s), {vis:.1f}% observable"
                f" | grades: {', '.join(allowed_grades)}"
                f" | sort: {sort_mode}{mag_note}",
                "#2e7d32",
            )

        except Exception as exc:
            self.set_status(f"&#10060; Error: {exc}", "#c62828")
            import traceback
            traceback.print_exc()
        finally:
            self.w_run.disabled = False


def launch(engine=None):
    """Instantiate and display the ReferenceStarPickerUI.

    Args:
        engine (sqlalchemy.engine.base.Engine, optional): Connected plandb
            engine. Defaults to None.

    Returns:
        ReferenceStarPickerUI: The instantiated UI object.
    """
    ui = ReferenceStarPickerUI(engine=engine)
    ui.display()
    return ui