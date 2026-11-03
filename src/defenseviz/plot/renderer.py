"""Physical-unit atlas renderer over evaluated DefenseViz result dictionaries.

Only presentation constants (spacing, fonts and a categorical palette) live here.
Scientific values and policy parameters always come from the result argument.
"""
from __future__ import annotations

from itertools import product
import json
import math
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager, ft2font
from matplotlib.colors import to_rgb, to_hex
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch
from matplotlib.textpath import TextToPath

INK = "#24393B"
TEAL = "#287C78"
GRID = "#DFE6E5"
PINK = "#FBE0EB"
PINK_LABEL = "#F3BCD2"
CATEGORY_COLORS = {
    "R-M": "#287C78", "Modification": "#B29BCB", "Nuclease": "#D9878C",
    "CRISPR-Cas": "#A8BF72", "Abortive infection": "#7FAED0",
    "Other mechanisms": "#C6AD75", "Unresolved": "#AAB5BC",
}
CALLER_COLORS = ("#D56B25", "#347DA7", "#7953A2", "#287C78", "#AD5878", "#7E8544", "#536A9C", "#A76D39")
MM_PER_INCH = 25.4
BAR_HEIGHT = 1.65
MAX_CALLERS = 8
MAX_COMPONENTS = 16
KNOWN_CALLER_COLORS = {"padloc": "#D56B25", "defensefinder": "#347DA7", "defensepredictor": "#7953A2"}


def _caller_color(caller, all_callers):
    return KNOWN_CALLER_COLORS.get(caller.lower(), CALLER_COLORS[all_callers.index(caller) % len(CALLER_COLORS)])


def _component_order(value):
    typ, role = value
    types = {"type i": 0, "type ii": 1, "type iii": 2, "type iv": 3}
    roles = {"R": 0, "M": 1, "S": 2, "RM": 3, "C": 4, "Aux": 5}
    return (types.get(str(typ).lower(), 4), str(typ), roles.get(role, 6), role)



def _int(value, field):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{field} must be a nonnegative integer")
    return value


def _chunks(items, size):
    return [items[i:i + size] for i in range(0, len(items), size)] or [[]]


def _wrap(value, width):
    """Preserve every character; do not silently truncate identifiers/labels."""
    lines = []
    for paragraph in str(value).splitlines() or [""]:
        lines.extend(textwrap.wrap(paragraph, width=width, break_long_words=True,
                                   break_on_hyphens=False, replace_whitespace=False) or [""])
    return "\n".join(lines)


class _Font:
    """Choose an installed font that covers labels, rather than export tofu."""
    def __init__(self):
        self.cache = {}
        preferred = ["DejaVu Sans", "Noto Sans", "Noto Sans CJK KR", "Noto Sans CJK SC", "Arial Unicode MS"]
        fonts = font_manager.fontManager.ttflist
        self.fonts = sorted(fonts, key=lambda f: (preferred.index(f.name) if f.name in preferred else 100,
                                                 "sans" not in f.name.lower(), f.name, f.fname))
        self.cmaps = {}
        self.used = set()

    def pick(self, value):
        needed = frozenset(ord(c) for c in str(value) if not c.isspace())
        if needed in self.cache:
            return self.cache[needed]
        for font in self.fonts:
            if font.fname not in self.cmaps:
                try:
                    self.cmaps[font.fname] = set(ft2font.FT2Font(font.fname).get_charmap())
                except (OSError, RuntimeError):
                    self.cmaps[font.fname] = set()
            if needed <= self.cmaps[font.fname]:
                self.used.add(font.name)
                self.cache[needed] = font.name
                return font.name
        bad = " ".join(f"U+{c:04X}" for c in sorted(needed) if not any(c in cm for cm in self.cmaps.values()))
        raise ValueError("No installed font covers a figure label. Install Noto Sans/CJK or a Unicode font. " + bad)


class _Canvas:
    def __init__(self, width, height, font):
        self.width, self.height, self.font = width, height, font
        self.fig = Figure(figsize=(width / MM_PER_INCH, height / MM_PER_INCH), dpi=120, facecolor="white")
        FigureCanvasAgg(self.fig)
        self.ax = self.fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(height, 0))
        self.ax.axis("off")
        self.texts = []

    def text(self, x, y, value, size=7.2, bold=False, color=INK, ha="left", va="center", rotation=0, role="label"):
        value = str(value)
        t = self.ax.text(x, y, value, fontsize=size, fontfamily=self.font.pick(value),
                         fontweight="bold" if bold else "normal", color=color, ha=ha, va=va,
                         rotation=rotation, rotation_mode="anchor", linespacing=1.18, clip_on=False, parse_math=False)
        self.texts.append((t, role, value))
        return t

    def rect(self, x, y, width, height, fill, edge=None, linewidth=.3):
        self.ax.add_patch(Rectangle((x, y), width, height, facecolor=fill,
                                   edgecolor=edge if edge else "none", linewidth=linewidth))

    def line(self, x1, y1, x2, y2, color=GRID, linewidth=.35):
        self.ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth, solid_capstyle="butt")

    def circle(self, x, y, radius, color, *, outline=False):
        self.ax.add_patch(Circle((x, y), radius, facecolor="none" if outline else color,
                                edgecolor=color if outline else "none", linewidth=.8))

    def measure(self):
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        result = []
        for text, role, label in self.texts:
            bb = text.get_window_extent(renderer)
            x1, y1 = self.ax.transData.inverted().transform((bb.x0, bb.y1))
            x2, y2 = self.ax.transData.inverted().transform((bb.x1, bb.y0))
            result.append({"role": role, "text": label, "bounds_mm": [float(x1), float(y1), float(x2), float(y2)],
                           "font_pt": text.get_fontsize()})
        return result


def _validate(result):
    if not isinstance(result, dict):
        raise ValueError("render result must be a mapping returned by integrate")
    samples = result.get("samples", [])
    ids = [s.get("sample_id") for s in samples]
    if any(not isinstance(s, str) or not s for s in ids) or len(set(ids)) != len(ids):
        raise ValueError("Figure samples require unique nonempty sample_id values")
    callers = result.get("callers", [])
    if any(not isinstance(c, str) or not c for c in callers) or len(set(callers)) != len(callers):
        raise ValueError("Figure callers require unique nonempty names")
    families = result.get("families", [])
    fkeys = [(f.get("family"), f.get("category", "Unresolved")) for f in families]
    if len(set(fkeys)) != len(fkeys) or any(not isinstance(f, str) or not f for f, _ in fkeys):
        raise ValueError("Figure family/category columns must be unique and named")
    for s in samples:
        if "highlight" in s and not isinstance(s["highlight"], bool):
            raise ValueError("Sample highlight must be a boolean")
    for key in ("cells", "sample_summary", "rm_components"):
        seen = set()
        for row in result.get(key, []):
            if row.get("sample_id") not in ids:
                raise ValueError(f"{key} references an unknown sample")
            if key == "cells":
                rk = (row["sample_id"], row["family"], row.get("category", "Unresolved"))
                if rk[1:] not in fkeys:
                    raise ValueError("Cell references an unknown family/category column")
                _int(row["cds_count"], "cds_count")
            elif key == "rm_components":
                rk = (row["sample_id"], row["rm_type"], row["role"])
                _int(row["count"], "component count")
                _int(row.get("supported_count", 0), "supported_count")
                if row.get("supported_count", 0) > row["count"]:
                    raise ValueError("supported_count cannot exceed component count")
            else:
                rk = row["sample_id"]
                for field in ("union_cds", "rm_pool_cds", "retained_cds", "unassessed_cds"):
                    _int(row.get(field, 0), field)
                for value in row.get("caller_counts", {}).values():
                    if value is not None:
                        _int(value, "caller_count")
            if rk in seen:
                raise ValueError(f"Duplicate {key} entry")
            seen.add(rk)
    return samples, callers, families


def _panel_a(c, result, x, y, width):
    """Workflow reflects supplied evidence/policy, never promises a fresh run."""
    c.text(x, y, "a   Saved evidence to comparative annotation", 10.8, True, role="panel_title")
    decisions = result.get("decisions", [])
    scored = [r for r in decisions if r.get("status") == "scored"]
    unassessed = sum(r.get("status") == "unassessed" for r in decisions)
    retained = sum(bool(r.get("decision", {}).get("keep")) for r in scored)
    sources = result.get("callers", [])
    policy = result.get("policy", {})
    if "parameters" in policy:
        policy = policy["parameters"]
    support = lambda name: sum(bool(r.get("decision", {}).get(name)) for r in scored)
    cutoff = policy.get("retain_cutoff")
    cutoff_text = str(cutoff) if cutoff is not None else "not supplied"
    max_axes = policy.get("min_capped_evidence_axes")
    cards = [
        ("1  Preserve input scope", TEAL,
         [f"{len(result.get('samples', [])):,} samples; {len(sources):,} caller(s)",
          "Native members / raw hits kept separate", "Missing caller output stays unassessed"]),
        ("2  Keep evidence separate", "#9474AE",
         [f"Profile support: {support('profile_support'):,} scored CDS",
          f"Motif/site support: {support('motif_support'):,}; context: {support('context_support'):,}",
          f"Structure support: {support('strict_structure_support'):,} saved records"]),
        ("3  DNMB policy (custom)" if result.get("policy_metadata", {}).get("is_selected_preset") is False else "3  Apply the DNMB policy", "#B96E7A",
         [f"{len(scored):,} CDS scored; threshold: {cutoff_text}",
          f"Capped evidence axes required: {max_axes if max_axes is not None else 'not supplied'}",
          "Sequence cap + motif + context + structure", "Subtract recorded penalties"]),
        ("4  Report annotation", "#668B92",
         [f"{retained:,} retained; {unassessed:,} unassessed CDS",
          "Gene counts and native systems stay separate", "Activity is unassessed"]),
    ]
    gap = 4.5
    cw = (width - 3 * gap) / 4
    card_height = 28
    top = y + 4.5
    card_bounds = []
    for i, (heading, color, lines) in enumerate(cards):
        left = x + i * (cw + gap)
        c.rect(left, top, cw, card_height, "#FAFBFB", color, .6)
        c.rect(left, top, 1.2, card_height, color)
        # The widest card may be compact on small datasets; preserve text by wrapping.
        chars = max(23, int((cw - 5) / 1.35))
        c.text(left + 3, top + 4, _wrap(heading, chars), 8.0, True, color=color, role="workflow_header")
        wrapped = "\n".join(_wrap(line, chars) for line in lines)
        c.text(left + 3, top + 9, wrapped, 6.6, va="top", role="workflow_body")
        if i < 3:
            c.ax.add_patch(FancyArrowPatch((left + cw + .3, top + card_height / 2),
                                           (left + cw + gap - .3, top + card_height / 2),
                                           arrowstyle="-|>", mutation_scale=5, linewidth=.6, color="#8AA29F"))
        card_bounds.append([left, top, left + cw, top + card_height])
    return {"bounds_mm": [x, y - 3, x + width, top + card_height], "cards_mm": card_bounds}


def _text_width(value, size, font, bold=False):
    prop = font_manager.FontProperties(family=font.pick(value), size=size,
                                      weight="bold" if bold else "normal")
    return TextToPath().get_text_width_height_descent(str(value), prop, False)[0] * 25.4 / 72


def _wrap_mm(value, maximum, size, font, bold=False):
    lines = []
    for paragraph in str(value).splitlines() or [""]:
        line = ""
        for char in paragraph:
            candidate = line + char
            if line and _text_width(candidate, size, font, bold) > maximum:
                space = line.rfind(" ")
                if space > 0:
                    lines.append(line[:space]); line = line[space + 1:] + char
                else:
                    lines.append(line); line = char
            else:
                line = candidate
        lines.append(line)
    return "\n".join(lines)


def _label(sample, font):
    label = str(sample.get("label") or sample["sample_id"])
    species = str(sample.get("species") or "")
    if species and species not in label:
        label += "\n" + species
    return _wrap_mm(label, 59., 7.3, font, sample.get("highlight", False))


def _rows(samples, maximum, font):
    """Bound physical height as well as number of rows; long names stay readable."""
    pages, current, height = [], [], 0.
    for sample in samples:
        rh = max(5.4, len(_label(sample, font).splitlines()) * 3.05 + 1.7)
        if current and (len(current) >= maximum or height + rh > 218):
            pages.append(current); current, height = [], 0.
        current.append((sample, rh)); height += rh
    if current:
        pages.append(current)
    return pages or [[]]


def _column_pages(callers, components, families, max_families):
    call_pages, rm_pages, family_pages = _chunks(callers, MAX_CALLERS), _chunks(components, MAX_COMPONENTS), _chunks(families, max_families)
    pages = []
    for i in range(max(len(call_pages), len(rm_pages), len(family_pages))):
        # Caller and component columns repeat when they fit on one page. When there
        # are many, their complete disjoint slices are named in the layout ledger.
        cp = call_pages[0] if len(call_pages) == 1 else call_pages[i] if i < len(call_pages) else []
        rp = rm_pages[0] if len(rm_pages) == 1 else rm_pages[i] if i < len(rm_pages) else []
        fp = family_pages[i] if i < len(family_pages) else []
        pages.append((cp, rp, fp))
    return pages


def _make_page(result, row_page, col_page, font, page_index, pages, title, scales):
    samples, row_heights = [r[0] for r in row_page], [r[1] for r in row_page]
    callers, components, families = col_page
    label_width, family_width = 64., 7.
    caller_width = max(15. + max(0, len(f"{scales['caller']:,}") - 3) * 1.5,
                       min(25., max([_text_width(name, 6.6, font, True) + 2 for name in result.get('callers', [])] or [15.])))
    summary_width = 13. + max(0, len(f"{scales['summary']:,}") - 3) * 1.5
    component_width = max(6.8, max([len(str(r['count'])) for r in result.get('rm_components', [])] or [1]) * 1.6 + 2)
    left, gap = 7., 4.
    cursor = left + label_width
    caller_start = cursor
    cursor += len(callers) * caller_width + (gap if callers else 0)
    summary_start = cursor
    cursor += 4 * summary_width + gap
    component_start = cursor
    cursor += len(components) * component_width + (gap if components else 0)
    family_start = cursor
    family_widths = [max(family_width, len(_wrap(f["family"], 26).splitlines()) * 2.8 + 2,
                         max([len(str(r["cds_count"])) for r in result.get("cells", []) if r["family"] == f["family"] and r.get("category", "Unresolved") == f.get("category", "Unresolved")] or [1]) * 1.5 + 2) for f in families]
    cursor += sum(family_widths)
    width = max(270., cursor + 7.)
    title_label = _wrap_mm(title, width - 42., 12.5, font, True)
    title_extra = (len(title_label.splitlines()) - 1) * 5.2
    caller_labels = {name: _wrap_mm(name, caller_width - 1.4, 6.6, font, True) for name in callers}
    header_height = max(6., max([len(name.splitlines()) for name in caller_labels.values()] or [1]) * 2.8 + 1.5)
    body_top = 81. + title_extra + header_height
    header_y = 81. + title_extra + (header_height - 2) / 2
    body_height = sum(row_heights) or 12.
    max_label_chars = max([len(line) for f in families for line in _wrap(f["family"], 26).splitlines()] +
                          [len(line) for typ, _ in components for line in _wrap(str(typ), 25).splitlines()] or [0])
    label_depth = max(8., max_label_chars * 1.5 + 5.)
    height = body_top + body_height + label_depth + 25.
    c = _Canvas(width, height, font)
    c.text(left, 7 + title_extra / 2, title_label, 12.5, True, role="title")
    c.text(width - 7, 7, f"{page_index} / {pages}", 7.5, ha="right", role="page_number")
    c.text(left, 13 + title_extra, "Saved-call evaluation and rendering | evidence support is not biological activity", 7.4, role="subtitle")
    panel_a = _panel_a(c, result, left, 22 + title_extra, width - 14)
    c.text(left, 64 + title_extra, "b   Comparative defense landscape", 10.8, True, role="panel_title")
    c.text(left, 70 + title_extra, f"{len(samples):,} sample rows on this page; all counts are distinct CDSs unless stated.", 7.0, role="caption")
    c.text(left, 80 + title_extra, "Samples", 8.2, True, role="table_header")
    if callers:
        c.text(caller_start + caller_width * len(callers) / 2, 77 + title_extra, "Native caller CDSs" if len(callers) * caller_width >= 35 else "Native\nCDSs", 7.5, True, ha="center", role="table_header")
    c.text(summary_start + 2 * summary_width, 77 + title_extra, "Integrated / R-M evidence", 7.5, True, ha="center", role="table_header")
    summary_keys = [("union_cds", "Union"), ("rm_pool_cds", "RM pool"), ("retained_cds", "Retained"), ("unassessed_cds", "Unass.")]
    for j, caller in enumerate(callers):
        c.text(caller_start + (j + .5) * caller_width, header_y, caller_labels[caller], 6.6, True,
               color=_caller_color(caller, result.get("callers", [])), ha="center", role="caller_header")
    for j, (_, label) in enumerate(summary_keys):
        c.text(summary_start + (j + .5) * summary_width, header_y, label, 6.8, True, ha="center", role="summary_header")
    if components:
        c.text(component_start + len(components) * component_width / 2, 77 + title_extra, "R-M candidate components" if len(components) * component_width >= 40 else "R-M\nCDSs", 7.5, True, ha="center", role="table_header")
    if families:
        c.text(family_start + sum(family_widths) / 2, 77 + title_extra, "Native caller CDSs by family" if sum(family_widths) >= 45 else "Family\nCDSs", 7.5, True, ha="center", role="table_header")
    family_centers, cc = [], family_start
    for family, fw in zip(families, family_widths):
        color = CATEGORY_COLORS.get(family.get("category"), CATEGORY_COLORS["Unresolved"])
        c.rect(cc, 81 + title_extra, fw, header_height - 2, color)
        family_centers.append(cc + fw / 2); cc += fw
    for j, (rm_type, role) in enumerate(components):
        c.rect(component_start + j * component_width, 81 + title_extra, component_width, header_height - 2, "#D6E8E4")
        c.text(component_start + (j + .5) * component_width, header_y, role, 6.4, True, ha="center", role="rm_role_header")
    summary_map = {r["sample_id"]: r for r in result.get("sample_summary", [])}
    cell_map = {(r["sample_id"], r["family"], r.get("category", "Unresolved")): r for r in result.get("cells", [])}
    comp_map = {(r["sample_id"], r["rm_type"], r["role"]): r for r in result.get("rm_components", [])}
    coverage = {(r["sample_id"], r["caller"]): r["status"] for r in result.get("coverage", [])}
    rows, bars, matrix, rm_cells = [], [], [], []
    y = body_top
    table_right = max(cursor, summary_start + 4 * summary_width)
    for sample, rh in row_page:
        sid = sample["sample_id"]; yc = y + rh / 2
        c.rect(left, y, table_right - left, rh, PINK if sample.get("highlight", False) else "#FFFFFF")
        if sample.get("highlight", False):
            c.rect(left, y, label_width - 1.4, rh, PINK_LABEL)
        c.text(left + 1.5, yc, _label(sample, font), 7.3, sample.get("highlight", False), role="sample_label")
        summary = summary_map.get(sid, {})
        for j, caller in enumerate(callers):
            x = caller_start + j * caller_width
            value = summary.get("caller_counts", {}).get(caller)
            state = coverage.get((sid, caller))
            if state in {"missing", "not_run"}:
                value = None
            color = _caller_color(caller, result.get("callers", []))
            bars.append(_bar(c, x, yc, caller_width, value, scales["caller"], color, sid, caller))
        for j, (key, _) in enumerate(summary_keys):
            x = summary_start + j * summary_width
            value = summary.get(key)
            color = "#B37591" if key == "unassessed_cds" else TEAL
            bars.append(_bar(c, x, yc, summary_width, value, scales["summary"], color, sid, key))
        for j, (rm_type, role) in enumerate(components):
            x = component_start + (j + .5) * component_width
            r = comp_map.get((sid, rm_type, role), {})
            n, supported = r.get("count", 0), r.get("supported_count", 0)
            c.line(x - component_width / 2, y, x - component_width / 2, y + rh)
            if n:
                c.text(x, yc, n, 7.1, ha="center", role="component_count")
                if supported:
                    c.circle(x, yc, 1.65, TEAL, outline=True)
            else:
                c.text(x, yc, "0", 6.5, color="#9AA5A7", ha="center", role="component_count")
            rm_cells.append({"sample_id": sid, "rm_type": rm_type, "role": role, "count": n,
                             "supported_count": supported, "bounds_mm": [x - component_width / 2, y, x + component_width / 2, y + rh]})
        for f, x, fw in zip(families, family_centers, family_widths):
            cell = cell_map.get((sid, f["family"], f.get("category", "Unresolved")), {})
            count = cell.get("cds_count", 0)
            state = cell.get("count_status")
            if state is None and not any(coverage.get((sid, name)) == "complete" for name in result.get("callers", [])):
                state = "unassessed"
            c.line(x - fw / 2, y, x - fw / 2, y + rh)
            if state == "unassessed":
                c.text(x, yc, "NA", 6.0, color="#7A8587", ha="center", role="family_missing")
            elif count:
                radius = min(1.60, .64 + .30 * math.log2(count + 1))
                color = CATEGORY_COLORS.get(f.get("category"), CATEGORY_COLORS["Unresolved"])
                pale = to_hex(tuple(.65 + .35 * v for v in to_rgb(color)))
                c.circle(x, yc, radius, pale)
                c.text(x, yc, count, 6.2, color=INK, ha="center", role="family_count")
            else:
                c.circle(x, yc, .20, "#DDE4E4")
            matrix.append({"sample_id": sid, "family": f["family"], "category": f.get("category", "Unresolved"),
                           "cds_count": count, "count_status": state or "observed",
                           "bounds_mm": [x - fw / 2, y, x + fw / 2, y + rh]})
        c.line(left, y + rh, table_right, y + rh)
        rows.append({"sample_id": sid, "highlight": sample.get("highlight", False), "bounds_mm": [left, y, table_right, y + rh]})
        y += rh
    if not row_page:
        c.text(left + 2, y + 6, "No samples supplied.", 8.5, color="#677C7E", role="empty_state")
    for j, (rm_type, role) in enumerate(components):
        c.text(component_start + (j + .5) * component_width, y + 2, _wrap(str(rm_type), 25), 6.4,
               ha="right", va="center", rotation=90, role="rm_type_label")
    for f, x in zip(families, family_centers):
        c.text(x, y + 2, _wrap(f["family"], 26), 6.7, ha="right", va="center", rotation=90, role="family_label")
    foot = body_top + body_height + label_depth
    c.text(left, foot, "Counts: observed distinct CDSs. Family rows may overlap; do not sum them as a genome-wide union.", 7.0, role="footnote")
    c.text(left, foot + 4.3, "NA: caller output missing/not run. Small gray dot: zero saved native CDS calls. Ring: one or more supported R-M components.", 7.0, role="footnote")
    c.text(left, foot + 8.6, "RM pool: scored/evaluable candidates plus native/reference R-M evidence. Retained: DNMB evidence threshold passed. Activity unassessed.", 7.0, role="footnote")
    categories = list(dict.fromkeys(f.get("category", "Unresolved") for f in result.get("families", [])))
    lx = left
    for category in categories:
        c.rect(lx, foot + 13.3, 2.2, 2.2, CATEGORY_COLORS.get(category, CATEGORY_COLORS["Unresolved"]))
        c.text(lx + 3.5, foot + 14.4, category, 6.8, role="category_legend")
        lx += 4.6 + len(category) * 1.35
    c.rect(lx, foot + 13.3, 2.2, 2.2, PINK_LABEL)
    c.text(lx + 3.5, foot + 14.4, "Highlighted sample", 6.8, role="highlight_legend")
    text_bounds = c.measure()
    overflow = [r for r in text_bounds if r["bounds_mm"][0] < -.2 or r["bounds_mm"][1] < -.2 or r["bounds_mm"][2] > width + .2 or r["bounds_mm"][3] > height + .2]
    if overflow:
        raise ValueError("Figure text exceeds its page: " + repr([(r["role"], r["text"]) for r in overflow[:4]]))
    layout = {"page": page_index, "size_mm": [width, height], "panel_a": panel_a,
              "panel_b": {"bounds_mm": [left, 60 + title_extra, table_right, foot + 18], "rows": rows,
                          "callers": callers, "families": families, "rm_components": [list(v) for v in components],
                          "bars": bars, "family_cells": matrix, "component_cells": rm_cells}, "text": text_bounds}
    return c, layout


def _bar(canvas, left, y, width, value, maximum, color, sample_id, field):
    # Identical physical thickness and padding for every caller/summary bar.
    reserve = max(6.3, len(f"{maximum:,}") * 1.5 + 1.5)
    track_width = width - reserve
    canvas.rect(left + .6, y - BAR_HEIGHT / 2, track_width, BAR_HEIGHT, "#ECF1F0")
    filled = 0
    if value is None:
        canvas.text(left + width - .6, y, "NA", 6.1, color="#7A8587", ha="right", role="bar_missing")
    else:
        filled = track_width * value / maximum if maximum else 0
        if filled:
            canvas.rect(left + .6, y - BAR_HEIGHT / 2, filled, BAR_HEIGHT, color)
        canvas.text(left + width - .6, y, f"{value:,}", 6.6, ha="right", role="bar_count")
    return {"sample_id": sample_id, "field": field, "value": value, "scale_max": maximum,
            "track_bounds_mm": [left + .6, y - BAR_HEIGHT / 2, left + .6 + track_width, y + BAR_HEIGHT / 2],
            "filled_width_mm": filled, "bar_height_mm": BAR_HEIGHT}


def render(result: dict, output_dir: Path, *, title="DefenseViz", max_rows_per_page=40,
           max_families_per_page=24) -> list[Path]:
    """Export dynamic vector panels a/b, SVG live text and measured layout JSON.

    Each row/family page slice repeats panel a and the relevant count columns.
    A record in ``coverage`` must mark a caller complete before a missing family
    call can be shown as zero. Input counts are never recomputed by the renderer.
    The caller owns atomic publication of ``output_dir``; existing files are
    refused to prevent accidental overwrite of previously reviewed figures.
    """
    for value, name in [(max_rows_per_page, "max_rows_per_page"), (max_families_per_page, "max_families_per_page")]:
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if not isinstance(title, str) or not title.strip():
        raise ValueError("title must be a nonempty string")
    samples, callers, families = _validate(result)
    output_dir = Path(output_dir)
    components = sorted(set((r["rm_type"], r["role"]) for r in result.get("rm_components", [])), key=_component_order)
    font = _Font()
    row_pages = _rows(samples, max_rows_per_page, font)
    col_pages = _column_pages(callers, components, families, max_families_per_page)
    page_specs = list(product(row_pages, col_pages))
    pdf_path = output_dir / "figure2_ab.pdf"
    svg_paths = [output_dir / ("figure2_ab.svg" if len(page_specs) == 1 else f"figure2_ab_page{i:03d}.svg") for i in range(1, len(page_specs) + 1)]
    layout_path = output_dir / "figure2_layout.json"
    outputs = [pdf_path, *svg_paths, layout_path]
    if any(p.exists() for p in outputs):
        raise FileExistsError("Figure output already exists; use a fresh output directory")
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = result.get("sample_summary", [])
    scales = {"caller": max([v for r in summaries for v in r.get("caller_counts", {}).values() if v is not None] or [0]),
              "summary": max([r.get(k, 0) for r in summaries for k in ("union_cds", "rm_pool_cds", "retained_cds", "unassessed_cds")] or [0])}
    metadata = {"Title": title, "Author": "DefenseViz", "Subject": "Saved-call annotation evidence; no biological activity inference", "CreationDate": None, "ModDate": None}
    layouts = []
    with matplotlib.rc_context({"svg.fonttype": "none", "svg.hashsalt": "DefenseViz", "pdf.fonttype": 42, "axes.unicode_minus": False}):
        with PdfPages(pdf_path, metadata=metadata) as pdf:
            for page_index, ((row_page, col_page), svg_path) in enumerate(zip(page_specs, svg_paths), 1):
                c, layout = _make_page(result, row_page, col_page, font, page_index, len(page_specs), title, scales)
                pdf.savefig(c.fig, dpi=120)
                c.fig.savefig(svg_path, format="svg", metadata={"Title": title, "Date": None, "Creator": "DefenseViz"})
                c.fig.clear()
                layouts.append(layout)
    layout_path.write_text(json.dumps({"schema_version": 1, "page_count": len(layouts), "pages": layouts,
                                      "units": "millimeters", "font_families": sorted(font.used),
                                      "scales": scales, "bar_height_mm": BAR_HEIGHT,
                                      "matrix_scope": result.get("cell_scope", "native_final_members"),
                                      "count_semantics": "distinct observed CDSs; overlapping families are not additive",
                                      "system_counts": "Native caller-specific system counts remain in machine-readable output; no pooled system inference",
                                      "activity_inferred": False}, indent=2, ensure_ascii=False) + "\n")
    return outputs
