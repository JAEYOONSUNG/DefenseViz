"""The approved five-card DNMB workflow, drawn as editable vector artwork.

Geometry follows the original 417.5-unit grid in its 447 × 53 mm viewport.
The Type-IV strip extends the visible artwork to 55.53 mm. Bounds uniformly
scale the complete template, including live Arial text and line widths.
This is a schematic: it neither executes nor assesses any biological workflow.
"""
from __future__ import annotations

import math
import re
from collections import defaultdict

from matplotlib.colors import to_rgba
from matplotlib.patches import PathPatch, Polygon
from matplotlib.path import Path
from matplotlib.text import Text
from matplotlib.textpath import TextPath, TextToPath
from matplotlib.transforms import Bbox

from .drawing import alias, display_label

NATIVE_WIDTH = 447.0
NATIVE_HEIGHT = 55.53
TEAL = "#287C78"
INK = "#213234"
COLORS = {"Nuclease": "#D9878C", "Modification": "#B29BCB",
          "Unknown": "#D2D4D6", "CRISPR-Cas": "#A8BF72",
          "Abortive infection": "#7FAED0"}
# Original Cairo/Arial whole-string baselines, relative to grid text anchors.
# Keeping these offsets avoids Matplotlib's different center/baseline metrics.
BASELINES = {6.7: .5520, 7.: .8245, 7.1: .8560, 7.2: .8876, 7.3: .9208,
             7.4: .9520, 7.5: .98485, 7.8: .7270, 8.: .7909, 8.1: .8224,
             8.3: .88705, 9.5: 1.09357, 11.5: 1.37983, 15.: 1.79186}
# grid::xsplinePoints(shape=1), converted back from the original physical
# viewport to the five-control-point ribbon's native offsets. No R runtime.
RIBBON_LOOP = [(-.4, 0.), (-.367084486, -.131417874), (-.218721872, -.701260126),
               (-.026304100, -1.378699708), (.141746279, -1.857720507), (.3, -2.15),
               (.480423223, -2.336722735), (.701103716, -2.447689961),
               (.955482529, -2.487695507), (1.220444927, -2.469842648),
               (1.483333333, -2.416666667), (1.763114487, -2.341378188),
               (2.074744870, -2.238571825), (2.394389906, -2.109986134),
               (2.681866522, -1.966923494), (2.916666667, -1.816666667),
               (3.112041768, -1.645323261), (3.303601345, -1.412460338),
               (3.510771077, -1.106900690), (3.666026371, -.856959140), (3.7, -.8)]


def _alpha(color, alpha):
    return to_rgba(color, alpha)


def _points(start, end, n):
    return [start + (end - start) * i / (n - 1) for i in range(n)]


class _InkText(Text):
    """Measure actual glyph ink for the original closely spaced baselines.

    Matplotlib's full font line box includes unused ascender/descender space.
    The approved source measured glyphs, so validate the same visible extent.
    The artist itself remains one ordinary editable text object in SVG/PDF.
    """
    def get_window_extent(self, renderer=None, dpi=None):
        if not self.get_text():
            return Bbox.from_bounds(0, 0, 0, 0)
        renderer = renderer or self.figure.canvas.get_renderer()
        factor = renderer.points_to_pixels(1) if dpi is None else dpi/72
        prop = self.get_fontproperties()
        bounds = TextPath((0, 0), self.get_text(), prop=prop, size=self.get_fontsize()).get_extents()
        advance = TextToPath().get_text_width_height_descent(self.get_text(), prop, False)[0]
        offset = {"left": 0., "center": advance/2, "right": advance}[self.get_horizontalalignment()]
        x, y = self.get_transform().transform(self.get_position())
        return Bbox.from_extents(x+(bounds.x0-offset)*factor, y+bounds.y0*factor,
                                 x+(bounds.x1-offset)*factor, y+bounds.y1*factor)


class _Artwork:
    def __init__(self, canvas, bounds):
        if len(bounds) != 4 or any(isinstance(v, bool) or not isinstance(v, (int, float))
                                   or not math.isfinite(v) for v in bounds):
            raise ValueError("Workflow bounds must be four finite physical coordinates")
        self.c = canvas
        self.x, self.y, width, height = bounds
        if width <= 0 or height <= 0:
            raise ValueError("Workflow bounds must have positive width and height")
        self.scale = min(width / NATIVE_WIDTH, height / NATIVE_HEIGHT)
        self.compact = False
        self.meta = {"template": "DNMB_five_card_447mm_v1", "bounds_mm": list(bounds),
                     "native_size_mm": [NATIVE_WIDTH, NATIVE_HEIGHT], "scale": self.scale,
                     "cards": [], "glyphs": [], "arrows": [], "labels": [],
                     "operation": "saved_evidence_schematic", "new_engine_runs": False,
                     "activity_assessed": False, "assay_executed": False,
                     "sample_counts_displayed": False, "raster_images": False}

    def xy(self, x, y):
        if self.compact:
            x = 8 + (x - 8) * 126 / 156
        return (self.x + x * 447 / 417.5 * self.scale,
                self.y + (y - 6) * 53 / 63 * self.scale)

    def text(self, x, y, value, size=8.2, bold=False, just="left"):
        xx, yy = self.xy(x, y)
        baseline = yy + BASELINES[size] * self.scale
        ref = {"panel": "a", "template": self.meta["template"]}
        artist = _InkText(xx, baseline, value, fontsize=size*self.scale,
                          family=self.c.fonts.pick(value, bold), weight="bold" if bold else "normal",
                          color=INK, ha=just, va="baseline", clip_on=False, parse_math=False)
        self.c.ax.add_artist(artist)
        self.c.text_artists.append((artist, {"text": value, "anchor_mm": [xx, baseline],
            "font_pt": size*self.scale, "alignment": just, "bold": bold, "role": "workflow_label",
            "owner_bounds_mm": None, "source_ref": ref, "color": INK, "rotation": 0,
            "bounds_kind": "glyph_ink"}))
        self.meta["labels"].append({"text": value, "native_anchor": [x, y],
                                    "anchor_mm": [xx, yy], "baseline_mm": baseline,
                                    "font_pt": size * self.scale, "bold": bold, "alignment": just})

    def line(self, x1, y1, x2, y2, color=TEAL, lw=.85, dashed=False):
        x1, y1 = self.xy(x1, y1); x2, y2 = self.xy(x2, y2)
        self.c.ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw * .75 * self.scale,
                       solid_capstyle="round", solid_joinstyle="round",
                       linestyle=(0, (2, 2)) if dashed else "-")

    def polyline(self, xs, ys, color=TEAL, lw=.85):
        points = [self.xy(x, y) for x, y in zip(xs, ys)]
        self.c.ax.plot(*zip(*points), color=color, linewidth=lw * .75 * self.scale,
                       solid_capstyle="round", solid_joinstyle="round")

    def polygon(self, xs, ys, fill=None, edge=TEAL, lw=.85):
        self.c.ax.add_patch(Polygon([self.xy(x, y) for x, y in zip(xs, ys)], closed=True,
            facecolor=fill if fill is not None else "none", edgecolor=edge or "none",
            linewidth=lw * .75 * self.scale, joinstyle="round"))

    def rect(self, x, y, w, h, fill=None, edge="#DEE5E6", lw=.55, dashed=False):
        xx, yy = self.xy(x, y); right, bottom = self.xy(x + w, y + h)
        self.c.rect(xx, yy, right - xx, bottom - yy, fill if fill is not None else "none",
                    edge, lw * .75 * self.scale)
        if dashed:
            self.c.ax.patches[-1].set_linestyle((0, (2, 2)))
        return [xx, yy, right - xx, bottom - yy]

    def circle(self, x, y, radius, fill=None, edge=TEAL, lw=.85):
        # R's circles use a physical mm radius, not the anisotropic native axes.
        xx, yy = self.xy(x, y)
        self.c.circle(xx, yy, radius * self.scale, fill, edge, lw * .75 * self.scale)

    def arc(self, x, y, rx, ry, start=0, end=math.pi, color=TEAL, lw=.85):
        angles = _points(start, end, 37)
        self.polyline([x + rx * math.cos(a) for a in angles],
                      [y + ry * math.sin(a) for a in angles], color, lw)

    def ellipse(self, x, y, rx, ry, fill, color=TEAL, lw=.85):
        angles = _points(0, 2 * math.pi, 73)
        self.polygon([x + rx * math.cos(a) for a in angles],
                     [y + ry * math.sin(a) for a in angles], fill, color, lw)

    def rounded(self, x, y, w, h, r, fill, color=TEAL, lw=.85):
        points = []
        for i, (cx, cy) in enumerate(((x+w-r, y+r), (x+w-r, y+h-r), (x+r, y+h-r), (x+r, y+r))):
            points += [(cx+r*math.cos(a), cy+r*math.sin(a))
                       for a in _points(-math.pi/2+i*math.pi/2, -math.pi/2+(i+1)*math.pi/2, 17)]
        self.polygon(*zip(*points), fill, color, lw)

    def arrow(self, xs, ys, color=TEAL, lw=.9, head=1.1, dashed=False, name=None):
        route = [self.xy(x, y) for x, y in zip(xs, ys)]
        lengths = [math.dist(a, b) for a, b in zip(route, route[1:])]
        if not lengths or min(lengths) <= 0:
            raise ValueError("Workflow arrow needs distinct consecutive points")
        segments = [route]
        if dashed:
            dash = 2 * lw * .75 * 25.4 / 72 * self.scale
            cumulative = [0.]
            for length in lengths:
                cumulative.append(cumulative[-1] + length)
            def at(distance):
                i = next((i for i, end in enumerate(cumulative[1:]) if distance <= end), len(lengths)-1)
                t = (distance-cumulative[i])/lengths[i]
                return tuple(a+(b-a)*t for a, b in zip(route[i], route[i+1]))
            segments = []
            start = 0.
            while start < cumulative[-1]:
                end = min(start+dash, cumulative[-1])
                segments.append([at(start), *[route[i] for i in range(1, len(route)-1)
                                             if start < cumulative[i] < end], at(end)])
                start += 2*dash
        dx, dy = (route[-1][i]-route[-2][i] for i in (0, 1))
        length = math.hypot(dx, dy); ux, uy = dx/length, dy/length
        size = head*self.scale; tip = route[-1]
        base = (tip[0]-size*math.cos(math.pi/6)*ux, tip[1]-size*math.cos(math.pi/6)*uy)
        triangle = [(base[0]-size*.5*uy, base[1]+size*.5*ux), tip,
                    (base[0]+size*.5*uy, base[1]-size*.5*ux)]
        vertices, codes = [], []
        for segment in [*[r + r[-2::-1] for r in segments], triangle+[triangle[0]]]:
            vertices += segment
            codes += [Path.MOVETO] + [Path.LINETO]*(len(segment)-2) + [Path.CLOSEPOLY]
        self.c.ax.add_patch(PathPatch(Path(vertices, codes), facecolor="none", edgecolor=color,
                                     lw=lw*.75*self.scale, capstyle="round", joinstyle="round"))
        self.meta["arrows"].append({"id": name, "route_mm": route, "head_mm": triangle,
                                     "dashed": dashed, "compound_path": True})

    def glyph(self, kind, x, y, w, h):
        left, top = self.xy(x, y); right, bottom = self.xy(x+w, y+h)
        self.meta["glyphs"].append({"kind": kind, "bounds_mm": [left, top, right-left, bottom-top],
                                    "illustrative": True})

    def database(self, x, y, color=TEAL):
        self.ellipse(x, y+6.6, 5.1, 1.45, "#E8F1F0", color)
        self.rect(x-5.1, y+1.4, 10.2, 5.2, "#E8F1F0", None)
        for xx in (x-5.1, x+5.1): self.line(xx, y+1.4, xx, y+6.6, color)
        for dy in (3.1, 4.9, 6.6): self.arc(x, y+dy, 5.1, 1.45, color=color, lw=.7)
        self.ellipse(x, y+1.4, 5.1, 1.45, "#F9FCFB", color)
        self.glyph("database_layered_cylinder", x-5.2, y-.15, 10.4, 8.3)

    def gene(self, x, y, label, color, w=8):
        self.polygon([x+a for a in (0, w-1.8, w, w-1.8, 0)],
                     [y+a for a in (-1.45, -1.45, 0, 1.45, 1.45)], _alpha(color, .18), color, .8)
        self.text(x+w/2-.25, y, label, 7.5, True, "center")
        self.glyph("CDS_directional_arrow", x, y-1.45, w, 2.9)

    def dna(self, x, y, w, color):
        u = _points(0, w, 121); v = [math.sin(a*math.pi/5)*1.3 for a in u]
        for i in range(5, 116, 10): self.line(x+u[i], y+v[i], x+u[i], y-v[i], _alpha(color, .7), .65)
        for sign in (-1, 1): self.polyline([x+a for a in u], [y+sign*a for a in v], color, 1.05)
        self.glyph("DNA_duplex", x, y-1.4, w, 2.8)

    def domain(self, x, y, w=12, color=TEAL):
        self.line(x-w/2-.3, y-1, x+w/2+.3, y-1, _alpha(color, .6), .7)
        for start, width, opacity in zip((-.5, -.13, .28), (.25, .29, .18), (.14, .32, .55)):
            self.rounded(x+start*w, y-2.2, width*w, 2.4, .5, _alpha(color, opacity), color, .75)
        for i in range(7):
            self.rect(x-w*.34+i*w*.105, y+1.4, w*.075, 1.3,
                      _alpha(color, .72) if i in (1, 3, 5) else "white", color, .45)
        self.glyph("domain_profile_and_motif_slots", x-w/2-.4, y-2.3, w+.8, 5.2)

    def components(self, x, y, color):
        self.line(x-4, y, x+4, y, _alpha(color, .7), .7)
        for i in (-1, 0, 1):
            self.circle(x+i*3.3, y, 1.1, "white", color)
            self.circle(x+i*3.3, y, .45, _alpha(color, .6), None)
        self.glyph("component_context", x-4.5, y-1.2, 9, 2.4)

    def funnel(self, x, y):
        self.line(x-5, y, x+5, y, "#C6D1D0", 2.2)
        self.line(x, y, x+5, y, _alpha(TEAL, .7), 2.2)
        for xx in (x-5, x+5): self.line(xx, y-.8, xx, y+.8, "#9AAEAB", .65)
        self.line(x, y-3, x, y+3)
        self.polygon([x-.65, x, x+.65], [y-3.1, y-2.2, y-3.1], TEAL, None)
        self.glyph("final_score_threshold_tick", x-5.2, y-3.3, 10.4, 6.5)

    def model(self, x, y):
        u = _points(0, 1, 121); xx = [x-4.9+4.5*a for a in u]; yy = [y+math.sin(6*math.pi*a)*1.45 for a in u]
        self.polyline(xx, yy, "#D5E6E2", 2.7); self.polyline(xx, yy, TEAL, .95)
        self.polyline([x+dx for dx, _ in RIBBON_LOOP], [y+dy for _, dy in RIBBON_LOOP], TEAL, .85)
        self.polygon([x+a for a in (2.5, 3.7, 3.7, 4.6, 3.1, 1.6, 2.5)],
                     [y+a for a in (-1, -1, 1, 1, 2.7, 1, 1)], "#93B8AC", TEAL, .8)
        self.glyph("protein_ribbon_schematic", x-5.1, y-2.9, 10, 5.8)

    def pair(self, x, y):
        left = [(x+a, y+b) for a, b in zip((-4.7, -3.5, -2.2, -3.7), (-1.7, -.8, 1.1, 2.1))]
        right = [(x+a, y+b) for a, b in zip((2.4, 4, 4.9, 3.5), (-1.7, -.7, 1.1, 2.3))]
        for points in (left, right): self.polyline(*zip(*points), "#9CB7B1", .75)
        a, b = self.xy(*left[2]), self.xy(*right[0])
        self.c.ax.plot([a[0], b[0]], [a[1], b[1]], color=TEAL,
                       lw=.7*.75*self.scale, linestyle=(0, (4, 4)))
        for tip, other in ((a, b), (b, a)):
            dx, dy = (tip[i]-other[i] for i in (0, 1)); length = math.hypot(dx, dy)
            ux, uy = dx/length, dy/length; size = .6*self.scale
            base = (tip[0]-size*math.cos(math.pi/6)*ux, tip[1]-size*math.cos(math.pi/6)*uy)
            points = [(base[0]-size*.5*uy, base[1]+size*.5*ux), tip,
                      (base[0]+size*.5*uy, base[1]-size*.5*ux)]
            self.c.ax.add_patch(PathPatch(Path(points, [Path.MOVETO, Path.LINETO, Path.LINETO]),
                                         facecolor="none", edgecolor=TEAL, lw=.7*.75*self.scale))
        for points, selected in ((left, 2), (right, 0)):
            for i, point in enumerate(points):
                self.circle(*point, .82 if i == selected else .5,
                            TEAL if i == selected else "#D9E8E3", "white" if i == selected else TEAL, .55)
        self.glyph("native_residue_pair_distance", x-5.3, y-2.4, 10.7, 6.7)

    def cas(self, x, y, color):
        self.line(x-5, y, x+5, y, color, .65)
        for i in (-1, 0, 1): self.rect(x+i*3.4-1.05, y-1.5, 2.1, 3, _alpha(color, .4), color, .75)
        self.glyph("CRISPR_array_context", x-5.1, y-1.6, 10.2, 3.2)

    def cell(self, x, y, color):
        self.rounded(x-4.5, y-2.75, 9, 5.5, 2.75, _alpha(color, .16), color)
        self.arc(x-1.2, y, 1.25, .95, 0, 1.7*math.pi, color, .7)
        self.circle(x+2, y, .6, color, None)
        self.glyph("cell_context", x-4.6, y-2.9, 9.2, 5.8)


def _workflow_parameters(result, options):
    """Resolve displayed criteria without allowing presentation to revise them."""
    def mapping(value, name):
        if value is None:
            return {}
        if not isinstance(value, dict):
            raise ValueError(f"Workflow {name} must be an object")
        return value

    selected = mapping(result.get("selected_policy"), "selected_policy")
    policy = mapping(result.get("policy"), "policy")
    layers = [("selected_policy.parameters", mapping(selected.get("parameters"), "selected_policy.parameters")),
              ("policy", policy), ("policy.parameters", mapping(policy.get("parameters"), "policy.parameters"))]
    if result.get("retain_cutoff") is not None:
        cutoff = result["retain_cutoff"]
        if isinstance(cutoff, list) and len(cutoff) == 1:
            cutoff = cutoff[0]
        layers.append(("retain_cutoff", {"retain_cutoff": cutoff}))
    layers.append(("configuration", mapping(options.get("parameters"), "parameter configuration")))
    values, origins = {}, {}
    for origin, layer in layers:
        for key in ("structure", "max_divergent_gap_bp", "retain_cutoff"):
            value = layer.get(key)
            if value is None:
                continue
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f"Workflow {key} must be an explicit nonnegative finite value")
            if key in values and value != values[key]:
                raise ValueError(f"Conflicting workflow {key} in {origins[key]} and {origin}; saved criteria cannot be overridden")
            if key not in values:
                values[key], origins[key] = value, origin
    return values, origins


def _workflow_providers(result, config):
    """Use declared identities; an empty standalone result is a template preview."""
    original = ["PADLOC", "DefenseFinder", "DefensePredictor", "REBASE"]
    if not result:
        providers, origin = original, "schematic_reference"
    else:
        declared = result.get("sources", result.get("callers", []))
        if not isinstance(declared, list):
            raise ValueError("Workflow declared sources must be a list")
        ids, excluded = [], set()
        for source in declared:
            if isinstance(source, str):
                sid, kind = source, None
            elif isinstance(source, dict):
                sid, kind = source.get("id"), source.get("kind")
            else:
                raise ValueError("Workflow sources require explicit source identifiers")
            if not isinstance(sid, str) or not sid.strip():
                raise ValueError("Workflow sources require explicit source identifiers")
            if kind in {"hmm_evidence", "hmm_only", "hmm_hit"}:
                excluded.add(sid)
            else:
                ids.append(sid)
        model = result.get("manuscript_atlas")
        if isinstance(model, dict) and "source_order" in model and "reference_sources" in model:
            native, references = model["source_order"], model["reference_sources"]
            if any(not isinstance(x, list) or any(not isinstance(s, str) or not s.strip() for s in x)
                   for x in (native, references)):
                raise ValueError("Workflow atlas providers require explicit source lists")
            providers, origin = native + references, "manuscript_atlas"
        else:
            providers, origin = ids, "declared_sources"
        providers = list(dict.fromkeys(s for s in providers if s not in excluded))
        if set(providers) == set(original):
            providers = original
    labels = [{"id": sid, "label": ("REBASE" if sid == "REBASE" and sid not in config.get("aliases", {})
                                      else alias(sid, config))} for sid in providers]
    if any(not isinstance(row["label"], str) or not row["label"].strip() for row in labels):
        raise ValueError("Workflow source aliases must be nonempty text")
    return labels, origin


def _draw_provider_label(a, providers):
    label = " | ".join(row["label"] for row in providers) or "Sources not declared"
    left, top = a.xy(9, 41)
    right, bottom = a.xy(43, 50)
    width = (right-left)/a.scale
    lines = a.c.fonts.wrap(label, width, 6.7).splitlines()
    if len(lines) > 2:
        raise ValueError("Workflow source labels exceed the template; provide shorter source aliases")
    positions = [44+3*i for i in range(len(lines))]
    for line, y in zip(lines, positions):
        a.text(26, y, line, 6.7, just="center")
        a.c.text_artists[-1][1]["owner_bounds_mm"] = [left, top, right-left, bottom-top]
    a.meta["provider_label_lines"] = lines
    a.meta["provider_label_bounds_mm"] = [left, top, right-left, bottom-top]


def _pt_display_label(row):
    """Compact only editorial display text; retain original family identities."""
    return re.sub(r"\s*\(PT\)\s*$", "", display_label(row), flags=re.IGNORECASE).strip()


def _workflow_pt_pairs(result, options, fonts):
    """Select authored display pairs, without claiming biological co-membership."""
    meta = {"origin": "observed_family_display_metadata", "selected": [], "omitted": [],
            "scope": "Illustrative labels from explicit display metadata; not system membership or activity."}
    if "pt_pairs" in options or "workflow_pt_pairs" in result:
        configured = "pt_pairs" in options
        pairs = options["pt_pairs"] if configured else result["workflow_pt_pairs"]
        meta["origin"] = "configuration" if configured else "result.workflow_pt_pairs"
        if not isinstance(pairs, list) or len(pairs) > 2 or any(not isinstance(p, (list, tuple)) or len(p) != 2
                or any(not isinstance(v, str) or not v.strip() for v in p) for p in pairs):
            raise ValueError("Workflow pt_pairs must contain at most two explicit label pairs")
        candidates = [{"labels": list(pair), "display_labels": [_pt_display_label({"label": value}) for value in pair],
                       "source_ref": {"origin": meta["origin"], "pair_index": i}} for i, pair in enumerate(pairs)]
    else:
        inventory = result.get("inventory_rows") or {}
        families = inventory.get("families", []) if isinstance(inventory, dict) else []
        if not isinstance(families, list):
            raise ValueError("Workflow inventory families must be a list")
        groups = defaultdict(list)
        for row in families:
            if not isinstance(row, dict):
                raise ValueError("Workflow inventory families must be objects")
            group = row.get("partner_group")
            if isinstance(group, str) and group.strip():
                groups[group].append(row)
        candidates = []
        def group_order(item):
            orders = [r.get("display_order") for r in item[1]]
            known = [v for v in orders if not isinstance(v, bool) and isinstance(v, (int, float)) and math.isfinite(v)]
            return min(known, default=math.inf), item[0]
        for group, rows in sorted(groups.items(), key=group_order):
            record = {"partner_group": group, "family_ids": [r.get("id") for r in rows],
                      "labels": [display_label(r) for r in rows]}
            observed = [r for r in rows if r.get("visible") is not False and
                        isinstance(r.get("N"), (int, float)) and not isinstance(r.get("N"), bool) and
                        math.isfinite(r["N"]) and r["N"] > 0 and float(r["N"]).is_integer()]
            if len(observed) != 2:
                meta["omitted"].append({**record, "reason": "requires_two_observed_visible_families",
                                        "observed_family_ids": [r.get("id") for r in observed]})
                continue
            orders = [r.get("partner_order") for r in observed]
            if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in orders) or orders[0] == orders[1]:
                meta["omitted"].append({**record, "reason": "missing_or_ambiguous_partner_order"})
                continue
            if not all(r.get("display_group") == "Modification" or r.get("category") == "Modification" or
                       isinstance(r.get("category_tags"), list) and "Modification" in r["category_tags"] for r in observed):
                meta["omitted"].append({**record, "reason": "not_declared_modification_display_pair"})
                continue
            ordered = sorted(observed, key=lambda r: r["partner_order"])
            if any(not isinstance(r.get("id"), str) or not r["id"].strip() for r in ordered):
                meta["omitted"].append({**record, "reason": "missing_family_identity"})
                continue
            if ordered[0]["id"] == ordered[1]["id"]:
                meta["omitted"].append({**record, "reason": "duplicate_family_identity"})
                continue
            if any(not isinstance(r.get("_display_label", r.get("label", r.get("id"))), str) for r in ordered):
                meta["omitted"].append({**record, "reason": "missing_display_label"})
                continue
            labels = [_pt_display_label(r) for r in ordered]
            if any(not label for label in labels):
                meta["omitted"].append({**record, "reason": "missing_display_label"})
                continue
            candidates.append({"partner_group": group, "family_ids": [r["id"] for r in ordered],
                "partner_orders": [r["partner_order"] for r in ordered], "observed_CDS": [r["N"] for r in ordered],
                "labels": [display_label(r) for r in ordered], "display_labels": labels,
                "source_ref": {"table": "inventory_rows.families", "family_ids": [r["id"] for r in ordered]}})
        if len(candidates) > 2:
            meta["omitted"].extend({**record, "reason": "more_than_two_display_pairs_no_subset_selected"} for record in candidates)
            candidates = []
    selected = []
    # Physical widths at the native scale; all geometry and typography scale
    # together. Preserve the 21/25-unit boxes and their original 7.3pt type.
    available = [(width-1.6)*447/417.5 for width in (21, 25)]
    for record in candidates:
        widths = [fonts.width(label, 7.3, True) for label in record["display_labels"]]
        if any(width > limit for width, limit in zip(widths, available)):
            if meta["origin"] != "observed_family_display_metadata":
                raise ValueError("Workflow pt_pairs labels exceed the template; provide shorter explicit display labels")
            meta["omitted"].append({**record, "reason": "labels_exceed_mini_boxes", "widths_mm": widths,
                                    "available_widths_mm": available})
            continue
        selected.append(record["display_labels"])
        meta["selected"].append({**record, "widths_mm": widths, "available_widths_mm": available})
    return selected, meta


def draw_manuscript_workflow(c, result, config, bounds):
    """Draw the five-card original in ``[left, top, width, height]`` mm.

    ``config['manuscript_workflow']`` may provide ``parameters`` (structure,
    max_divergent_gap_bp, retain_cutoff) and up to two explicit ``pt_pairs``.
    Saved ``result.selected_policy.parameters`` and ``result.policy`` are
    authoritative; config may fill missing criteria or repeat matching values.
    Conflicting values are rejected. Missing criteria retain saved-review labels.
    Providers come from the enriched model or declared source identities. An
    empty result renders the reference schematic. Without explicit pt_pairs,
    observed family partner display metadata may supply two mini-box pairs.
    No biological pairing or assay outcome is inferred.
    """
    options = config.get("manuscript_workflow", {})
    if not isinstance(options, dict):
        raise ValueError("manuscript_workflow must be an object")
    parameters, parameter_origins = _workflow_parameters(result, options)
    providers, provider_origin = _workflow_providers(result, config)
    pairs, pair_meta = _workflow_pt_pairs(result, options, c.fonts)
    a = _Artwork(c, bounds)
    a.meta.update(parameters={k: parameters.get(k) for k in ("structure", "max_divergent_gap_bp", "retain_cutoff")},
                  parameter_origins=parameter_origins, providers=providers, provider_origin=provider_origin,
                  pt_pairs=pairs, pt_pair_selection=pair_meta, illustrative_validation_route=True)
    text, line, arrow = a.text, a.line, a.arrow
    def criterion(key, prefix, suffix, absent):
        return f"{prefix}{parameters[key]:g}{suffix}" if parameters.get(key) is not None else absent
    def card(number, x, w, fill, color, title, badge_x=None, dashed=False):
        box = a.rect(x, 22, w, 45, fill, color, .7 if dashed else .65, dashed)
        bx = x+3 if badge_x is None else badge_x
        a.rect(bx-.65, 25.1, 3, 3.8, _alpha("#C7CCD0" if number == 2 else color, .25),
               "#C7CCD0" if number == 2 else color, .6)
        text(bx, 27, f"{number}  {title}", 9.5, True)
        a.meta["cards"].append({"number": number, "title": title, "bounds_mm": box})

    text(3/(447/417.5), 10, "a", 15, True)
    text(10/(447/417.5), 10, "Evidence criteria for DNMB-guided barrier prioritization", 11.5, True)
    a.compact = True
    card(1, 8, 156, _alpha(TEAL, .035), TEAL, "R-M evidence review")
    a.database(26, 30)
    text(26, 40, "Integrated R-M evidence", 7.2, True, "center")
    _draw_provider_label(a, providers)
    line(32, 34, 44, 34, lw=.8)
    for yy in (31.5, 40.5):
        line(44, 34, 44, yy, lw=.8); arrow([44, 47], [yy, yy], lw=.8)
        line(71, yy, 74, yy, lw=.8); line(74, yy, 74, 34, lw=.8)
    a.rect(47, 29, 24, 5, "#E3EDEB", "#68A19C", .6)
    text(59, 31.5, "Intersection", 6.7, True, "center")
    a.rect(47, 38, 24, 5, "#F5F7F7", "#9EAEB0", .6)
    text(59, 40.5, "Remainder", 6.7, just="center")
    arrow([74, 83], [34, 34], name="rm_pool_to_evidence")
    a.domain(90, 34)
    text(90, 40, "Reference / HMM / motif", 7.2, True, "center")
    text(90, 44, "InterPro / eggNOG review", 6.7, just="center")
    arrow([103, 117], [34, 34], name="rm_evidence_to_context")
    for x, label, color in ((122, "R", COLORS["Nuclease"]), (132, "M", COLORS["Modification"]), (142, "S", COLORS["Unknown"])):
        a.gene(x, 34, label, color)
    line(122, 37, 150, 37, lw=.5)
    text(136, 40, "Cognate context", 7.2, True, "center")
    text(136, 44, criterion("max_divergent_gap_bp", "Neighborhood <=", " bp", "Neighborhood review"), 6.7, just="center")
    arrow([151, 158, 158, 148], [34, 34, 55, 55], name="rm_context_to_structure")
    a.model(136, 55)
    text(136, 61, "Structures", 7.5, True, "center")
    text(136, 65, "AFDB reuse / ESMFold2-Fast", 6.7, just="center")
    arrow([124, 98], [55, 55], name="rm_structure_to_geometry")
    a.pair(84, 55)
    text(84, 59, "Type-specific sites", 7.1, True, "center")
    text(84, 62, "Geometry + confidence", 7.1, True, "center")
    text(84, 65, criterion("structure", "Structure evidence +", "/0", "Structure evidence (saved)"), 6.7, just="center")
    arrow([72, 44], [55, 55], name="rm_geometry_to_retention")
    a.funnel(32, 55)
    text(32, 61, "Assigned CDS", 7.5, True, "center")
    text(32, 51.3, criterion("retain_cutoff", "Score >=", "", "Score (saved)"), 6.7, just="center")
    a.rect(8, 68, 156, 4, "#F1F3F4", "#A9B2B7", .5)
    text(10, 70, "Type IV candidates", 6.7, True)
    arrow([49, 55], [70, 70], "#65777B", .7, .7, name="typeiv_to_review")
    text(58, 70, "Exact-domain component review", 6.7)
    arrow([107, 113], [70, 70], "#65777B", .7, .7, name="typeiv_to_annotation")
    text(116, 70, "R / A annotation", 6.7)

    a.compact = False
    card(2, 137, 56, "#F8FAF4", "#829764", "Proposed validation", dashed=True)
    a.database(145, 30, "#829764")
    text(155, 34, "MethREfinder", 8.1, True)
    text(165, 42, "Heterologously express the", 8.1, just="center")
    text(165, 45.5, "context-linked M(S) subset", 8.1, just="center")
    arrow([165, 165], [47.4, 49], "#829764", .7, .5, name="proposed_expression_to_assay")
    text(165, 51, "Methylation-sensitive restriction assay", 8.1, just="center")
    text(165, 54.6, "Available commercial enzymes +", 8.1, just="center")
    text(165, 58.2, "matching recognition sequences", 8.1, just="center")
    arrow([165, 165], [60, 61.5], "#829764", .7, .5, name="proposed_assay_to_assessment")
    text(165, 64, "R-M evasion assessment", 8.1, True, "center")

    modification, nuclease = COLORS["Modification"], COLORS["Nuclease"]
    card(3, 196, 58, _alpha(modification, .07), modification, "Modification-associated review")
    a.dna(199, 34, 12, modification); a.circle(203, 31, .6, modification, None)
    text(215, 34, "Modification candidates", 7.3, True)
    arrow([205, 205], [37, 39], modification, .85, .8, name="modification_to_domains")
    a.domain(205, 42.2, 8, modification)
    text(215, 41, "Nuclease-domain / HMM", 7.4, True)
    text(215, 45, "+ available motif evidence", 7.1)
    arrow([205, 205], [47, 50], modification, .85, .8, name="domains_to_context")
    a.components(205, 54, modification)
    text(215, 54, "Component / context review", 7.2, True)
    arrow([254.3, 256.5], [54, 54], modification, .85, .8, name="context_to_nuclease")
    if pairs:
        left_box = a.rect(199, 58.3, 21, 7.6, _alpha(modification, .13), modification, .65)
        right_box = a.rect(226, 58.3, 25, 7.6, _alpha(nuclease, .13), nuclease, .65)
        a.meta["pt_pair_label_boxes_mm"] = [left_box, right_box]
        for i, (yy, (left, right)) in enumerate(zip([62.05] if len(pairs) == 1 else [60.2, 63.9], pairs)):
            text(209.5, yy, left, 7.3, True, "center")
            c.text_artists[-1][1].update(owner_bounds_mm=left_box,
                source_ref={"panel": "a", **pair_meta["selected"][i]["source_ref"], "pair_side": "left"})
            line(220.6, yy, 225.4, yy, "#65777B", .7, True)
            text(238.5, yy, right, 7.3, True, "center")
            c.text_artists[-1][1].update(owner_bounds_mm=right_box,
                source_ref={"panel": "a", **pair_meta["selected"][i]["source_ref"], "pair_side": "right"})

    card(4, 257, 50, _alpha(nuclease, .06), nuclease, "Nuclease evidence", badge_x=259)
    a.dna(261, 38, 17, nuclease); a.domain(295, 37.5, 11, nuclease)
    text(259, 49, "Effector / domain / motif", 8.3, True)
    text(259, 56, "Restriction function reviewed", 7.8)
    card(5, 310, 73, None, "#AAB7BC", "Complementary defense priorities")
    a.cas(318, 35, COLORS["CRISPR-Cas"]); text(327, 33, "Cas", 8., True)
    text(327, 37, "Programmable; CRISPR-array-dependent", 7.)
    text(313, 43, "Vector-targeting compatibility: review required", 7., True)
    a.cell(318, 51, COLORS["Abortive infection"]); text(326, 51, "Abortive infection", 7.8)
    text(313, 59, "Lower priority for domestication", 7.5, True)
    for i, (left, right) in enumerate(((134, 137), (193, 196), (254, 257), (307, 310))):
        arrow([left+.25, right-.25], [30, 30], "#526C72", 1.2, dashed=i < 2, name=f"stage_{i+1}_to_{i+2}")
    c.ledger["manuscript_workflow"] = a.meta
    return a.meta
