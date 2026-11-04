"""Physical-unit, live-text drawing primitives for publication views."""
from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager, ft2font
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Rectangle
from matplotlib.textpath import TextToPath

INK = "#24393B"
TEAL = "#278481"
GRID = "#DFE6E5"
PURPLE = "#7953A2"
ASSOCIATED = "#CFE6D3"
EFFECTOR = "#206B47"
CATEGORY_COLORS = {"R-M": TEAL, "Modification": "#B09AC9", "Nuclease": "#DE8497",
                   "Cas": "#A4BC70", "CRISPR-Cas": "#A4BC70", "Abortive infection": "#80B0D3",
                   "Other": "#C2A66C", "Other mechanisms": "#C2A66C", "Unresolved": "#B9C0C2"}
SOURCE_COLORS = {"DefenseFinder": "#8C2F50", "PADLOC": "#244A70",
                 "DefensePredictor": "#585D66", "REBASE": "#66713C"}
ALIASES = {"DefenseFinder": "DF", "PADLOC": "PAD", "DefensePredictor": "DP", "REBASE": "REB"}


def stable_color(name):
    import colorsys
    hue = int(hashlib.sha256(name.encode("utf-8")).hexdigest()[:8], 16) / 0xFFFFFFFF
    return "#" + "".join(f"{round(c * 255):02X}" for c in colorsys.hls_to_rgb(hue, .43, .42))


def source_color(name, config):
    return config.get("source_colors", {}).get(name, SOURCE_COLORS.get(name, stable_color(name)))


def category_color(name, config):
    return config.get("category_colors", {}).get(name, CATEGORY_COLORS.get(name, stable_color(name)))


def alias(name, config):
    if name in config.get("aliases",{}):return config["aliases"][name]
    if name in ALIASES:return ALIASES[name]
    # Preserve the full identity in the caption/manifest; a stable key fits headers.
    return name if len(name)<=16 else "S"+hashlib.sha256(name.encode("utf-8")).hexdigest()[:4].upper()


def display_label(row):
    """Editorial qualifiers are display-only; IDs and biological metadata stay intact."""
    label=str(row.get("label",row.get("id","")))
    return re.sub(r"\s*\((?:pooled|pred\.)\)","",label,flags=re.IGNORECASE).strip()


def status_kind(state):
    """Classify explicit status, never infer absence from a missing count."""
    if isinstance(state, dict):
        state = state.get("status", state.get("state"))
    if state in ("complete", "observed", "available", "assessed", "supported", "observed_saved_assessment"):
        return "complete"
    if state in ("partial", "observed_lower_bound", "partial_coverage", "lower_bound"):
        return "partial"
    if state in ("not_applicable", "structurally_not_applicable"):
        return "not_applicable"
    return "unassessed"


def number(value, state):
    kind = status_kind(state)
    if kind == "not_applicable":
        return "-"
    if value is None or kind == "unassessed":
        return "NA"
    return f"{value:,}" + ("+" if kind == "partial" else "")


class Fonts:
    def __init__(self, preferred="Arial", fallbacks=("Liberation Sans", "DejaVu Sans")):
        self.preferred = preferred
        self.fallbacks = list(fallbacks)
        names = {f.name for f in font_manager.fontManager.ttflist}
        self.base = next((f for f in [preferred, *fallbacks] if f in names), None)
        if self.base is None:
            raise ValueError("No configured publication font is installed; install Arial or DejaVu Sans.")
        self.cache, self.used = {}, {}
        self.height_cache = {}

    def pick(self, text, bold=False):
        key = (str(text), bold)
        if key in self.cache:
            return self.cache[key]
        needed = {ord(c) for c in str(text) if not c.isspace()}
        names = [self.base, *self.fallbacks, "Noto Sans", "Noto Sans CJK KR", "Arial Unicode MS"]
        names += sorted({f.name for f in font_manager.fontManager.ttflist})
        for name in dict.fromkeys(names):
            try:
                prop = font_manager.FontProperties(family=name, weight="bold" if bold else "normal")
                path = font_manager.findfont(prop, fallback_to_default=False)
                if needed <= set(ft2font.FT2Font(path).get_charmap()):
                    self.cache[key] = name
                    self.used[(name, path)] = {"family": name, "path": path,
                        "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}
                    return name
            except (ValueError, OSError, RuntimeError):
                continue
        raise ValueError(f"No installed font covers publication label {text!r}; install a Unicode font.")

    def width(self, text, pt, bold=False):
        prop = font_manager.FontProperties(family=self.pick(text, bold), size=pt,
                                          weight="bold" if bold else "normal")
        return max((TextToPath().get_text_width_height_descent(line, prop, False)[0]
                    for line in str(text).splitlines()), default=0) * 25.4 / 72

    def height(self, text, pt, bold=False):
        family=self.pick(text,bold);key=(family,pt,bold,len(str(text).splitlines()))
        if key not in self.height_cache:
            fig=Figure(figsize=(1,1),dpi=150);canvas=FigureCanvasAgg(fig)
            artist=fig.text(.1,.5,str(text),fontsize=pt,family=family,
                            weight="bold" if bold else "normal",linespacing=1.12,parse_math=False)
            canvas.draw()
            self.height_cache[key]=artist.get_window_extent(canvas.get_renderer()).height/150*25.4
            fig.clear()
        return self.height_cache[key]

    def wrap(self, text, width, pt, bold=False):
        lines = []
        for paragraph in str(text).splitlines() or [""]:
            line = ""
            for char in paragraph:
                test = line + char
                if line and self.width(test, pt, bold) > width:
                    cut = line.rfind(" ")
                    if cut > 0:
                        lines.append(line[:cut]); line = line[cut + 1:] + char
                    else:
                        lines.append(line); line = char
                else:
                    line = test
            lines.append(line)
        return "\n".join(lines)


class Canvas:
    def __init__(self, width, height, fonts, section, page_id):
        self.width, self.height, self.fonts = width, height, fonts
        self.fig = Figure(figsize=(width / 25.4, height / 25.4), dpi=150, facecolor="white")
        FigureCanvasAgg(self.fig)
        self.ax = self.fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(height, 0))
        self.ax.axis("off")
        self.text_artists = []
        self.ledger = {"page_id": page_id, "section": section, "size_mm": [width, height],
                       "texts": [], "bars": [], "glyphs": [], "cells": [], "bands": [], "notes": []}

    def text(self, x, y, value, pt=7.8, *, ha="left", bold=False, color=INK,
             role="label", owner=None, ref=None, rotation=0, italic=False):
        value = str(value)
        t = self.ax.text(x, y, value, fontsize=pt, family=self.fonts.pick(value, bold),
                         weight="bold" if bold else "normal", color=color, ha=ha, va="center",
                         style="italic" if italic else "normal",
                         rotation=rotation, rotation_mode="anchor", linespacing=1.12,
                         clip_on=False, parse_math=False)
        self.text_artists.append((t, {"text": value, "anchor_mm": [x, y], "font_pt": pt,
            "alignment": ha, "bold": bold, "italic": italic, "role": role, "owner_bounds_mm": owner,
            "source_ref": ref, "color": color, "rotation": rotation}))
        return t

    def rect(self, x, y, w, h, fill, edge=None, lw=.3):
        self.ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge or "none", linewidth=lw))

    def line(self, x1, y1, x2, y2, color=GRID, lw=.3):
        self.ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw, solid_capstyle="butt")

    def circle(self, x, y, radius, fill, edge=None, lw=.5):
        self.ax.add_patch(Circle((x, y), radius, facecolor=fill or "none", edgecolor=edge or "none", linewidth=lw))

    def bar(self, x, y, count, denominator, state, color, *, width=12., height=1.4, ref=None, segments=None):
        kind = status_kind(state)
        self.rect(x, y - height / 2, width, height, "#EDF1F1")
        filled = 0.
        if count is not None and kind in ("complete", "partial") and denominator:
            filled = width * count / denominator
            if filled > width + 1e-8:
                raise ValueError("Bar count exceeds its declared row union")
            if segments is None:
                self.rect(x, y - height / 2, filled, height, color)
            else:
                cursor = x
                for amount, shade in segments:
                    if amount:
                        segment_width = width * amount / denominator
                        self.rect(cursor, y - height / 2, segment_width, height, shade)
                        cursor += segment_width
        self.ledger["bars"].append({"bounds_mm": [x, y-height/2, width, height], "value": count,
            "denominator": denominator, "status": kind, "filled_width_mm": filled, "color": color,
            "source_ref": ref, "segments": segments})

    def measure(self):
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        result = []
        for artist, meta in self.text_artists:
            bbox = artist.get_window_extent(renderer)
            x1, y1 = self.ax.transData.inverted().transform((bbox.x0, bbox.y1))
            x2, y2 = self.ax.transData.inverted().transform((bbox.x1, bbox.y0))
            result.append(dict(meta, font_family=artist.get_fontfamily()[0],
                               bounds_mm=[float(x1), float(y1), float(x2-x1), float(y2-y1)]))
        self.ledger["texts"] = result
        return result

    def validate(self):
        texts = self.measure()
        errors = []
        for t in texts:
            x,y,w,h = t["bounds_mm"]
            if min(x,y) < -.05 or x+w > self.width+.05 or y+h > self.height+.05:
                errors.append({"kind":"page_overflow", "text":t["text"], "bounds":t["bounds_mm"]})
            if t["owner_bounds_mm"]:
                a,b,c,d=t["owner_bounds_mm"]
                if x<a-.08 or y<b-.08 or x+w>a+c+.08 or y+h>b+d+.08:
                    errors.append({"kind":"cell_overflow", "text":t["text"], "bounds":t["bounds_mm"]})
        for i,a in enumerate(texts):
            if a["role"] == "venn_count":
                continue
            ax,ay,aw,ah=a["bounds_mm"]
            for b in texts[i+1:]:
                if b["role"] == "venn_count":
                    continue
                bx,by,bw,bh=b["bounds_mm"]
                if min(ax+aw,bx+bw)-max(ax,bx)>.12 and min(ay+ah,by+bh)-max(ay,by)>.12:
                    errors.append({"kind":"text_overlap","texts":[a["text"],b["text"]]})
        for bar in self.ledger["bars"]:
            bx,by,bw,bh=bar["bounds_mm"]
            for t in texts:
                x,y,w,h=t["bounds_mm"]
                if min(x+w,bx+bw)-max(x,bx)>.08 and min(y+h,by+bh)-max(y,by)>.08:
                    errors.append({"kind":"text_bar_overlap","text":t["text"],"bar_ref":bar["source_ref"]})
        for layer in self.ledger.get("overlaps",[]):
            if layer["mode"]!="venn4":continue
            from .overlap import _quadratic_extrema
            transform=layer["transform"];sx,sy,ox,oy=(transform[k] for k in ("sx","sy","ox","oy"))
            for t in texts:
                ref=t.get("source_ref") or {}
                if t["role"]!="venn_count" or ref.get("endpoint")!=layer["endpoint"]:continue
                x,y,w,h=t["bounds_mm"]
                box=((x-ox)/sx,(oy-y-h)/sy,(x+w-ox)/sx,(oy-y)/sy)
                for bit,ellipse in zip(ref["mask"],layer["normalized_ellipses"]):
                    lo,hi=_quadratic_extrema(box,ellipse)
                    if (bit=="1" and hi>=1) or (bit=="0" and lo<=1):
                        errors.append({"kind":"venn_region_crossing","text":t["text"],"mask":ref["mask"]})
        self.ledger["validation"]={"errors":errors}
        if errors:
            raise ValueError("Publication layout validation failed: " + repr(errors[:5]))
        return self.ledger
