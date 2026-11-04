"""Live text with the approved R/grid cap-height anchor and real glyph bounds.

R/grid centers single-line text on its capital-height rectangle, not on a
font-wide ascent/descent line box. The baseline shift rotates with the label.
Only measurement uses outlines; SVG/PDF output remains ordinary editable text.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import math
from pathlib import Path as FilePath

from matplotlib import font_manager
from matplotlib.path import Path
from matplotlib.text import Text
from matplotlib.textpath import TextToPath
from matplotlib.transforms import Affine2D, Bbox

MM_PER_POINT = 25.4 / 72


@lru_cache(maxsize=8192)
def _literal_path(font_path, value):
    """Unhinted outlines at TextToPath's native 100 pt, with math disabled."""
    prop = font_manager.FontProperties(fname=font_path)
    vertices, codes = TextToPath().get_text_path(prop, value, ismath=False)
    return Path(vertices, codes) if len(vertices) else None


@lru_cache(maxsize=128)
def _cap_height(font_path):
    return _literal_path(font_path, "M").get_extents().height / TextToPath.FONT_SCALE


class ManuscriptText(Text):
    """A normal text artist whose extent is its visible, rotated glyph ink."""

    def get_window_extent(self, renderer=None, dpi=None):
        if renderer is None:
            # get_figure(root=...) was added after the supported 3.7 release.
            # A normal Figure/Canvas is sufficient here on every supported API.
            renderer = self.get_figure()._get_renderer()
        font_path = font_manager.findfont(self.get_fontproperties())
        x, y = self.get_transform().transform(self.get_position())
        if not self.get_text().strip():
            return Bbox.from_bounds(x, y, 0, 0)
        # Use the same per-line baseline and horizontal placement as Text.draw.
        # This matters for italic bearings, kerning, alignment and DPI hinting.
        _, lines, _ = self._get_layout(renderer)
        factor = renderer.points_to_pixels(1) if dpi is None else dpi / 72
        rotation = self.get_rotation()
        transform = Affine2D().scale(self.get_fontsize() / TextToPath.FONT_SCALE * factor).rotate_deg(rotation)
        bounds = []
        for line in lines:
            # Text._get_layout returns (text, metrics, x, y) through 3.10,
            # and (text, metrics, (x, y)) starting in 3.11. Both offsets are
            # already aligned/rotated; do not recompute or discard them.
            value = line[0]
            dx, dy = line[2:] if len(line) == 4 else line[2]
            path = _literal_path(font_path, value)
            if path is not None:
                bounds.append(path.get_extents(transform + Affine2D().translate(x + dx, y + dy)))
        return Bbox.union(bounds) if bounds else Bbox.from_bounds(x, y, 0, 0)


def draw_manuscript_text(c, x, y, label, pt=8.82, *, ha="left", bold=False,
                         italic=False, rotation=0, role="label", ref=None,
                         color="black", owner=None):
    """Draw one live line at an R/grid center anchor in top-left millimeters.

    ``x,y`` are the original textGrob anchor, including for rotated labels.
    Explicit lines keep the approved layout's individually measured baselines.
    Metadata records both that anchor and the derived artist baseline.
    """
    value = str(label)
    if "\n" in value or "\r" in value:
        raise ValueError("Manuscript text requires explicit single-line anchors")
    if ha not in {"left", "center", "right"}:
        raise ValueError("Text alignment must be left, center or right")
    if not all(math.isfinite(float(v)) for v in (x, y, pt, rotation)) or pt <= 0:
        raise ValueError("Text coordinates/rotation must be finite and point size positive")
    family = c.fonts.pick(value, bold)
    prop = font_manager.FontProperties(family=family, size=pt,
        weight="bold" if bold else "normal", style="italic" if italic else "normal")
    font_path = font_manager.findfont(prop, fallback_to_default=False)
    # Include the actual italic/bold face, not merely the family's regular file.
    if (family, font_path) not in c.fonts.used:
        c.fonts.used[(family, font_path)] = {"family": family, "path": font_path,
            "sha256": hashlib.sha256(FilePath(font_path).read_bytes()).hexdigest()}
    shift = _cap_height(font_path) * pt * MM_PER_POINT / 2
    angle = math.radians(rotation)
    baseline = (x + shift * math.sin(angle), y + shift * math.cos(angle))
    artist = ManuscriptText(*baseline, value, fontproperties=prop, color=color,
        ha=ha, va="baseline", rotation=rotation, rotation_mode="anchor",
        clip_on=False, parse_math=False, usetex=False)
    c.ax.add_artist(artist)
    c.text_artists.append((artist, {"text": value, "anchor_mm": [x, y],
        "baseline_mm": list(baseline), "font_pt": pt, "alignment": ha,
        "bold": bold, "italic": italic, "role": role, "owner_bounds_mm": owner,
        "source_ref": ref, "color": color, "rotation": rotation,
        "bounds_kind": "glyph_ink", "anchor_kind": "grid_cap_height_center"}))
    return artist
