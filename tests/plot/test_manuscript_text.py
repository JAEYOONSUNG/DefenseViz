import xml.etree.ElementTree as ET

import matplotlib as mpl
import numpy as np
import pytest

from defenseviz.plot.publication.drawing import Canvas, Fonts
from defenseviz.plot.publication.manuscript_text import (
    draw_manuscript_text, MM_PER_POINT,
)


def canvas(preferred="Arial"):
    return Canvas(110, 100, Fonts(preferred), "text_test", 1)


@pytest.mark.parametrize("label,pt,x,y,rotation,svg_baseline,bold", [
    ("R-M subunits and contexts", 9.24, 82.0162, 71.7, 0, (232.49,206.55),True),
    ("RMS",8.82,74.1412,84.8422,55,(212.75,242.31),False),
    ("MS",8.82,68.6287,85.7528,55,(197.12,244.89),False),
    ("Cas",9.6,75.0062,82.8594,90,(782.98-200*72/25.4,234.88),True),
])
def test_original_grid_anchors_match_live_svg_baselines(label,pt,x,y,rotation,svg_baseline,bold):
    c=canvas()
    if c.fonts.base != "Arial":pytest.skip("Exact historical Arial metric comparison")
    t=draw_manuscript_text(c,x,y,label,pt,ha="center",bold=bold,rotation=rotation)
    # Reference SVG is rounded to .01 pt; font outline and Cairo rounding differ
    # by < .004 mm. Positions are original approved anchors, not fitted offsets.
    assert t.get_position()==pytest.approx(np.array(svg_baseline)*MM_PER_POINT,abs=.004)


@pytest.mark.parametrize("rotation",[0,55,90,-90])
@pytest.mark.parametrize("ha",["left","center","right"])
def test_bounds_enclose_rotated_actual_glyphs_without_full_line_box(rotation,ha):
    c=canvas();t=draw_manuscript_text(c,55,50,"MS",8.82,ha=ha,rotation=rotation)
    bounds=c.measure()[0]["bounds_mm"]
    assert bounds[2]>0 and bounds[3]>0
    # For 90 degrees, a cap-only line has the same narrow ink dimension as the
    # unrotated line's height; unused descender padding must not be introduced.
    d=canvas();u=draw_manuscript_text(d,55,50,"MS",8.82,ha=ha)
    plain=d.measure()[0]["bounds_mm"]
    if abs(rotation)==90:
        assert bounds[2]==pytest.approx(plain[3],abs=1e-9)
        assert bounds[3]==pytest.approx(plain[2],abs=1e-9)
    assert c.ledger['texts'][0]['anchor_mm']==[55,50]
    assert t.get_parse_math() is False


def test_descenders_italic_overhang_and_real_collisions_are_measured():
    c=canvas();draw_manuscript_text(c,20,20,"M",12)
    draw_manuscript_text(c,40,20,"gy",12,italic=True)
    first,second=c.measure()
    assert second['bounds_mm'][1]+second['bounds_mm'][3]>first['bounds_mm'][1]+first['bounds_mm'][3]
    draw_manuscript_text(c,20,20,"M",12)
    with pytest.raises(ValueError,match="text_overlap"):c.validate()


def test_original_close_heading_and_rotated_groups_do_not_false_collide():
    c=canvas()
    draw_manuscript_text(c,82.0162,71.7,"R-M subunits and contexts",9.24,ha='center',bold=True)
    draw_manuscript_text(c,68.6287,75.5,"I",8.8,ha='center',bold=True)
    draw_manuscript_text(c,68.6287,85.7528,"MS",8.82,ha='center',rotation=55)
    draw_manuscript_text(c,74.1412,84.8422,"RMS",8.82,ha='center',rotation=55)
    c.validate()


def test_literal_dollars_stay_editable_text_and_italic_font_is_recorded(tmp_path):
    c=canvas();t=draw_manuscript_text(c,20,20,"$R$",12,italic=True,rotation=55)
    c.validate();out=tmp_path/'literal.svg'
    with mpl.rc_context({'svg.fonttype':'none'}):c.fig.savefig(out)
    root=ET.parse(out).getroot();texts=[e for e in root.iter() if e.tag.endswith('text')]
    assert any(''.join(e.itertext())=='$R$' for e in texts)
    assert t.get_usetex() is False
    # Windows Arial uses ariali.ttf; filename words do not identify its face.
    from matplotlib import font_manager
    resolved=font_manager.findfont(t.get_fontproperties())
    record=next(row for row in c.fonts.used.values() if row['path']==resolved)
    face=font_manager.get_font(record['path']).style_name.lower()
    assert 'italic' in face or 'oblique' in face


def test_fallback_metrics_are_used_and_invalid_inputs_rejected():
    c=canvas('FontThatDoesNotExist')
    t=draw_manuscript_text(c,20,20,'Genome A',9.24)
    assert t.get_fontfamily()[0]==c.fonts.base
    assert c.measure()[0]['font_family']==c.fonts.base
    for kw in ({'pt':0},{'rotation':float('nan')},{'ha':'bad'}):
        with pytest.raises(ValueError):draw_manuscript_text(c,20,20,'test',**kw)
    with pytest.raises(ValueError):draw_manuscript_text(c,20,20,'line\nline')


@pytest.mark.parametrize("tuple_shape",[3,4])
@pytest.mark.parametrize("rotation",[0,55,90])
def test_supported_layout_tuple_shapes_preserve_glyph_positions(monkeypatch,tuple_shape,rotation):
    c=canvas();t=draw_manuscript_text(c,55,50,"RMS",8.82,rotation=rotation,ha="right")
    c.fig.canvas.draw();renderer=c.fig.canvas.get_renderer()
    expected=t.get_window_extent(renderer).bounds
    real_layout=t._get_layout
    def layout_compat(renderer):
        box,lines,rest=real_layout(renderer)
        converted=[]
        for line in lines:
            value,metrics=line[:2]
            dx,dy=line[2:] if len(line)==4 else line[2]
            converted.append((value,metrics,dx,dy) if tuple_shape==4 else (value,metrics,(dx,dy)))
        return box,converted,rest
    monkeypatch.setattr(t,"_get_layout",layout_compat)
    assert t.get_window_extent(renderer).bounds==pytest.approx(expected,abs=1e-10)


def test_renderer_lookup_supports_legacy_get_figure_signature(monkeypatch):
    c=canvas();t=draw_manuscript_text(c,30,30,"Groups",7.3)
    c.fig.canvas.draw()
    expected=t.get_window_extent(c.fig.canvas.get_renderer()).bounds
    layout=t._get_layout(c.fig.canvas.get_renderer())
    # On a new Matplotlib, its own _get_layout calls get_figure(root=True).
    # Hold its already-computed layout to isolate this helper's renderer lookup.
    monkeypatch.setattr(t,"_get_layout",lambda renderer:layout)
    monkeypatch.setattr(t,"get_figure",lambda:c.fig)
    assert t.get_window_extent().bounds==pytest.approx(expected,abs=1e-10)
