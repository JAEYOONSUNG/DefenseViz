"""Optional Letter composition with measured continuation pages.

The atlas is drawn at its established native size and uniformly exported to
Letter landscape. Export never alters measurements, source memberships or the
fixed sizes of the separate-page layout. Exported strokes are kept at or above
a physical floor; each page records its export scale and measured smallest
text and stroke. Continuation type tables are arranged to use their own Letter
page rather than inheriting the atlas width.
"""
from copy import deepcopy

from matplotlib.colors import to_rgba
from matplotlib.patches import Patch

from .atlas import (_chunks, _type_cell, _type_case_label_layout, _type_case_label, _tool_header,
                    _pack_type_legend, _draw_type_legend, _type_legend_items,
                    EXTERNAL_CASE_FILL, PRIMARY_CASE_FILL)
from .comparison import comparison_pages, source_order
from .drawing import STROKE_FLOOR_PT, alias, status_kind, stroke_artists
from .manuscript_atlas import manuscript_atlas_pages
from .manuscript_text import draw_manuscript_text

LETTER = (279.4, 215.9)
TYPE_TITLE = 'CDS inclusion by assigned R-M type'
TITLE_PT, LEGEND_PT, LEGEND_STEP = 11.5, 9.4, 3.6
LEGEND_OFFSET, TABLE_OFFSET = 3.8, 3.4      # title anchor -> first legend line -> tables
HEADING_CLEARANCE, TABLE_CLEARANCE = 4., 2.5
TABLE_GAP, MARGIN_X, MARGIN_TOP, MARGIN_BOTTOM = 6., 10., 5., 1.6
# Presentation preference for continuation pages (numbers ~9 pt at Letter),
# not a data threshold. Pages fall back to the best scale >= composite_min_scale.
CONTINUATION_TARGET_SCALE = .85


def _type_plans(result, fonts, config):
    rows = result.get('type_tables', [])
    samples = {s['sample_id']: s for s in result.get('samples', [])}
    selected = config.get('type_samples')
    ids = [sid for sid in (samples if selected is None else selected)
           if any(r['sample_id'] == sid for r in rows)]
    selection = 'explicit' if selected is not None else 'all_typed_default'
    if selected is None:
        highlighted = [sid for sid in ids if samples[sid].get('highlight')]
        if highlighted:
            ids = highlighted
            selection = 'highlighted_default'
    plans = []
    if not ids:
        return plans
    shown = [r for r in rows if r['sample_id'] in set(ids)]
    for sources in _chunks(source_order(result, config), 3):
        fields = ['Type'] + (['Pool'] if any('pool_CDS' in r for r in rows) else []) + sources + ['DNMB', 'R', 'M', 'S', 'A']
        base = {'Type': 10., 'Pool': 12., 'DNMB': 15., 'R': 13., 'M': 16., 'S': 14., 'A': 13.3333333333}
        def measured_value(row, field):
            value = _type_cell(row, field)['display']
            if field == 'Type' and value.startswith('Type '): value = value[5:]
            return fonts.width(value, 10.8, True) + (1. if field == 'Type' else 3.)
        # Columns are sized to the tables that are actually displayed.
        widths = {field: max(base.get(field, 15.), fonts.width(alias(field, config), 8.6, True)+3,
                     max((measured_value(r, field) for r in shown), default=0)) for field in fields}
        casew = sum(widths.values())
        width_budget = LETTER[0]/config['composite_min_scale']
        perpage = min(3, max(1, int((width_budget-20.+6.)/(casew+6.))))
        header = _tool_header(fonts, fields, widths, 8.5)
        for cases in _chunks(ids, perpage):
            labels = {sid: fonts.wrap(samples[sid].get('label', sid), casew-4, 11.4, True) for sid in cases}
            # The label layout reserves 2 mm vertical padding; 1 mm keeps the
            # measured label clear of the band. The reduced reservation stays
            # in the styled layout so its own fit check keeps explicit italics.
            styled = {}
            for sid in cases:
                layout = _type_case_label_layout(fonts, casew, samples[sid], 11.4)
                styled[sid] = None if layout is None else dict(layout, height=layout['height']-1.)
            strain = {sid: max(5., styled[sid]['height'] if styled[sid] is not None
                               else fonts.height(labels[sid], 11.4, True)+1.) for sid in cases}
            plans.append(dict(cases=cases, fields=fields, widths=widths, case_width=casew, strain=strain,
                              labels=labels, styled=styled, samples=samples, rows=rows, selection=selection,
                              tool_header=header))
    return plans


def _row_strain(plan, sids):
    return max(plan['strain'][sid] for sid in sids)


def _row_height(plan, sids):
    """Tables in one row share the tallest label band and row count."""
    rowcount = max(sum(r['sample_id'] == sid for r in plan['rows']) for sid in sids)
    return _row_strain(plan, sids) + 3. + 4.8 + rowcount*5.6


def _legend(fonts, plan, cases, width):
    return _pack_type_legend(fonts, _type_legend_items(cases, plan['samples']), width, LEGEND_PT)


def _heading_extent(fonts, legend, title):
    widths = [fonts.width(title, TITLE_PT, True)]
    widths += [sum(fonts.width(text, LEGEND_PT) + (3.2 if fill else 0.) for text, fill in line)
               + fonts.width(' · ', LEGEND_PT)*(len(line)-1) for line in legend]
    return max(widths)


def _table_offset(legend):
    return LEGEND_OFFSET + (len(legend)-1)*LEGEND_STEP + TABLE_OFFSET


def _draw_heading(c, x, y, legend, plan, cases, continued):
    """Panel letter, title and measured legend lines at a cap-height anchor."""
    draw_manuscript_text(c, x-7, y, 'c', 15, bold=True)
    draw_manuscript_text(c, x, y, TYPE_TITLE + (' (continued)' if continued else ''), TITLE_PT, bold=True)
    _draw_type_legend(c, x, y+LEGEND_OFFSET, legend, LEGEND_PT, LEGEND_STEP,
                      lambda xx, yy, text, pt: draw_manuscript_text(c, xx, yy, text, pt))
    c.ledger['type_sample_ids'] = c.ledger.get('type_sample_ids', []) + list(cases)
    c.ledger['type_sample_selection'] = plan['selection']
    c.ledger['type_source_order'] = [f for f in plan['fields'] if f not in ('Type', 'Pool', 'DNMB', 'R', 'M', 'S', 'A')]
    c.ledger['type_continued'] = continued
    c.ledger['type_tool_header'] = plan['tool_header'][0]
    c.ledger['type_context_note'] = 'Counts are annotation assignments. Context does not establish biological activity.'


def _draw_case(c, plan, sid, x, top, config, strainh):
    """One type table; values come only from the validated scoped results."""
    t = lambda xx, yy, value, pt=10.8, **kw: draw_manuscript_text(c, xx, yy, value, pt, **kw)
    casew = plan['case_width']
    sample = plan['samples'][sid]
    external = sample.get('cohort') in ('external', 'external_validation')
    c.rect(x, top, casew, strainh, EXTERNAL_CASE_FILL if external else PRIMARY_CASE_FILL)
    _type_case_label(c, x, top, casew, strainh, sample, plan['labels'][sid], plan['styled'][sid], 11.4)
    groupy = top+strainh
    rolesw = sum(plan['widths'][f] for f in ('R', 'M', 'S', 'A'))
    c.rect(x, groupy, casew-rolesw, 3., '#F1F4F5')
    c.rect(x+casew-rolesw, groupy, rolesw, 3., '#E4E8EB')
    positions = []
    cursor = x
    for field in plan['fields']:
        cw = plan['widths'][field]
        positions.append((field, cursor, cw))
        cursor += cw
    label, center = plan['tool_header']
    t(x+center, groupy+1.5, label, 8.5, bold=True, ha='center')
    t(x+casew-rolesw/2, groupy+1.5, 'DNMB components', 8.5, bold=True, ha='center')
    for field, left, cw in positions:
        c.rect(left, groupy+3., cw, 4.8, '#ECEFF1' if field in ('R', 'M', 'S', 'A') else '#F4F7F6')
        t(left+cw/2, groupy+5.4, alias(field, config), 8.6, bold=True, ha='center')
    for i, row in enumerate(r for r in plan['rows'] if r['sample_id'] == sid):
        yy = groupy+7.8+i*5.6
        total = str(row.get('rm_type')).lower() in ('all', 'total')
        for field, left, cw in positions:
            c.rect(left, yy, cw, 5.6, '#DFE4E8' if total else '#F1F3F4' if field in ('R', 'M', 'S', 'A') else '#FAFBFB' if i%2 else 'white')
            cell = _type_cell(row, field)
            text, value, state, subset = cell['display'], cell['value'], cell['state'], cell['subset']
            if field == 'Type' and text.startswith('Type '):
                text = text[5:]
            if status_kind(state) == 'not_applicable':
                c.rect(left, yy, cw, 5.6, '#CCD3D8')
            t(left+cw/2 if field == 'Type' else left+cw-1.8, yy+2.8, text,
                bold=total, ha='center' if field == 'Type' else 'right', owner=[left, yy, cw, 5.6],
                role='type_label' if field == 'Type' else 'type_numeric',
                ref=dict(sample_id=sid, type=row['rm_type'], field=field, value=value, included_subset=subset))
            c.ledger['cells'].append(dict(kind='type_cell', sample_id=sid, rm_type=row['rm_type'],
                field=field, value=value, state=state, included_subset=subset,
                subset_assessment=cell['subset_assessment'], display=text, external_case=external,
                bounds_mm=[left, yy, cw, 5.6], numeric_font_pt=10.8))
        c.line(x, yy+5.6, x+casew, yy+5.6)


def _resize_source(c, width, height):
    c.width, c.height = width, height
    c.fig.set_size_inches(width/25.4, height/25.4)
    c.ax.set_xlim(0, width)
    c.ax.set_ylim(height, 0)
    c.ledger['size_mm'] = [width, height]


def _ink_boxes(c):
    """Measured text ink plus every drawn patch/line extent, in native mm."""
    boxes = [t['bounds_mm'] for t in c.measure()]
    renderer = c.fig.canvas.get_renderer()
    inverse = c.ax.transData.inverted()
    for artist in list(c.ax.patches)+list(c.ax.lines):
        if not artist.get_visible():
            continue
        if isinstance(artist, Patch) and to_rgba(artist.get_facecolor())[3] == 0 and (
                to_rgba(artist.get_edgecolor())[3] == 0 or not artist.get_linewidth()):
            continue
        box = artist.get_window_extent(renderer)
        (x0, y0), (x1, y1) = inverse.transform((box.x0, box.y0)), inverse.transform((box.x1, box.y1))
        boxes.append([min(x0, x1), min(y0, y1), abs(x1-x0), abs(y1-y0)])
    return boxes


def _bottom(boxes, left, right):
    return max((b[1]+b[3] for b in boxes if b[0] < right and b[0]+b[2] > left), default=0.)


def _merged_layout(c, plan, fonts):
    """Place panel c from measured ink in its own horizontal ranges.

    The heading clears ink above its own text span and the tables clear ink
    above their own span, so a tall legend in a different column does not
    force unused vertical whitespace. Clearances are fixed physical spacings.
    """
    cases = plan['cases']
    tables_width = len(cases)*plan['case_width']+(len(cases)-1)*TABLE_GAP
    width = max(c.width, 2*MARGIN_X+tables_width)
    legend = _legend(fonts, plan, cases, width-2*MARGIN_X)
    boxes = _ink_boxes(c)
    heading_right = MARGIN_X+_heading_extent(fonts, legend, TYPE_TITLE)
    offset = _table_offset(legend)
    y = max(_bottom(boxes, MARGIN_X-7, heading_right)+HEADING_CLEARANCE,
            _bottom(boxes, MARGIN_X, MARGIN_X+tables_width)+TABLE_CLEARANCE-offset)
    top = y+offset
    # A legend outside panel c's horizontal span still belongs to the page.
    # Compact placement must never shrink the canvas past existing ink.
    height = max(top+_row_height(plan, cases), max((b[1]+b[3] for b in boxes), default=0.))+MARGIN_BOTTOM
    return dict(y=y, top=top, width=width, height=height, legend=legend)


def letter_export(c, minimum_scale):
    """Uniform vector export; retain native geometry and measured final ink."""
    c.validate()
    width, height = c.width, c.height
    scale = min(LETTER[0]/width, LETTER[1]/height, 1.)
    if scale < minimum_scale-1e-9:
        raise ValueError('Letter composition exceeds its minimum export scale; reduce rows/families per page or select separate pages')
    source = deepcopy(c.ledger)
    ox, oy = (LETTER[0]-width*scale)/2, (LETTER[1]-height*scale)/2
    def point(v):
        return [ox+v[0]*scale, oy+v[1]*scale]
    def box(v):
        return [*point(v), v[2]*scale, v[3]*scale]
    # Text point sizes are physical; geometry continues to use native axes.
    # Paths and data-coordinate circles scale automatically with the axes.
    for artist, meta in c.text_artists:
        artist.set_fontsize(artist.get_fontsize()*scale)
        for key in ('anchor_mm', 'baseline_mm'):
            if key in meta: meta[key] = point(meta[key])
        if meta.get('owner_bounds_mm'): meta['owner_bounds_mm'] = box(meta['owner_bounds_mm'])
        meta['font_pt'] *= scale
    for line in c.ax.lines:
        line.set_markersize(line.get_markersize()*scale)
    for artist in c.ax.texts:
        arrow = getattr(artist, 'arrow_patch', None)
        if arrow is not None:
            arrow.set_mutation_scale(arrow.get_mutation_scale()*scale)
    # Scale every visible stroke, but never below a printable physical width.
    raised = 0
    for artist, native in list(stroke_artists(c.ax)):
        scaled = native*scale
        if scaled < STROKE_FLOOR_PT:
            scaled = STROKE_FLOOR_PT
            raised += 1
        artist.set_linewidth(scaled)
    c.fig.set_size_inches(LETTER[0]/25.4, LETTER[1]/25.4)
    c.ax.set_position([ox/LETTER[0], oy/LETTER[1], width*scale/LETTER[0], height*scale/LETTER[1]])
    c.width, c.height = LETTER
    c.physical_export = True
    # Preserve all native specialized geometry in source_layout. These standard
    # records, used for count and area QA, now use physical export millimeters.
    for kind in ('bars', 'glyphs', 'cells', 'bands'):
        for item in c.ledger[kind]:
            if 'bounds_mm' in item: item['bounds_mm'] = box(item['bounds_mm'])
            for key in ('outer_radius_mm', 'inner_radius_mm', 'filled_width_mm', 'track_width_mm'):
                if key in item: item[key] *= scale
            if 'numeric_font_pt' in item: item['numeric_font_pt'] *= scale
    for symbol in c.ledger.get('rm_symbols', []):
        symbol['center_mm'] = point(symbol['center_mm'])
        symbol['radius_mm'] *= scale
    for key in ('panel_bounds', 'fraction_layout', 'taxonomy_legend_layout', 'section_headings', 'manuscript_workflow'):
        c.ledger.pop(key, None)
    c.ledger.update(size_mm=list(LETTER), source_layout=source,
        coordinate_spaces={'page_geometry':'physical_export_mm', 'source_layout':'native_source_mm'},
        export_transform=dict(scale=scale, offset_mm=[ox, oy], source_size_mm=[width, height],
                              physical_size_mm=list(LETTER), operation='uniform_vector_export',
                              stroke_floor_pt=STROKE_FLOOR_PT, strokes_raised_to_floor=raised))
    return c


def _grid_pages(cases, plan, ncol, nrow, fonts):
    """Row-major pages and their native extents for one candidate grid."""
    pages = []
    for page_cases in _chunks(cases, ncol*nrow):
        rows = _chunks(page_cases, ncol)
        used_cols = max(len(r) for r in rows)
        grid_w = used_cols*plan['case_width']+(used_cols-1)*TABLE_GAP
        content_w = 2*MARGIN_X+grid_w
        legend = _legend(fonts, plan, page_cases, grid_w)
        offset = _table_offset(legend)
        heights = [_row_height(plan, r) for r in rows]
        content_h = MARGIN_TOP+offset+sum(heights)+TABLE_GAP*(len(rows)-1)+MARGIN_BOTTOM
        scale = min(LETTER[0]/content_w, LETTER[1]/content_h, 1.)
        pages.append(dict(cases=page_cases, rows=rows, heights=heights, legend=legend, offset=offset,
                          content=[content_w, content_h], grid_width=grid_w, scale=scale))
    return pages


def _choose_grid(cases, plan, fonts, minimum_scale):
    candidates = []
    for ncol in range(1, min(3, len(cases))+1):
        for nrow in range(1, -(-len(cases)//ncol)+1):
            pages = _grid_pages(cases, plan, ncol, nrow, fonts)
            candidates.append((min(p['scale'] for p in pages), len(pages), ncol, nrow, pages))
    feasible = [c for c in candidates if c[0] >= minimum_scale-1e-9]
    if not feasible:
        raise ValueError('Letter composition exceeds its minimum export scale; reduce rows/families per page or select separate pages')
    preferred = [c for c in feasible if c[0] >= CONTINUATION_TARGET_SCALE]
    if preferred:
        best = min(preferred, key=lambda c: (c[1], -c[0], -c[2]))
    else:
        best = max(feasible, key=lambda c: (c[0], -c[1], c[2]))
    return best


def _continuation_pages(plans, fonts, config, type_config, factory, continued):
    """Remaining type tables, grouped by source columns, on their own Letter pages."""
    groups = []
    for plan in plans:
        key = tuple(plan['fields'])
        if groups and groups[-1][0] == key:
            merged = groups[-1][1]
            merged['cases'] = merged['cases']+plan['cases']
            for name in ('labels', 'styled', 'strain'):
                merged[name] = {**merged[name], **plan[name]}
        else:
            groups.append((key, dict(plan)))
    pages = []
    for _, plan in groups:
        scale, _, ncol, nrow, layout = _choose_grid(plan['cases'], plan, fonts, type_config['composite_min_scale'])
        for page in layout:
            s = page['scale']
            width, height = LETTER[0]/s, LETTER[1]/s
            c = factory(width, height, 'type_benchmark')
            x0 = (width-page['content'][0])/2+MARGIN_X
            y = MARGIN_TOP
            _draw_heading(c, x0, y, page['legend'], plan, page['cases'], continued)
            top = y+page['offset']
            for row, row_height in zip(page['rows'], page['heights']):
                for col, sid in enumerate(row):
                    _draw_case(c, plan, sid, x0+col*(plan['case_width']+TABLE_GAP), top, type_config,
                               _row_strain(plan, row))
                top += row_height+TABLE_GAP
            c.ledger['type_grid'] = dict(columns=ncol, rows_per_page=nrow, rows=len(page['rows']),
                                         content_mm=page['content'], scale_choice=scale,
                                         target_scale=CONTINUATION_TARGET_SCALE)
            pages.append(letter_export(c, type_config['composite_min_scale']))
            continued = True
    return pages


def composite_pages(result, type_result, fonts, config, type_config, factory):
    pages = []
    atlas = []
    if 'atlas' in config['sections']:
        # Native width is the publication template. Extra families/rows continue;
        # the figure is never expanded to make the final text arbitrarily small.
        atlas_config = {**config, 'atlas_page_size': [447., min(config['atlas_page_size'][1], LETTER[1]/config['composite_min_scale'])]}
        atlas = manuscript_atlas_pages(result, fonts, atlas_config, factory)
    plans = _type_plans(type_result, fonts, type_config) if 'type_benchmarks' in config['sections'] else []
    started = False
    if atlas and plans:
        last, plan = atlas[-1], plans[0]
        placed = _merged_layout(last, plan, fonts)
        if min(LETTER[0]/placed['width'], LETTER[1]/placed['height']) >= config['composite_min_scale']:
            _resize_source(last, placed['width'], placed['height'])
            _draw_heading(last, MARGIN_X, placed['y'], placed['legend'], plan, plan['cases'], False)
            for i, sid in enumerate(plan['cases']):
                _draw_case(last, plan, sid, MARGIN_X+i*(plan['case_width']+TABLE_GAP), placed['top'], type_config,
                           _row_strain(plan, plan['cases']))
            last.ledger['section'] = 'atlas_composite'
            last.ledger['type_placement'] = dict(title_anchor_mm=placed['y'], tables_top_mm=placed['top'],
                heading_clearance_mm=HEADING_CLEARANCE, table_clearance_mm=TABLE_CLEARANCE,
                bottom_margin_mm=MARGIN_BOTTOM, method='measured_ink_by_horizontal_span')
            plans.pop(0)
            started = True
    pages.extend(letter_export(c, config['composite_min_scale']) for c in atlas)
    if plans:
        pages.extend(_continuation_pages(plans, fonts, config, type_config, factory, started))
    if 'comparison' in config['sections']:
        pages.extend(comparison_pages(result, fonts, config, factory))
    return pages
