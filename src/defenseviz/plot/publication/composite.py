"""Optional two-page manuscript composition with measured continuation pages.

The atlas is drawn at its established source size, then uniformly exported to
Letter landscape. This explicit publication export never alters measurements,
source memberships, or the fixed type sizes used by the separate-page layout.
"""
from copy import deepcopy

from .atlas import _chunks, _type_display, _type_case_label_layout, _type_case_label
from .comparison import comparison_pages, source_order
from .drawing import alias, status_kind
from .manuscript_atlas import manuscript_atlas_pages
from .manuscript_text import draw_manuscript_text

LETTER = (279.4, 215.9)


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
    for sources in _chunks(source_order(result, config), 3):
        fields = ['Type'] + (['Pool'] if any('pool_CDS' in r for r in rows) else []) + sources + ['DNMB', 'R', 'M', 'S', 'A']
        base = {'Type': 10., 'Pool': 12., 'DNMB': 15., 'R': 13., 'M': 16., 'S': 14., 'A': 13.3333333333}
        def measured_value(row, field):
            value = _type_display(row, field)[0]
            if field == 'Type' and value.startswith('Type '): value = value[5:]
            return fonts.width(value, 10.8, True) + (1. if field == 'Type' else 3.)
        widths = {field: max(base.get(field, 15.), fonts.width(alias(field, config), 8.6, True)+3,
                     max((measured_value(r, field) for r in rows), default=0)) for field in fields}
        casew = sum(widths.values())
        width_budget = LETTER[0]/config['composite_min_scale']
        perpage = min(3, max(1, int((width_budget-20.+6.)/(casew+6.))))
        for cases in _chunks(ids, perpage):
            labels = {sid: fonts.wrap(samples[sid].get('label', sid), casew-4, 11.4, True) for sid in cases}
            styled = {sid: _type_case_label_layout(fonts, casew, samples[sid], 11.4) for sid in cases}
            strainh = max(5., max(styled[sid]['height'] if styled[sid] is not None else fonts.height(labels[sid], 11.4, True)+2 for sid in cases))
            rowcount = max(sum(r['sample_id'] == sid for r in rows) for sid in cases)
            plans.append(dict(cases=cases, fields=fields, widths=widths, case_width=casew,
                width=max(447., 20.+len(cases)*casew+(len(cases)-1)*6.),
                height=8.+strainh+3.+4.8+rowcount*5.6+2., strain_height=strainh,
                labels=labels, styled=styled, samples=samples, rows=rows, selection=selection))
    return plans


def _draw_types(c, plan, y, config):
    """Append table rows only from the validated independently scoped results."""
    t = lambda x, yy, value, pt=10.8, **kw: draw_manuscript_text(c, x, yy, value, pt, **kw)
    t(3, y, 'c', 15, bold=True)
    t(10, y, 'CDS inclusion by assigned R-M type', 11.5, bold=True)
    t(10, y+3.8, 'R(n): incomplete context · M/S(n): included orphans · —: unavailable or not applicable', 9.4)
    top = y+8.
    casew = plan['case_width']
    strainh = plan['strain_height']
    c.ledger['type_sample_ids'] = plan['cases']
    c.ledger['type_sample_selection'] = plan['selection']
    c.ledger['type_source_order'] = [f for f in plan['fields'] if f not in ('Type', 'Pool', 'DNMB', 'R', 'M', 'S', 'A')]
    for ci, sid in enumerate(plan['cases']):
        sample = plan['samples'][sid]
        x = 10.+ci*(casew+6.)
        c.rect(x, top, casew, strainh, '#F3EFE4' if sample.get('cohort') in ('external', 'external_validation') else '#F8E9EF')
        _type_case_label(c, x, top, casew, strainh, sample, plan['labels'][sid], plan['styled'][sid], 11.4)
        groupy = top+strainh
        rolesw = sum(plan['widths'][f] for f in ('R', 'M', 'S', 'A'))
        c.rect(x, groupy, casew-rolesw, 3., '#F1F4F5')
        c.rect(x+casew-rolesw, groupy, rolesw, 3., '#E4E8EB')
        t(x+(casew-rolesw)/2, groupy+1.5, 'Shared-pool CDSs', 8.5, bold=True, ha='center')
        t(x+casew-rolesw/2, groupy+1.5, 'DNMB components', 8.5, bold=True, ha='center')
        positions = []
        cursor = x
        for field in plan['fields']:
            cw = plan['widths'][field]
            positions.append((field, cursor, cw))
            c.rect(cursor, groupy+3., cw, 4.8, '#ECEFF1' if field in ('R', 'M', 'S', 'A') else '#F4F7F6')
            t(cursor+cw/2, groupy+5.4, alias(field, config), 8.6, bold=True, ha='center')
            cursor += cw
        for i, row in enumerate(r for r in plan['rows'] if r['sample_id'] == sid):
            yy = groupy+7.8+i*5.6
            total = str(row.get('rm_type')).lower() in ('all', 'total')
            for field, left, cw in positions:
                c.rect(left, yy, cw, 5.6, '#DFE4E8' if total else '#F1F3F4' if field in ('R', 'M', 'S', 'A') else '#FAFBFB' if i%2 else 'white')
                text, value, state, subset = _type_display(row, field)
                if field == 'Type' and text.startswith('Type '):
                    text = text[5:]
                if status_kind(state) == 'not_applicable':
                    c.rect(left, yy, cw, 5.6, '#CCD3D8')
                t(left+cw/2 if field == 'Type' else left+cw-1.8, yy+2.8, text,
                    bold=total, ha='center' if field == 'Type' else 'right', owner=[left, yy, cw, 5.6],
                    role='type_label' if field == 'Type' else 'type_numeric',
                    ref=dict(sample_id=sid, type=row['rm_type'], field=field, value=value, included_subset=subset))
                c.ledger['cells'].append(dict(kind='type_cell', sample_id=sid, rm_type=row['rm_type'],
                    field=field, value=value, state=state, included_subset=subset, display=text,
                    bounds_mm=[left, yy, cw, 5.6], numeric_font_pt=10.8))
            c.line(x, yy+5.6, x+casew, yy+5.6)
    c.ledger['type_context_note'] = 'Counts are annotation assignments. Context does not establish biological activity.'


def _resize_source(c, width, height):
    c.width, c.height = width, height
    c.fig.set_size_inches(width/25.4, height/25.4)
    c.ax.set_xlim(0, width)
    c.ax.set_ylim(height, 0)
    c.ledger['size_mm'] = [width, height]


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
    for patch in c.ax.patches:
        patch.set_linewidth(patch.get_linewidth()*scale)
    for line in c.ax.lines:
        line.set_linewidth(line.get_linewidth()*scale)
        line.set_markersize(line.get_markersize()*scale)
    for artist in c.ax.texts:
        arrow = getattr(artist, 'arrow_patch', None)
        if arrow is not None:
            arrow.set_mutation_scale(arrow.get_mutation_scale()*scale)
            arrow.set_linewidth(arrow.get_linewidth()*scale)
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
                              physical_size_mm=list(LETTER), operation='uniform_vector_export'))
    return c


def composite_pages(result, type_result, fonts, config, type_config, factory):
    pages = []
    atlas = []
    if 'atlas' in config['sections']:
        # Native width is the publication template. Extra families/rows continue;
        # the figure is never expanded to make the final text arbitrarily small.
        atlas_config = {**config, 'atlas_page_size': [447., min(config['atlas_page_size'][1], LETTER[1]/config['composite_min_scale'])]}
        atlas = manuscript_atlas_pages(result, fonts, atlas_config, factory)
    plans = _type_plans(type_result, fonts, type_config) if 'type_benchmarks' in config['sections'] else []
    if atlas and plans:
        last, plan = atlas[-1], plans[0]
        text = last.measure()
        inkbottom = max((t['bounds_mm'][1]+t['bounds_mm'][3] for t in text), default=last.height)
        # A few mm of white separation is sufficient; reserved trailing whitespace
        # is reusable, but the actual measured taxonomy legend is never cropped.
        y = inkbottom+4.
        height = y+plan['height']
        width = max(last.width, plan['width'])
        if min(LETTER[0]/width, LETTER[1]/height) >= config['composite_min_scale']:
            _resize_source(last, width, height)
            _draw_types(last, plan, y, type_config)
            last.ledger['section'] = 'atlas_composite'
            plans.pop(0)
    pages.extend(letter_export(c, config['composite_min_scale']) for c in atlas)
    for plan in plans:
        c = factory(plan['width'], plan['height']+10., 'type_benchmark')
        _draw_types(c, plan, 5., type_config)
        pages.append(letter_export(c, config['composite_min_scale']))
    if 'comparison' in config['sections']:
        pages.extend(comparison_pages(result, fonts, config, factory))
    return pages
