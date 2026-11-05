import hashlib

import pytest

from defenseviz.io.project import Reader


def workbook(tmp_path, formula, *, header_formula=False):
    openpyxl = pytest.importorskip('openpyxl')
    w = openpyxl.Workbook()
    s = w.active
    s.title = 'genes'
    s.append(['locus_tag', 'product', 'link', 'score'])
    s.append(['g1', 'native protein', formula, 4.5])
    if header_formula:
        s['C1'] = formula
    path = tmp_path / 'genes.xlsx'
    w.save(path)
    return path


def test_literal_hyperlink_opt_in_preserves_cells_and_provenance(tmp_path):
    formula = '==HYPERLINK("https://example.org/x?a=1&b=2","A ""quoted"" label α")'
    path = workbook(tmp_path, formula)
    reader = Reader(tmp_path / 'project.json')
    rows, headers = reader.table(path, sheet='genes', literal_hyperlinks=True)
    assert rows == [{'locus_tag': 'g1', 'product': 'native protein',
                     'link': 'A "quoted" label α', 'score': 4.5}]
    assert headers == {'locus_tag', 'product', 'link', 'score'}
    assert reader.table_transforms == [dict(
        path='genes.xlsx', sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        sheet='genes', cell='C2', formula=formula,
        url='https://example.org/x?a=1&b=2', label='A "quoted" label α',
        operation='literal_hyperlink_label', formula_evaluated=False)]
    reader.table(path, sheet='genes', literal_hyperlinks=True)
    assert len(reader.table_transforms) == 1


def test_hyperlinks_are_rejected_without_explicit_opt_in(tmp_path):
    path = workbook(tmp_path, '=HYPERLINK("https://example.org","label")')
    with pytest.raises(ValueError, match='formulas'):
        Reader(tmp_path / 'project.json').table(path, sheet='genes')


@pytest.mark.parametrize('formula', [
    '=1+1', '=SUM(A2:B2)', '=HYPERLINK(A2,"label")',
    '=HYPERLINK("https://example.org",A2)',
    '=HYPERLINK("file:///etc/passwd","label")',
    '=HYPERLINK("javascript:alert(1)","label")',
    '=HYPERLINK("https://example.org"&A2,"label")',
    '=HYPERLINK("https://example.org","label")&"suffix"',
])
def test_opt_in_never_evaluates_nonliteral_formulas(tmp_path, formula):
    path = workbook(tmp_path, formula)
    with pytest.raises(ValueError, match='formulas'):
        Reader(tmp_path / 'project.json').table(path, sheet='genes', literal_hyperlinks=True)


def test_formula_header_is_never_reinterpreted(tmp_path):
    path = workbook(tmp_path, '=HYPERLINK("https://example.org","label")', header_formula=True)
    with pytest.raises(ValueError, match='formulas'):
        Reader(tmp_path / 'project.json').table(path, sheet='genes', literal_hyperlinks=True)
