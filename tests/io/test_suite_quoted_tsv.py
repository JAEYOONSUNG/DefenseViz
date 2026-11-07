"""Quoted legacy tables and literal native annotations have different contracts."""
import csv
import hashlib
import io

import pytest

from defenseviz.io.suite import _native_tsv_table, load_suite
from tests.io.test_suite import native


def parse(text):
    return _native_tsv_table(text, 'synthetic.tsv')


def test_literal_quotes_are_not_stripped_or_repaired():
    rows, _, meta = parse('query\tnote\ng:a\t"literal"\ng:b\tgene "A"\ng:c\tunclosed "\ng:d\t""two""\n')
    assert [r['note'] for r in rows] == ['"literal"', 'gene "A"', 'unclosed "', '""two""']
    assert meta['dialect'] == 'literal'
    assert meta['selection_reason'] == 'valid_literal_table_preserved'


def test_uniform_quoted_header_decodes_escaped_quotes_and_header_only_zero():
    rows, headers, meta = parse('"query"\t"note"\n"g:a"\t"gene ""A"""\n')
    assert rows == [{'query': 'g:a', 'note': 'gene "A"'}]
    assert headers == {'query', 'note'}
    assert meta['selection_reason'] == 'uniformly_quoted_header'
    assert parse('"query"\t"note"\n')[0] == []
    assert parse('query\tnote\n')[0] == []


@pytest.mark.parametrize('quoting', [csv.QUOTE_ALL, csv.QUOTE_MINIMAL])
def test_encoded_tab_newlines_and_empty_trailing_cells(quoting):
    records = [['g:a', 'line1\n# inside note\n"end"\tvalue', ''], ['g:b', 'plain', '']]
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, delimiter='\t', quoting=quoting, lineterminator='\r\n')
    writer.writerow(['query', 'note', 'tail'])
    writer.writerows(records)
    rows, _, meta = parse('# metadata\r\n\r\n' + stream.getvalue())
    assert [[r['query'], r['note'], r['tail']] for r in rows] == records
    assert meta['dialect'] == 'quoted'
    assert meta['preamble_lines'] == 2
    assert meta['multiline_records'] == [{'row': 2, 'first_line': 4, 'last_line': 6}]


@pytest.mark.parametrize('text', [
    '"query"\t"note"\n"g:a"\t"unterminated\n',
    '"query"\t"note"\n"g:a"junk\t"note"\n',
    '"query"\t"note"\ng:a\tgene"A\n',
])
def test_declared_quoted_errors_never_fall_back_to_literal(text):
    with pytest.raises(ValueError):
        parse(text)


@pytest.mark.parametrize('newline', ['\r', '\n', '\r\n'])
@pytest.mark.parametrize('bad_row, line, reason', [
    ('"g:a"\t"text \\"escaped\\""', 4, 'text after closing field quote'),
    ('"g:a"\tgene"A', 4, 'bare quote'),
    ('"g:a"\t"first{nl}second"trailing', 5, 'text after closing field quote'),
    ('"g:a"\t"unterminated', 4, 'unterminated quoted field'),
    ('"g:a"\t"extra"\t"field"', 4, 'column count'),
])
def test_quoted_diagnostics_identify_source_and_physical_line(newline, bad_row, line, reason):
    path = 'synthetic/조합 input.tsv'
    text = newline.join(['# provenance', '', '"query"\t"note"', bad_row.format(nl=newline)])
    with pytest.raises(ValueError) as error:
        _native_tsv_table(text, path)
    assert path in str(error.value)
    assert f'physical line {line}:' in str(error.value)
    assert reason in str(error.value)


def test_duplicate_quoted_header_diagnostic_includes_preamble_offset():
    with pytest.raises(ValueError, match=r'physical line 3:.*duplicate column header'):
        parse('# metadata\n\n"query"\t"query"\n')


def test_ambiguous_structural_quoting_cannot_create_phantom_genes():
    # Both literal rows and a quoted multiline note have two columns. Choosing
    # the literal parse would invent a gene named '# note'; choosing quoted
    # would silently merge two potentially real literal rows.
    text = 'query\tnote\ng:a\t"line1\n# note\tend"\ng:b\tplain\n'
    with pytest.raises(ValueError, match='ambiguous structural quoting'):
        parse(text)


@pytest.mark.parametrize('newline', ['\r', '\n', '\r\n'])
def test_first_column_ambiguity_detected_with_every_native_newline(newline):
    # A field-opening quote at the beginning of a CR-only record must be
    # detected just like LF/CRLF. Regex MULTILINE '^' alone misses that boundary.
    text = newline.join(['query\tnote', '"g:a\tx', 'hi"\ty', ''])
    with pytest.raises(ValueError, match='ambiguous structural quoting'):
        parse(text)


@pytest.mark.parametrize('text, reason', [
    ('query\tnote\ng:a\tx\textra\n', 'column count'),
    ('query\tnote\ng:a\tunquoted\nnote\n', 'column count'),
    ('query\tquery\ng:a\tg:b\n', 'duplicate column header'),
    ('"query"\t"query"\n', 'duplicate column header'),
    ('query\tnote\ng:a\t\x00\n', 'NUL byte'),
    ('', 'empty file'),
])
def test_invalid_tables_remain_invalid(text, reason):
    with pytest.raises(ValueError, match=reason):
        parse(text)


def test_real_import_preserves_semantic_fields_and_exact_raw_source_hash(tmp_path):
    root = native(tmp_path)
    path = root / 'dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv'
    with path.open(newline='') as stream:
        records = list(csv.reader(stream, delimiter='\t'))
    previous_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    source_before = load_suite(root)
    with path.open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream, delimiter='\t', quoting=csv.QUOTE_ALL)
        writer.writerows(records)
    expected = hashlib.sha256(path.read_bytes()).hexdigest()
    field_limit = csv.field_size_limit()
    after = load_suite(root)
    assert csv.field_size_limit() == field_limit
    assert after['input_sha256'][str(path)] == expected
    assert not after['suite']['unmapped_records']
    def rebind(value):
        if isinstance(value, dict):
            return {k: rebind(v) for k, v in value.items()}
        if isinstance(value, list):
            return [rebind(v) for v in value]
        return expected if value == previous_hash else value
    for key in ['calls', 'saved_decisions']:
        assert after[key] == rebind(source_before[key])
    curation = [r for r in after['curation_records'] if r['source_ref']['path'].endswith(path.name)]
    assert curation and all(r['source_ref']['sha256'] == expected for r in curation)
    transformed = [r for r in after['suite']['table_transforms'] if r['path'].endswith(path.name)]
    assert len(transformed) == 1
    assert transformed[0]['dialect'] == 'quoted'
    assert transformed[0]['sha256'] == expected
    assert transformed[0]['logical_data_rows'] == len(records) - 1


def test_multiline_gene_note_is_parsed_without_workbook_fallback(tmp_path):
    root = native(tmp_path)
    path = root / 'dnmb_module_rebasefinder/rebasefinder_input.tsv'
    with path.open(newline='') as stream:
        records = list(csv.reader(stream, delimiter='\t'))
    records[1][-1] = 'native description\nsecond line "quoted"'
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerows(records)
    data = load_suite(root)
    assert len(data['genes']) == 3
    assert not any(r['code'] == 'malformed_gene_tsv_fallback' for r in data['suite']['conflicts'])
    meta = next(r for r in data['suite']['table_transforms'] if r['path'].endswith(path.name))
    assert meta['dialect'] == 'quoted'
    assert meta['multiline_records'] == [{'row': 2, 'first_line': 2, 'last_line': 3}]
