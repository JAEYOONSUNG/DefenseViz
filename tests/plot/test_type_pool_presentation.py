import json

import pytest

from defenseviz.plot.publication import render_publication
from tests.plot.test_publication import fixture


def test_explicit_case_order_and_supplied_taxon_typography(tmp_path):
    data=fixture(samples=2,families=1)
    for sample in data['samples']:
        sample['display']={'taxon_label':'Synthetic species','strain_label':sample['label'].split()[-2]+' '+sample['label'].split()[-1]}
    render_publication(data,tmp_path,sections=['type_benchmarks'],type_samples=['sample/1','sample/0'])
    page=json.loads((tmp_path/'publication_layout.json').read_text())['pages'][0]
    assert page['sample_ids']==['sample/1','sample/0']
    taxa=[t for t in page['texts'] if t['text']=='Synthetic species']
    assert len(taxa)==2 and all(t['italic'] and t['bold'] and t['font_pt']==9 for t in taxa)
    strains=[t for t in page['texts'] if t['text'].startswith('isolate ')]
    assert len(strains)==2 and all(not t['italic'] and t['font_pt']==9 for t in strains)
    assert not page['validation']['errors']


def test_repeated_case_ids_are_not_silently_duplicated(tmp_path):
    with pytest.raises(ValueError,match='distinct sample'):
        render_publication(fixture(),tmp_path,type_samples=['sample/0','sample/0'])
