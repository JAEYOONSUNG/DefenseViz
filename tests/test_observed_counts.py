"""A computed annotation inventory is not a declaration of biological absence."""
from copy import deepcopy

import pytest

from defenseviz.automatic_systems import summarize_assessments


def row():
    return dict(sample_id='s',gene_id='g',in_integrated_universe=True,
                review_execution='completed',endpoints={
                    axis:dict(status='unresolved',coverage='unassessed')
                    for axis in ('effector','associated','motif')})


def test_completed_review_counts_observed_zero_without_resolving_function():
    source=[row()];before=deepcopy(source)
    result=summarize_assessments(source,[['s','g']],inventory_status='complete')
    assert source==before
    assert result['effector_CDS']==0 and result['unresolved_CDS']['effector']==1
    assert result['metric_status']['effector']=='unassessed'
    assert result['observed_metric_status']['effector']=='complete'
    assert result['observed_count_scope']=='qualified_annotations_in_reviewed_observed_integrated_members'


@pytest.mark.parametrize('inventory,execution',[('partial','completed'),('unassessed','completed'),('complete','pending'),('complete',None)])
def test_incomplete_inventory_or_execution_cannot_supply_observed_zero(inventory,execution):
    item=row();item['review_execution']=execution
    result=summarize_assessments([item],[['s','g']],inventory_status=inventory)
    assert result['observed_metric_status']['effector']=='unassessed'


def test_missing_member_is_not_a_zero_count():
    with pytest.raises(ValueError,match='unaccounted'):
        summarize_assessments([row()],[['s','missing']],inventory_status='complete')


def test_successfully_empty_observed_inventory_is_countable():
    result=summarize_assessments([],[],inventory_status='complete')
    assert result['CDS_count']==result['effector_CDS']==0
    assert result['observed_metric_status']['effector']=='complete'
