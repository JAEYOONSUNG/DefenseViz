"""Reference metadata cannot acquire R evidence through legacy name fallback."""
from pathlib import Path
import json
import shutil
import subprocess

import pytest


HELPER = Path(__file__).parents[1] / "src/defenseviz/data/native/native_metadata.R"
RSCRIPT = shutil.which("Rscript")


def run_r(tmp_path, body):
    if not RSCRIPT:
        pytest.skip("Rscript is required for native metadata fixtures")
    script = tmp_path / "check.R"
    script.write_text(
        f"source({json.dumps(str(HELPER))})\n"
        "make_ref <- function(names,types,roles=rep('R',length(names)),family=rep('Type_II',length(names))) "
        "data.frame(enzyme_name=names,enz_type=types,subunit=roles,rm_type=family,"
        "sequence=rep(paste0(strrep('A',51),'<>'),length(names)),seq_length=rep(53L,length(names)),"
        "rec_seq=rep('AAAA',length(names)),is_gold_standard=rep(FALSE,length(names)))\n"
        + body
    )
    result = subprocess.run([RSCRIPT, "--vanilla", str(script)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("token,expected", [
    ("I", "Type_I"), ("II", "Type_II"), ("III", "Type_III"), ("IV", "Type_IV"),
    ("IIA", "Type_II"), ("IIG", "Type_II"), ("IIIB", "Type_III"),
    ("IX", None), ("IIIX", None), ("V", None), ("I/II", None),
])
def test_exact_type_token(tmp_path, token, expected):
    assertion = "is.na(p$rm_type)" if expected is None else f"identical(p$rm_type,{json.dumps(expected)})"
    run_r(tmp_path, f"p <- .dnmb_project_reference(make_ref('M1.Example',{json.dumps('Type '+token+' methyltransferase')}))$reference\nstopifnot({assertion})")


def test_explicit_numbered_prefix_corrects_stale_valid_role(tmp_path):
    run_r(tmp_path, """
ref <- make_ref(c('M1.X','M2.X','R1.X','S2.X','Unknown'),
                rep('Unspecified',5),c('R','RM','M','R','unknown'))
z <- .dnmb_project_reference(ref)
stopifnot(identical(z$reference$subunit,c('M','M','R','S',NA_character_)),
          all(z$reference$seq_length==51),all(z$audit$previous_sequence_length==53),
          all(z$audit$terminal_marker_removed))
""")


def test_fused_description_wins_but_conflict_is_not_default_r(tmp_path):
    run_r(tmp_path, """
ref <- make_ref(c('M.Fused','M1.Conflict','Unspecified'),
 c('Type IIG restriction enzyme/methyltransferase','Type II restriction enzyme','methyltransferase and specificity'),
 c('M','R','R'))
z <- .dnmb_project_reference(ref)
stopifnot(identical(z$reference$subunit,c('RM',NA_character_,NA_character_)),
 identical(z$audit$role_conflict,c(FALSE,TRUE,TRUE)))
""")


def test_fasta_identity_uses_writer_raw_filter_and_suffixes(tmp_path):
    run_r(tmp_path, """
ref <- make_ref(c('same','same','same'),rep('Type II methyltransferase',3))
ref$sequence <- c(paste0(strrep('A',25),strrep(' ',30)),paste0(strrep('A',50),'<>'),strrep('M',51))
f <- tempfile();writeLines(c('>same_1',ref$sequence[1],'>same_2',ref$sequence[2],'>same_3',ref$sequence[3]),f)
z <- .dnmb_project_reference(ref)
stopifnot(identical(.dnmb_reference_index(z$reference,f),setNames(1:3,paste0('same_',1:3))),z$reference$seq_length[2]==50L)
writeLines(c('>same_1',ref$sequence[1],'>same_2',ref$sequence[2],'>same_3',strrep('R',51)),f)
stopifnot(inherits(try(.dnmb_reference_index(z$reference,f),silent=TRUE),'try-error'))
writeLines(c('>same_2',ref$sequence[1],'>same_1',ref$sequence[2],'>same_3',ref$sequence[3]),f)
stopifnot(inherits(try(.dnmb_reference_index(z$reference,f),silent=TRUE),'try-error'))
""")


def test_supplemental_numbered_m_and_exact_duplicate_metadata(tmp_path):
    run_r(tmp_path, """
ref <- make_ref(c('M1.Example','Duplicate','Duplicate'),c('Type III methyltransferase','Type II restriction enzyme','Type III methyltransferase'))
ref$is_gold_standard <- c(FALSE,FALSE,TRUE)
ref <- .dnmb_project_reference(ref)$reference
idx <- setNames(1:3,c('M1.Example','Duplicate_2','Duplicate_3'))
b <- data.frame(query_id=c('M','D','D'),rebase_enzyme=c('M1.Example','Duplicate_2','Duplicate_3'),pct_identity=c(.9,.8,.8),bitscore=c(100,90,80),evalue=1e-10,length=51)
z <- .dnmb_exact_best_blast(b,ref,idx)
stopifnot(z$subunit[z$query_id=='M']=='M',z$rm_type[z$query_id=='M']=='Type III',
 z$best_rebase_match[z$query_id=='D']=='Duplicate_3',z$subunit[z$query_id=='D']=='M')
b$rebase_enzyme[1]<-'missing'
stopifnot(inherits(try(.dnmb_exact_best_blast(b,ref,idx),silent=TRUE),'try-error'))
""")


def test_context_search_closes_new_neighbors_after_initial_search(tmp_path):
    run_r(tmp_path, """
searched<-c('initial'); calls<-list()
z<-.dnmb_close_context_search(data.frame(query=c('initial','neighbor1')),
 function(x)x$query,function()searched,
 function(ids,round){calls[[round]]<<-ids;searched<<-c(searched,ids)},
 function(x)if('neighbor2' %in% x$query)x else rbind(x,data.frame(query='neighbor2')),3L)
stopifnot(z$inventory$completed,z$inventory$rounds==2L,z$inventory$translated_candidates==3L,
          identical(calls,list('neighbor1','neighbor2')),length(z$inventory$unsearched_candidates)==0L)
""")


def test_context_search_missing_completion_fails_and_empty_is_valid(tmp_path):
    run_r(tmp_path, """
candidate_ids<-function(x)x$query
z<-.dnmb_close_context_search(data.frame(query=character()),candidate_ids,function()character(),
 function(...)stop('unexpected search'),identity,1L)
stopifnot(z$inventory$completed,z$inventory$rounds==0L,z$inventory$translated_candidates==0L)
bad<-try(.dnmb_close_context_search(data.frame(query='not_done'),candidate_ids,function()character(),
 function(...)NULL,identity,1L),silent=TRUE)
stopifnot(inherits(bad,'try-error'))
""")
