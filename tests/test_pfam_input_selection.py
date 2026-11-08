"""Exercise the actual R adapter selector without running native detectors."""
from pathlib import Path
import shutil
import subprocess

import pytest


ADAPTER = Path(__file__).parents[1] / "src/defenseviz/data/native/full_rm.R"
RSCRIPT = shutil.which("Rscript")


def run_selector(tmp_path, body):
    if not RSCRIPT:
        pytest.skip("Rscript is required for the native input-selector tests")
    script = tmp_path / "selector.R"
    script.write_text('''
expressions <- parse(commandArgs(TRUE)[1])
definition <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name(".defenseviz_select_pfam_input")), expressions)
stopifnot(length(definition)==1L)
environment <- new.env(parent=baseenv())
eval(definition[[1]], envir=environment)
select <- environment$.defenseviz_select_pfam_input
make <- function(columns, values=rep("Domain_name", length(columns))) {
  as.data.frame(setNames(as.list(values), columns), check.names=FALSE)
}
''' + body)
    result = subprocess.run([RSCRIPT, "--vanilla", str(script), str(ADAPTER)],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("columns,expected,source", [
    ('"GO.annotations_Pfam", "EggNOG_PFAMs", "Signature.accession_Pfam"', "EggNOG_PFAMs", "eggnog_pfam"),
    ('"Pathways.annotations_Pfam", "Signature.accession_Pfam", "EggNOG_PFAMs"', "Signature.accession_Pfam", "interpro_pfam_accession"),
    ('"Signature description_Pfam", "Signature accession_Pfam"', "Signature accession_Pfam", "interpro_pfam_accession"),
    ('"EggNOG_PFAMs", "Signature.accession_Pfam", "PFAMs"', "PFAMs", "legacy_pfam"),
])
def test_supported_schema_order_and_legacy_preference(tmp_path, columns, expected, source):
    run_selector(tmp_path, f'''
x <- select(make(c({columns})))
stopifnot(identical(x$column, "{expected}"), identical(x$source, "{source}"),
          x$selected_nonplaceholder_cells==1L)
''')


def test_decoy_content_cannot_select_an_annotation_column(tmp_path):
    run_selector(tmp_path, '''
x <- select(make(c("GO.annotations_Pfam", "Pathways.annotations_Pfam",
                   "Signature.description_Pfam", "Score_Pfam", "EggNOG_PFAMs"),
                 c("PF12345", "N6_N4_Mtase", "PF54321", "PF67890", "Domain_name")))
stopifnot(x$column=="EggNOG_PFAMs", length(x$ignored_pfam_columns)==4L,
          identical(x$supported_columns_in_input_order, "EggNOG_PFAMs"))
''')


def test_empty_selected_source_does_not_switch_to_populated_alternative(tmp_path):
    run_selector(tmp_path, '''
x <- select(make(c("EggNOG_PFAMs", "Signature.accession_Pfam"), c("-", "PF12345")))
stopifnot(x$column=="EggNOG_PFAMs", x$selected_nonplaceholder_cells==0L)
d <- data.frame(EggNOG_PFAMs=c(NA, "", "NA", "-", "None", "  ", "Domain_name"), check.names=FALSE)
x <- select(d)
stopifnot(x$selected_nonplaceholder_cells==1L)
''')


def test_no_supported_source_returns_explicit_na_for_native_keyword_path(tmp_path):
    run_selector(tmp_path, '''
for (columns in list(c("product"), c("GO.annotations_Pfam", "custom_pfam_hits"))) {
  x <- select(make(columns))
  stopifnot(is.na(x$column), is.na(x$source), length(x$supported_columns_in_input_order)==0L,
            x$selected_nonplaceholder_cells==0L)
}
''')


@pytest.mark.parametrize("columns", [
    '"Signature accession_Pfam", "Signature.accession_Pfam"',
    '"EggNOG_PFAMs", "EggNOG_PFAMs"',
    '"PFAMs", "PFAMs.1"',
    '"EggNOG_PFAMs", "Signature.accession_Pfam.2"',
])
def test_ambiguous_aliases_and_repaired_duplicate_names_fail(tmp_path, columns):
    run_selector(tmp_path, f'''
error <- tryCatch(select(make(c({columns}))), error=identity)
stopifnot(inherits(error, "error"), grepl("Ambiguous", conditionMessage(error)))
''')
