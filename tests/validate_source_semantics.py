"""Regression guard for source meaning, distinct from arithmetic correctness.

Run after audit_source_semantics.py. This validates provenance and the concrete
counterexamples that reject the old guide-significance claim, while permitting
the separately supported, qualified published-catalogue comparison.
"""
import base64
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import struct
import sys
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
from fetch_inputs import ROTTNER_ARCHIVE, resolve_rottner_archive

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-cache', help='External author-code cache directory; validation never downloads')
    args = parser.parse_args()
    audit_path = ROOT / "evidence/source_semantics_audit.json"
    audit = json.loads(audit_path.read_text())
    assert audit["audit_script_sha256"] == sha(ROOT / "code/audit_source_semantics.py"), "Stale semantic audit: rerun source audit"
    for path, expected in {**audit["source_hashes"], **audit["derived_evidence_hashes"]}.items():
        assert sha(ROOT / path) == expected, f"Semantic provenance changed: {path}"
    checks = []
    def check(name, condition):
        assert condition, name
        checks.append(name)
    archive = resolve_rottner_archive(args.source_cache, offline=True)
    check('remote_author_archive_hash', sha(archive) == audit['remote_source_hashes'][ROTTNER_ARCHIVE['id']] == ROTTNER_ARCHIVE['sha256'])
    author = audit['archived_author_code']
    check('remote_author_source_identity', author['source_url'] == ROTTNER_ARCHIVE['url'] and author['storage'] == 'external_cache_not_redistributed')
    with zipfile.ZipFile(archive) as source:
        check('author_selection_script_hash', hashlib.sha256(source.read(author['selection_script_member'])).hexdigest() == author['selection_script_sha256'])
    distributions = {"p.low": (72618,69562), "p.high": (69562,72618), "p.twosided": (110998,31182), "FDR": (99564,42616)}
    s8 = audit["xml"]["Table S8-Screen_all sgRNA"]["columns"]
    for column, (zeros, ones) in distributions.items():
        item = s8[column]
        check("stored_numeric_binary_"+column, item["stored_values"] == {"0": zeros, "1": ones} and item["cell_types"] == {"n": 142180})
    # Independently decode the actual workbook's saved M bytes, not merely a
    # generated statement that the query was inspected.
    with zipfile.ZipFile(ROOT / audit["source_file"]) as book:
        data = base64.b64decode(ET.fromstring(book.read("customXml/item1.xml")).text)
        version, size = struct.unpack("<II", data[:8])
        assert version == 0
        with zipfile.ZipFile(io.BytesIO(data[8:8+size])) as package:
            mbytes = package.read("Formulas/Section1.m")
    pipeline = audit["saved_import_pipeline"]
    check("saved_M_hash_verified_against_workbook", hashlib.sha256(mbytes).hexdigest() == pipeline["raw_M_sha256"])
    mtext = mbytes.decode("utf-8-sig")
    query_name = '#"RRApaired_screen sgrna_summary (1)"'
    raw_query = mtext.split("shared " + query_name + " = let",1)[1].split(";",1)[0]
    for column in distributions:
        check("saved_S8_integer_conversion_"+column, '{"'+column+'", Int64.Type}' in raw_query and pipeline["queries"][query_name]["column_types"][column] == "Int64.Type")
    types = pipeline["queries"]["sgRNA_hits"]["column_types"]
    check("S7_column_specific_precision_preserved", types["p.low"] == "Int64.Type" and all(types[x] == "type number" for x in ["p.high","p.twosided","FDR"]))
    check("S8_LFC_counts_continuous_types", all(pipeline["queries"][query_name]["column_types"][x] == "type number" for x in ["LFC","control_count","treatment_count"]))
    diag = audit["diagnostics"]
    linked = diag["S7_S8_linkage"]
    check("2384_unambiguous_links_and_rounding_matches", linked["unambiguous_rows"] == 2384 and linked["ambiguous_or_missing"] == [] and linked["FDR_S8_equals_nearest_integer_S7_count"] == 2384)
    with (ROOT / "evidence/S7_S8_semantic_linkage.csv").open(newline="",encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    counterexample = [r for r in rows if r["gene"] == "TOMM6" and r["S7_row"] == "2385" and r["S8_row"] == "91896"]
    check("TOMM6_zero_is_not_q_below_010", len(counterexample) == 1 and float(counterexample[0]["S7_FDR"]) == .42555 and float(counterexample[0]["S8_stored_FDR"]) == 0 and counterexample[0]["replicate"] == "r0" and counterexample[0]["coordinate"] == "chr6:41755542-41755561")
    check("2421_nonpublished_stored_rule_passes", diag["nonpublished_primary_positive_stored_zero_rule_count"] == 2421 and len(set(diag["nonpublished_primary_positive_stored_zero_rule_genes"])) == 2421)
    candidate_genes = {"C12orf29","COA5","CREB3","EXOSC3","NCOA6","PDS5A","PEX13","RAD52","SLC2A11","SLC8B1","SSR2","TM2D2","UBE2C","ZNF616"}
    subset = diag["candidate_S7_selected_subset_evidence"]
    s7 = {r["gene"]: r for r in subset["candidates"]}
    check("selected_S7_membership", set(s7) == candidate_genes | {"NUPL1","TOMM6"})
    for gene in sorted(candidate_genes | {"NUPL1"}):
        r = s7[gene]
        check("selected_S7_two_paired_guides_"+gene, r["source_s7_positive_fdr_lt010_r0"] == r["source_s7_positive_fdr_lt010_r1"] == r["source_s7_positive_fdr_lt010_sameguide"] == 2 and len(r["source_s7_positive_fdr_lt010_sameguide_keys"].split(';')) == 2)
        check("selected_S7_incomplete_"+gene, all(r["S7_rows_"+rep] == 2 and r["S8_rows_"+rep] == 4 and not r["S7_covers_all_S8_rows_"+rep] for rep in ["r0","r1"]))
        check("no_selected_S7_third_guide_claim_"+gene, not r["source_s7_three_positive_each_replicate"] and not r["source_s7_three_identical_positive_guides"])
    check("no_exact_point_one_ties", subset["all_S7_exact_FDR_010_ties"] == 0)
    check("TOMM6_selected_positive_counts", [s7["TOMM6"][k] for k in ["source_s7_positive_fdr_lt010_r0","source_s7_positive_fdr_lt010_r1","source_s7_positive_fdr_lt010_sameguide"]] == [1,2,1])
    prev = diag["opposing_stored_zero_guide_prevalence"]["legacy14_candidates"]
    check("opposing_encoded_signs_remain_visible", prev["n"] == 14 and prev["any_zero_negative_guide_either_replicate"] == 10 and prev["at_least_one_zero_negative_guide_each_replicate"] == 7)
    catalogue = diag["independent_revised_catalogue_checks"]
    check("final_revised_comparisons_not_deferred", catalogue["comparison_status"] == "raw_reconstruction_completed")
    for name, expected in [("primary", [12184,285,479,14]), ("no_expression",[17276,390,514,15])]:
        r = catalogue["analyses"][name]
        check("qualified_catalogue_membership_"+name, [r[x] for x in ["n","A","L","overlap"]] == expected and r["output_membership_identical"] and r["output_A_L_identical"])
    check("all427_axes_match", catalogue["candidate_axes"]["total"] == 427 and catalogue["candidate_axes"]["output_membership_and_axes_identical"])
    permission = audit["claim_permissions"]
    check("V1_guide_significance_claim_REJECTED", permission["V1_S8_FDR_below_010_independent_guide_significance"] == "REJECTED")
    check("genomewide_q_reconstruction_NOT_claimed", permission["complete_genomewide_continuous_q_values_reconstructed"] is False)
    check("selected_subset_scope_only", permission["selected_Table_S7_continuous_values"] == "SOURCE_SUPPORTED_SELECTED_SUBSET_ONLY")
    check("catalogue_permitted_with_qualification", permission["published_Table_S1_catalogue_comparison"] == "PERMITTED_WITH_QUALIFICATION")
    result = {"status": "PASS", "checks_passed": len(checks), "checks": checks,
        "V1_guide_significance_claim": "REJECTED", "complete_genomewide_q_reconstruction": False,
        "catalogue_disposition": "PERMITTED_WITH_QUALIFICATION", "semantic_audit_sha256": sha(audit_path), "test_script_sha256": sha(Path(__file__)),
        "scope": "Regression guard for source interpretation and provenance; no biological validation or recovery of missing probabilities is claimed."}
    (ROOT / "evidence/source_semantics_validation.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps({k:result[k] for k in ["status","checks_passed","V1_guide_significance_claim","complete_genomewide_q_reconstruction","catalogue_disposition"]}))

if __name__ == "__main__": main()
