"""Bind editable narrative to computed facts and tables; refuse factual drift."""
import copy, hashlib, json, re
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
def load(name):return json.loads((ROOT/name).read_text(encoding='utf-8'))
def table(caption,headers,rows,widths,font=8.5):
    return {'type':'table','caption':caption,'headers':headers,'rows':rows,'widths':widths,'font_size':font}

def require_intervals_include_one(rows,context):
    rows=list(rows)
    if not rows or any(not (r['ci_lower']<=1<=r['ci_upper']) for r in rows):
        raise ValueError(f'{context}: interval excludes 1; revise the unresolved/interval-includes-1 narrative before building')

def prepare():
    dependencies=['manuscript/narrative.json','results/results.json','results/candidate_accounting.json',
      'results/candidate_table_427.csv','results/conditional_information.json','results/continuous_summary.json',
      'results/stratum_influence.csv','results/legacy_protocol/results.json','evidence/source_semantics_audit.json',
      'docs/source_design_comparison.json','docs/candidate_annotations.json','evidence/candidate_S7_selected_guide_evidence.csv','evidence/identifier_mapping_sensitivity.json']
    narrative=load(dependencies[0]);result=load('results/results.json');assert result['version']==2
    analyses={r['name']:r for r in result['analyses']};p=analyses['primary']
    c=pd.read_csv(ROOT/'results/candidate_table_427.csv');a=load('results/candidate_accounting.json')
    info=load('results/conditional_information.json');continuous=load('results/continuous_summary.json')
    influence=pd.read_csv(ROOT/'results/stratum_influence.csv')
    semantic=load('evidence/source_semantics_audit.json');diag=semantic['diagnostics'];link=diag['S7_S8_linkage']
    selected=c[c.primary_eligible&c.published_content_label.eq('Decrease')]
    comparators=c[c.primary_eligible&c.A&c.published_content_label.eq('Increase')]
    assert len(comparators)==2 and set(comparators.gene)=={'RAB14','LRRC61'}, 'Higher-content comparator membership drift'
    assert comparators.screen_evidence_status.eq('Published higher-content evidence').all()
    comparators=comparators.set_index('gene').loc[['RAB14','LRRC61']].reset_index()
    recovered=pd.read_csv(ROOT/'evidence/candidate_S7_selected_guide_evidence.csv').set_index('gene')
    candidate_index=c.set_index('gene')
    recovered_fields=['source_s7_positive_fdr_lt010_r0','source_s7_positive_fdr_lt010_r1','source_s7_positive_fdr_lt010_sameguide']
    for gene in set(candidate_index.index)&set(recovered.index):
        assert all(candidate_index.loc[gene,key]==recovered.loc[gene,key] for key in recovered_fields), ('Candidate selected-guide join drift',gene)
    assert (recovered.loc[list(selected.gene)+['NUPL1'],['source_s7_positive_fdr_lt010_r0','source_s7_positive_fdr_lt010_r1','source_s7_positive_fdr_lt010_sameguide']]==2).all().all()
    assert not recovered.loc[list(selected.gene)+['NUPL1'],['S7_covers_all_S8_rows_r0','S7_covers_all_S8_rows_r1']].any().any()
    assert len(selected)==p['overlap']==len(a['primary_lower_content_candidates'])
    assert len(c)==a['total'] and c.primary_eligible.sum()==a['primary_eligible']==p['A']
    assert int(c.in_U0.sum())==a['common_coverage']
    assert p['overlap']==p['table_2x2'][0][0] and sum(map(sum,p['table_2x2']))==p['n']
    assert info['observed_overlap']==p['overlap'] and info['n_genes']==p['n']
    assert result['primary_membership_and_estimates_unchanged']
    assert link['unambiguous_rows']==link['FDR_S8_equals_nearest_integer_S7_count']
    fdr_values=semantic['xml']['Table S8-Screen_all sgRNA']['columns']['FDR']['stored_values']
    n0=int(fdr_values['0']);n1=int(fdr_values['1'])
    tomm=link['S7_FDR_greater_than_0_1_with_S8_zero'][0]
    np1=c[c.gene.eq('NUPL1')].iloc[0]
    assert np1.screen_evidence_status=='Published lower-content evidence' and not np1.primary_eligible
    probs={x['common_odds_ratio']:x['criterion_probability'] for x in info['criterion_probabilities']}
    mapping=load('evidence/identifier_mapping_sensitivity.json');mapped=mapping['conditional']
    require_intervals_include_one(analyses.values(),'Revised catalog analyses')
    require_intervals_include_one(influence.to_dict('records'),'Leave-one-stratum-out analyses')
    require_intervals_include_one([mapped],'Identifier-mapping sensitivity')
    fmt=lambda r:f"{r['conditional_mle']:.4f} (exact95% CI {r['ci_lower']:.4f}-{r['ci_upper']:.4f})"
    facts={
      'title':narrative['title'],
      'mapping_sentence':f"A separate one-to-one HGNC/Ensembl mapping sensitivity includes {mapping['n']:,} genes and {mapping['A']} exposed targets, with {mapping['overlap']} overlaps: OR {mapped['conditional_mle']:.4f} (exact95% CI {mapped['ci_lower']:.4f}-{mapped['ci_upper']:.4f}), p={mapped['p_enrichment']:.4f}. NUPL1 enters as NUP58, while the SLC24A6/SLC8B1 many-to-one collision is excluded rather than merged. Every mapping and exclusion is recorded; this post-review sensitivity does not alter the exact-symbol primary.",
      'primary_abstract':f"Among {p['n']:,} eligible genes, {p['overlap']} of {p['A']} survival-favoring targets ({100*p['fraction_A1']:.2f}%) had published lower-content calls, versus {p['table_2x2'][1][0]} of {p['n']-p['A']:,} other genes ({100*p['fraction_A0']:.2f}%). The conditional odds ratio was {fmt(p)}, with one-sided enrichment p={p['p_enrichment']:.4f}; enrichment remains unresolved.",
      'accounting_abstract':f"Across all {a['total']} source survival-favoring genes, {a['screen_evidence_status_counts']['Published lower-content evidence']} have published lower-content and {a['screen_evidence_status_counts']['Published higher-content evidence']} higher-content evidence; {a['primary_eligible']} meet primary eligibility. NUPL1 retains lower-content evidence despite missing the exact expression-symbol match.",
      'semantic_finding':f"All {n0+n1:,} S8 FDR values are stored as zero ({n0:,}) or one ({n1:,}), not merely displayed that way. All {link['unambiguous_rows']:,} uniquely linked S7 rows have S8 FDR values equal to nearest-integer rounding of the continuous S7 value. For a matched {tomm['gene']} guide, S7 FDR {tomm['S7_FDR']} becomes S8 zero. The linked rows demonstrate lost FDR precision.",
      'rule_discrepancy':f"Among primary genes absent from Table S1, {diag['nonpublished_primary_positive_stored_zero_rule_count']:,} nevertheless satisfy the former positive-LFC/stored-zero rule in at least two guides per replicate.",
      'nupl1_sentence':f"NUPL1 illustrates the distinction: survival beta {np1.survival_beta:.5f}, survival FDR {np1.survival_fdr:.5f}, a published Decrease call and {int(np1.covered_insulin_guides)} covered guides are available. Its primary exclusion is {np1.eligibility_reason.replace('_',' ')}. The evidence field preserves the lower-content call independently of this adjustment-variable exclusion.",
      'primary_results':f"The published-only primary comparison has {p['overlap']}/{p['A']} lower-content calls among survival-favoring genes versus {p['table_2x2'][1][0]}/{p['n']-p['A']:,} among the remaining eligible genes. The conditional odds ratio is {fmt(p)}; one-sided enrichment p={p['p_enrichment']:.4f}. Removing the invalid stored-FDR audit leaves the original primary membership, strata and numerical estimate unchanged. The unstratified conditional estimate is {p['unstratified']['conditional_mle']:.4f} ({p['unstratified']['ci_lower']:.4f}-{p['unstratified']['ci_upper']:.4f}).",
      'accounting_results':f"Of {a['total']} source survival-favoring genes, {a['common_coverage']} have common coverage and {a['primary_eligible']} enter the primary comparison. Across the full resource, {a['screen_evidence_status_counts']['Published lower-content evidence']} have published lower-content evidence, {a['screen_evidence_status_counts']['Published higher-content evidence']} higher-content evidence, {a['screen_evidence_status_counts']['No published directional call']} have an insulin measurement without a published directional call, and {a['screen_evidence_status_counts']['Screen measurement unavailable']} lack an exact insulin-screen gene measurement. Only {len(selected)} of the {a['screen_evidence_status_counts']['Published lower-content evidence']} lower-content candidates are primary-eligible.",
      'guide_concordance':f"All {len(selected)} primary candidates and NUPL1 have two positive guides with FDR <0.10 in each replicate in the recovered S7 subset, with the same two guide identities paired. Only two of four guides per replicate are represented there; the missing probabilities remain unknown. This supports the available selected-guide evidence, not full reconstruction of hit selection or a third significant guide. Only {a['primary_survival_fdr_lt005_lower_content']} of the {len(selected)} primary candidates meet source survival FDR <0.05. CREB3 and PDS5A have three positive stored-zero guides per replicate but only two identical identities meeting that provisional rule in both. {diag['opposing_stored_zero_guide_prevalence']['legacy14_candidates']['any_zero_negative_guide_either_replicate']} of {len(selected)} have an opposing negative-LFC stored-zero guide in at least one replicate; stored zero does not establish its significance.",
      'sensitivity_sentence':f"The direct survival-favoring versus survival-harming estimate is {fmt(analyses['positive_vs_negative'])}; it does not establish direction-specific enrichment. Without the expression restriction, the published catalog includes TOMM6 and yields {analyses['no_expression']['n']:,} genes, {analyses['no_expression']['A']} exposed targets and {analyses['no_expression']['overlap']} overlaps, with OR {fmt(analyses['no_expression'])}. All {len(analyses)} revised catalog intervals include1.",
      'precision_sentence':f"Only {p['overlap']} exposed outcome-positive genes and {info['n_informative_strata']} informative strata among {info['n_strata']} retained strata support the primary association. Holding observed margins fixed, the null expected overlap is {info['null_expected_overlap']:.2f}; the original lower-95%-bound criterion first succeeds at {info['minimum_overlap_for_criterion']} overlaps. Its conditional probability is {100*probs[1.5]:.1f}% if the common OR is1.5 and {100*probs[2.0]:.1f}% if it is2.0.",
      'continuous_sentence':f"We compare mean per-gene median LFC within each stratum and weight differences by the observed share of exposed genes, m_s/{p['A']}. This targets the exposure-group stratum composition rather than weighting strata equally. "+ ' '.join(f"For {r['replicate_label']}, the contrast is {r['exposure_standardized_mean_difference']:.4f} log2 units (95% gene-bootstrap percentile interval {r['gene_bootstrap_percentile_025']:.4f} to {r['gene_bootstrap_percentile_975']:.4f})." for r in continuous['summary']),
      'influence_sentence':f"Omitting each informative stratum in turn gives conditional OR estimates from {influence.conditional_mle.min():.4f} to {influence.conditional_mle.max():.4f}; all {len(influence)} intervals include1. This is an influence diagnostic with changing fixed margins, not a set of independent replication tests."
    }
    design=load('docs/source_design_comparison.json')
    tables={'source_design':table('Table 1. Source-design comparison [1,2].',design['headers'],design['rows'],[1.05,2.70,2.75],8.2),
      'contingency':table('Table 2. Primary published-catalog counts.', ['Survival grouping','Published lower content','No published lower-content call','Total'],
         [['Favoring',*p['table_2x2'][0],p['A']],['Other eligible',*p['table_2x2'][1],p['n']-p['A']]],[1.65,1.45,2.4,1.0],9),
      'candidates':table('Table 3. Primary candidates and recovered selected-guide evidence. S7 covers two of four guides per replicate.',
         ['Gene','Survival beta','Survival FDR','Median LFC r0 / r1','S7 positive FDR <0.10: r0 / r1 / same'],
         [[r.gene,f'{r.survival_beta:.3f}',f'{r.survival_fdr:.5g}',f'{r.median_lfc_r0:.3f} / {r.median_lfc_r1:.3f}',f"{int(recovered.loc[r.gene,'source_s7_positive_fdr_lt010_r0'])} / {int(recovered.loc[r.gene,'source_s7_positive_fdr_lt010_r1'])} / {int(recovered.loc[r.gene,'source_s7_positive_fdr_lt010_sameguide'])}"] for _,r in selected.iterrows()],
         [1.05,.8,.85,1.65,2.15],8.3)}
    tables['comparators']=table('Table 7. Higher-content comparison candidates from the same resource.',
       ['Gene','Survival beta','Survival FDR','Published content call'],
       [[r.gene,f'{r.survival_beta:.5f}',f'{r.survival_fdr:.6f}',r.published_content_label] for _,r in comparators.iterrows()],
       [1.15,1.35,1.35,2.65],8.8)
    annotations=load('docs/candidate_annotations.json')['rows'];annot={x['gene']:x for x in annotations}
    assert set(a['primary_lower_content_candidates'])|{'NUPL1'}==set(annot)
    bio_genes=['SSR2','NCOA6','CREB3','SLC8B1']
    tables['biology']=table('Table 4. Independent evidence that informs follow-up, with model-specific limits.',
       ['Gene / process','Direct evidence','Interpretation limit'],
       [[gene+' / '+({'SSR2':'ER proinsulin import','NCOA6':'Transcriptional coactivation','CREB3':'ER stress response','SLC8B1':'Mitochondrial sodium/calcium exchange'}[gene]),annot[gene]['direct_beta_cell_evidence']+f' [{6+i}]',annot[gene]['inference_limit']] for i,gene in enumerate(bio_genes)],
       [1.25,3.05,2.2],8.4)
    old={r['name']:r for r in load('results/legacy_protocol/results.json')['analyses']}
    tables['legacy_guides']=table('Table S1. Historical stored-zero diagnostics; not stronger-significance analyses.',
       ['Rule','N / exposed / overlaps','Conditional OR (95% CI)','One-sided enrichment p'],
       [[label,f"{old[name]['n']} / {old[name]['A']} / {old[name]['conflicts']}",f"{old[name]['conditional_mle']:.3f} ({old[name]['ci_lower']:.3f}-{old[name]['ci_upper']:.3f})",f"{old[name]['p_enrichment']:.3f}"] for name,label in [('three_directional_guides','Three per replicate; filtered universe'),('three_guides_fixed_universe','Three per replicate; fixed universe')]],
       [2.0,1.4,1.8,1.3],8.2)
    refs=[
      '1. Benaglio P, Zhu H, Okino ML, et al. Type 1 diabetes risk genes mediate pancreatic beta cell survival in response to proinflammatory cytokines. Cell Genomics. 2022;2(12):100214. https://doi.org/10.1016/j.xgen.2022.100214',
      '2. Rottner AK, Ye Y, Navarro-Guerrero E, et al. A genome-wide CRISPR screen identifies CALCOCO2 as a regulator of beta cell function influencing type 2 diabetes risk. Nature Genetics. 2023;55(1):54-65. Published online 21 December 2022. https://doi.org/10.1038/s41588-022-01261-2',
      '3. Okino ML, Zhu H, Corban S, et al. Functional genomics reveals mediators of beta cell survival in ER stress and type 2 diabetes risk. bioRxiv [preprint]. 2026;2026.03.30.715154. Posted 2 April 2026. https://doi.org/10.64898/2026.03.30.715154',
      '4. Howe LJ, Aulchenko YS, Davey Smith G, et al. Evaluating transportability of in vitro cellular models to in vivo human phenotypes using gene perturbation data. Nature Communications. 2026;17:513. First published online 13 December 2025. https://doi.org/10.1038/s41467-025-67199-1',
      '5. TGDStanford. TGDStanford/Rottner_etal_2022: Rottner_etal_2022. Version v1.0.0 [software]. Zenodo; 19 October 2022. https://doi.org/10.5281/zenodo.7226348',
      '6. Huang Y, Xu X, Arvan P, Liu M. Deficient endoplasmic reticulum translocon-associated protein complex limits the biosynthesis of proinsulin and insulin. The FASEB Journal. 2021;35(5):e21515. https://doi.org/10.1096/fj.202002774R',
      '7. Scoville DW, Cyphert HA, Liao L, et al. MLL3 and MLL4 methyltransferases bind to the MAFA and MAFB transcription factors to regulate islet beta-cell function. Diabetes. 2015;64(11):3772-3783. https://doi.org/10.2337/db15-0281',
      '8. Cnop M, Abdulkarim B, Bottu G, et al. RNA sequencing identifies dysregulation of the human pancreatic islet transcriptome by the saturated fatty acid palmitate. Diabetes. 2014;63(6):1978-1993. https://doi.org/10.2337/db13-1383',
      '9. Nita II, Hershfinkel M, Fishman D, et al. The mitochondrial Na+/Ca2+ exchanger upregulates glucose dependent Ca2+ signalling linked to insulin secretion. PLoS ONE. 2012;7(10):e46649. https://doi.org/10.1371/journal.pone.0046649',
      '10. Jiang H, Zheng S, Qian Y, et al. Restored UBE2C expression in islets promotes beta-cell regeneration in mice by ubiquitinating PER1. Cellular and Molecular Life Sciences. 2023;80(8):226. https://doi.org/10.1007/s00018-023-04868-8',
      "11. van Tienhoven R, O'Meally D, Scott TA, et al. Genetic protection from type 1 diabetes resulting from accelerated insulin mRNA decay. Cell. 2025;188(9):2407-2416.e9. https://doi.org/10.1016/j.cell.2025.02.018"]
    content=copy.deepcopy(narrative);content['authoritative_source']=False
    content['generated_from']='manuscript/narrative.json plus hash-bound scientific tables; edit narrative, not generated content'
    pages=[]
    def substitute(text):
        def rep(m):
            if m.group(1) not in facts:raise ValueError('Unbound manuscript fact '+m.group(1))
            return facts[m.group(1)]
        text=re.sub(r'\{\{(\w+)\}\}',rep,text)
        # Restore conventional readable spacing in inherited manuscript phrasing.
        for old,new in [('exact95%','exact 95%'),('exact95','exact 95'),('equal-tailed95%','equal-tailed 95%'),('all427','all 427'),('complete15','complete 15'),('the14','the 14'),('least20','least 20'),('least10','least 10'),('coverage3','coverage 3'),('versus6+','versus 6+'),('coverage3','coverage 3'),('versus4+','versus 4+'),('above1','above 1'),('include1','include 1'),('ratio1','ratio 1'),('OR1.5','OR 1.5'),('is1.5','is 1.5'),('is2.0','is 2.0'),('version1','version 1'),('Version2','Version 2')]:
            text=text.replace(old,new)
        return text
    for page in content['pages']:
        items=[]
        for item in page:
            if item['type']=='generated_table':items.append(tables[item['key']]);continue
            if item['type']=='generated_references':
                items.append({'type':'heading','text':'References'})
                items.extend({'type':'reference','text':text} for text in refs);continue
            for key in ['text','caption','alt']:
                if key in item:item[key]=substitute(item[key])
            items.append(item)
        pages.append(items)
    content['pages']=pages
    (ROOT/'manuscript/paper_content.json').write_text(json.dumps(content,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    report={'status':'PASS','method':'Generate every result table and headline fact from scientific outputs; fixed label and membership assertions block inconsistency.',
       'checks':['All427 candidates accounted with independent evidence and eligibility axes','Primary2x2 and percentages tied to model result','Candidate rows tied to candidate CSV','Higher-content comparators RAB14 and LRRC61 tied to candidate CSV','Information diagnostic counts tied to primary result','NUPL1 retains published evidence outside primary','Annotations cover all14 plus NUPL1','All placeholders resolved'],
       'generated_facts':facts,'dependencies':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in dependencies},
       'generated_content_sha256':hashlib.sha256((ROOT/'manuscript/paper_content.json').read_bytes()).hexdigest()}
    (ROOT/'evidence/manuscript_consistency.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return content

if __name__=='__main__':prepare()
