"""Independent downstream audit and SciPy inference for the screen comparison.

Does not import analyze.py, exact_stats.py, pandas, or the other verifier.
Source-parsing independence is supplied separately by independent_verify.py.
Here the public gene ledger is treated as an input, all analysis selections
and strata are rebuilt independently, and SciPy supplies the conditional PMFs.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import bisect
import csv
import hashlib
import json
import math
import platform
import sys

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.stats import hypergeom, nchypergeom_fisher

ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 1e-8


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def number(row, key):
    value = row.get(key, '')
    return float(value) if value else float('nan')


def truth(row, key):
    return row.get(key) == 'True'


def quantile(values, fraction):
    ordered = sorted(values)
    position = (len(ordered)-1)*fraction
    lo, hi = math.floor(position), math.ceil(position)
    return ordered[lo] + (ordered[hi]-ordered[lo])*(position-lo)


def audit_pass(row, threshold):
    label = row['published_content_label']
    if label not in ('Decrease', 'Increase'):
        return True
    sign = 'positive' if label == 'Decrease' else 'negative'
    return all(number(row, f'{sign}_support_{rep}') >= threshold for rep in ('r0', 'r1'))


def eligible(row, expression, directional=2):
    if not truth(row, 'survival_eligible') or not truth(row, 'insulin_covered'):
        return False
    if not audit_pass(row, directional):
        return False
    return expression is None or (
        not truth(row, 'expression_ambiguous')
        and math.isfinite(number(row, 'baseMean'))
        and number(row, 'baseMean') >= expression
    )


def independently_rebuild(ledger, name):
    expr, directional, fdr = 10, 2, .1
    if name == 'no_expression': expr = None
    elif name == 'expression_ge1': expr = 1
    elif name == 'expression_ge100': expr = 100
    elif name == 'survival_fdr005': fdr = .05
    elif name == 'three_directional_guides': directional = 3
    universe = [row for row in ledger if eligible(row, expr, directional)]
    bounds = sorted(set(quantile([math.log2(number(row,'baseMean')+1) for row in universe], q/5)
                        for q in range(6))) if expr is not None else []
    groups = [str(x) for x in range(1, len(bounds))] if bounds else ['0']
    labels = ['|'.join(items) for items in product(groups, ['3-5','6+'], ['3','4+'])]
    strata = {label: {'n':0,'m':0,'k':0,'a':0} for label in labels}
    genes = {}
    for row in universe:
        beta, survival_fdr = number(row,'survival_beta'), number(row,'survival_fdr')
        positive = beta > 0 and survival_fdr < fdr
        negative = beta < 0 and survival_fdr < .1
        if name == 'positive_vs_negative' and not (positive or negative):
            continue
        exposure = negative if name == 'survival_harming' else positive
        outcome_threshold = 3 if name == 'three_guides_fixed_universe' else directional
        outcome = row['published_content_label'] == 'Decrease' and audit_pass(row,outcome_threshold)
        expression_group = str(bisect.bisect_left(bounds[1:-1],math.log2(number(row,'baseMean')+1))+1) if bounds else '0'
        guide_group = '6+' if number(row,'survival_guides') >= 6 else '3-5'
        insulin_group = '4+' if number(row,'covered_insulin_guides') >= 4 else '3'
        label = '|'.join([expression_group,guide_group,insulin_group])
        genes[row['gene']] = (exposure,outcome,label)
        counts = strata[label]
        counts['n'] += 1
        counts['m'] += int(exposure)
        counts['k'] += int(outcome)
        counts['a'] += int(exposure and outcome)
    return genes, strata, bounds


def scipy_conditional(strata):
    observed = sum(row['a'] for row in strata.values())
    lower = sum(max(0,row['m']+row['k']-row['n']) for row in strata.values())
    upper = sum(min(row['m'],row['k']) for row in strata.values())
    terms = []
    for row in strata.values():
        n,m,k = row['n'],row['m'],row['k']
        lo,hi = max(0,m+k-n),min(m,k)
        if lo != hi:
            terms.append((n,m,k,np.arange(lo,hi+1)))
    def distribution(log_theta, central=False):
        result = np.array([1.0])
        for n,m,k,support in terms:
            p = (hypergeom.pmf(support,n,k,m) if central else
                 nchypergeom_fisher.pmf(support,n,k,m,math.exp(log_theta)))
            p /= p.sum()
            result = np.convolve(result,p)
        return result/result.sum()
    total_support = np.arange(lower,upper+1)
    def mean_difference(log_theta):
        return float(total_support @ distribution(log_theta))-observed
    def lower_limit_fn(log_theta):
        return float(distribution(log_theta)[observed-lower:].sum())-.025
    def upper_limit_fn(log_theta):
        return float(distribution(log_theta)[:observed-lower+1].sum())-.025
    if not terms:
        estimate, lo, hi = None, 0.0, 'Infinity'
    else:
        estimate = 0.0 if observed == lower else ('Infinity' if observed == upper else math.exp(brentq(mean_difference,-20,20,xtol=1e-12)))
        lo = 0.0 if observed == lower else math.exp(brentq(lower_limit_fn,-20,20,xtol=1e-12))
        hi = 'Infinity' if observed == upper else math.exp(brentq(upper_limit_fn,-20,20,xtol=1e-12))
    return {
        'conditional_mle': estimate,'ci_lower':lo,'ci_upper':hi,
        'p_enrichment':float(distribution(0,central=True)[observed-lower:].sum()),
        'support_lower':lower,'support_upper':upper,'observed_total':observed,
        'n_informative_strata':len(terms),
    }


def main():
    ledger = read_csv(ROOT/'results/gene_eligibility_ledger.csv')
    lookup = {row['gene']:row for row in ledger}
    results = json.loads((ROOT/'results/results.json').read_text())['analyses']
    definitions = json.loads((ROOT/'results/stratum_definitions.json').read_text())
    comparisons = []
    violations = []
    file_hashes = {}
    max_numeric_difference = 0.0
    for result in results:
        name = result['name']
        expected_genes,expected_strata,bounds = independently_rebuild(ledger,name)
        actual_genes = {row['gene']:(truth(row,'A'),truth(row,'L'),row['stratum'])
                        for row in read_csv(ROOT/f'results/{name}_genes.csv')}
        actual_strata = {row['stratum']:{key:int(row[key]) for key in ('n','m','k','a')}
                         for row in read_csv(ROOT/f'results/{name}_strata.csv')}
        gene_mismatches = sum(expected_genes.get(gene)!=actual_genes.get(gene)
                              for gene in set(expected_genes)|set(actual_genes))
        stratum_mismatches = sum(expected_strata.get(s)!=actual_strata.get(s)
                                 for s in set(expected_strata)|set(actual_strata))
        recorded_bounds = definitions['log2_baseMean_plus1_boundaries'][name]
        boundary_match = len(bounds)==len(recorded_bounds) and all(abs(a-b)<1e-12 for a,b in zip(bounds,recorded_bounds))
        independent = scipy_conditional(expected_strata)
        numeric_differences = {key:abs(independent[key]-result[key])
                               for key in ('conditional_mle','ci_lower','ci_upper','p_enrichment')}
        max_numeric_difference = max(max_numeric_difference,*numeric_differences.values())
        counts = {key:sum(s[key] for s in expected_strata.values()) for key in ('n','m','k','a')}
        n,m,k,a = (counts[key] for key in ('n','m','k','a'))
        table = [[a,m-a],[k-a,n-m-k+a]]
        table_match = table == result['table_2x2']
        fraction_match = math.isclose(a/m,result['fraction_A1'],abs_tol=1e-14) and math.isclose((k-a)/(n-m),result['fraction_A0'],abs_tol=1e-14)
        naive_odds = a*(n-m-k+a)/((m-a)*(k-a))
        fisher_p = float(hypergeom.sf(a-1,n,k,m))
        unstratified_match = math.isclose(naive_odds,result['unstratified']['sample_odds_ratio'],rel_tol=1e-12) and math.isclose(fisher_p,result['unstratified']['p_enrichment'],rel_tol=1e-10)
        pass_checks = gene_mismatches==stratum_mismatches==0 and boundary_match and table_match and fraction_match and unstratified_match and max(numeric_differences.values()) < TOLERANCE
        if not pass_checks: violations.append(name)
        comparisons.append({'name':name,'passed':pass_checks,'gene_count':len(expected_genes),
                            'gene_selection_label_or_stratum_mismatches':gene_mismatches,
                            'stratum_margin_mismatches':stratum_mismatches,'linear_quantile_boundaries_match':boundary_match,
                            'joint_table_matches':table_match,'fractions_match':fraction_match,'unstratified_sample_OR_and_Fisher_tail_match':unstratified_match,
                            'independent_scipy':independent,'absolute_numeric_differences':numeric_differences})
        for filename in [f'{name}_genes.csv',f'{name}_strata.csv']:
            file_hashes[f'results/{filename}']=hashlib.sha256((ROOT/'results'/filename).read_bytes()).hexdigest()

    candidates = read_csv(ROOT/'results/candidate_table_427.csv')
    expected_candidate_genes = {row['gene'] for row in ledger if truth(row,'survival_eligible') and number(row,'survival_beta')>0 and number(row,'survival_fdr')<.1}
    candidate_errors = []
    class_counts = Counter()
    classes = ['Supported cross-screen conflict','Survival-favoring with higher-content evidence','Content effect unresolved','Comparison unavailable or inconsistent']
    for row in candidates:
        source = lookup[row['gene']]
        klass = 'Comparison unavailable or inconsistent'
        if eligible(source,10):
            klass = 'Content effect unresolved'
            if source['published_content_label']=='Decrease': klass=classes[0]
            elif source['published_content_label']=='Increase': klass=classes[1]
        if row['candidate_class'] != klass:candidate_errors.append(row['gene'])
        class_counts[klass]+=1
    candidate_members_match = expected_candidate_genes == {row['gene'] for row in candidates} and len(candidates)==len(expected_candidate_genes)
    candidate_order = [(classes.index(row['candidate_class']),float(row['survival_fdr']),row['gene']) for row in candidates]
    candidate_sort_match = candidate_order == sorted(candidate_order)

    primary_genes,_,_ = independently_rebuild(ledger,'primary')
    diagnostic_groups = defaultdict(lambda: {False:[],True:[]})
    for gene,(exposure,outcome,stratum) in primary_genes.items():
        for rep in ('r0','r1'):
            diagnostic_groups[(stratum,rep)][exposure].append(number(lookup[gene],f'median_lfc_{rep}'))
    diagnostics = read_csv(ROOT/'results/threshold_free_stratum_diagnostics.csv')
    diagnostic_errors=[]
    signs=defaultdict(dict)
    for row in diagnostics:
        stratum,rep=row['stratum'],row['replicate']
        values=diagnostic_groups[(stratum,rep)]
        for exposure in (0,1):
            arr=values[bool(exposure)]
            if len(arr)!=int(row[f'n_A{exposure}']):diagnostic_errors.append(f'{stratum}/{rep}/n{exposure}')
            for label,q in [('q25',.25),('median',.5),('q75',.75)]:
                got=row[f'{label}_A{exposure}']
                if (not arr and got!='') or (arr and not math.isclose(quantile(arr,q),float(got),abs_tol=1e-12)):
                    diagnostic_errors.append(f'{stratum}/{rep}/{label}{exposure}')
        if values[False] and values[True]:
            difference=quantile(values[True],.5)-quantile(values[False],.5)
            if not math.isclose(difference,float(row['median_difference_A1_minus_A0']),abs_tol=1e-12):diagnostic_errors.append(f'{stratum}/{rep}/difference')
            signs[stratum][rep]=1 if difference>0 else (-1 if difference<0 else 0)
    direction_summary={rep:dict(Counter(s[rep] for s in signs.values())) for rep in ('r0','r1')}
    opposite_sign_strata=[name for name,s in signs.items() if s['r0']*s['r1']<0]
    inputs=json.loads((ROOT/'evidence/analysis_inputs.json').read_text())
    analysis_input_hashes_match=all(hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256'] for item in inputs['files'])
    for path in ['results/results.json','results/gene_eligibility_ledger.csv','results/candidate_table_427.csv','results/threshold_free_stratum_diagnostics.csv','code/analyze.py','code/exact_stats.py','config.json','evidence/analysis_inputs.json','tests/run_statistical_crosscheck.py']:
        file_hashes[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    passed=not violations and not candidate_errors and candidate_members_match and candidate_sort_match and not diagnostic_errors and analysis_input_hashes_match
    evidence={
        'status':'passed' if passed else 'failed','reviewed_at_utc':datetime.now(timezone.utc).isoformat(),
        'scope':'Independent downstream analysis selection and stratification from the public gene ledger, numerical inference, candidate classification and descriptive diagnostics. Raw source reconstruction is a separate evidence/independent_verification.json task.',
        'independence':'Does not import analyze.py, exact_stats.py, pandas or independent_verify.py. SciPy Fisher noncentral hypergeometric PMFs and Brent roots; central hypergeometric PMFs for null tails; explicit linear quantiles and bisect for strata.',
        'python_version':platform.python_version(),'numpy_version':np.__version__,'scipy_version':scipy.__version__,
        'numeric_tolerance':TOLERANCE,'maximum_absolute_numeric_difference':max_numeric_difference,
        'analyses':comparisons,'failed_analyses':violations,
        'candidate_check':{'membership_matches':candidate_members_match,'classification_mismatches':candidate_errors,'sorting_matches':candidate_sort_match,'class_counts':dict(class_counts)},
        'threshold_free_check':{'quantile_or_count_mismatches':diagnostic_errors,'strata_with_both_exposure_groups':len(signs),'median_contrast_sign_counts':direction_summary,'opposite_sign_strata':opposite_sign_strata},
        'analysis_input_hashes_match':analysis_input_hashes_match,
        'file_sha256':file_hashes,'command':'python tests/run_statistical_crosscheck.py',
    }
    (ROOT/'evidence/statistical_crosscheck.json').write_text(json.dumps(evidence,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({key:evidence[key] for key in ['status','maximum_absolute_numeric_difference','failed_analyses','candidate_check','threshold_free_check','analysis_input_hashes_match']},indent=2))
    return 0 if passed else 1


if __name__=='__main__':
    sys.exit(main())
