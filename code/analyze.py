"""Reproduce the beta-cell screen comparison from hash-pinned source tables."""
from __future__ import annotations
import argparse, hashlib, itertools, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from exact_stats import exact_conditional

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'data/raw'; EVID=ROOT/'evidence'; OUT=ROOT/'results'
CONFIG=json.loads((ROOT/'config.json').read_text())
INS_BOOK=RAW/'41588_2022_1261_MOESM4_ESM.xlsx'

def write_json(path, obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False,default=lambda x: x.item() if isinstance(x,np.generic) else str(x))+'\n',encoding='utf-8')

def verify_hashes():
    manifest=json.loads((EVID/'input_manifest.json').read_text())
    for item in manifest['inputs']:
        path=ROOT/item['path']
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=item['sha256']: raise RuntimeError(f'Input hash changed: {path.name}')

def symbol(series):
    return series.astype('string').str.strip().replace('',pd.NA)

def read_sources():
    verify_hashes()
    s=pd.read_csv(RAW/'Supplementary_Data_09.CRISPR-KO-screen_results.csv',keep_default_na=False)
    s['survival_source_row']=np.arange(2,len(s)+2)
    s=s.rename(columns={'Gene':'survival_source_symbol','sgRNA number':'survival_guides','Cytokine survival beta':'survival_beta','Cytokine survival fdr':'survival_fdr'})
    s['gene']=symbol(s.survival_source_symbol)
    for c in ['survival_guides','survival_beta','survival_fdr']: s[c]=pd.to_numeric(s[c],errors='coerce')
    s['survival_unique']=~s.gene.duplicated(keep=False)&s.gene.notna()
    s['survival_finite']=np.isfinite(s[['survival_beta','survival_fdr','survival_guides']]).all(axis=1)
    s['survival_eligible']=s.survival_unique&s.survival_finite&(s.survival_guides>=3)&~s.gene.str.startswith('hsa-',na=False)
    if not s.survival_unique.all(): raise RuntimeError('Nonunique survival identity requires adjudication')

    book=pd.ExcelFile(INS_BOOK,engine='openpyxl')
    h=pd.read_excel(book,sheet_name='Table S1-Screen hits ',usecols='A:C')
    h['hit_source_row']=np.arange(2,len(h)+2)
    h=h.rename(columns={'Gene':'insulin_hit_source_symbol','Mean log fold change':'published_mean_lfc','Effect on insulin content':'published_content_label'})
    h['gene']=symbol(h.insulin_hit_source_symbol)
    if h.gene.isna().any() or h.gene.duplicated().any(): raise RuntimeError('Nonunique insulin hit identity')
    if not h.published_content_label.isin(['Decrease','Increase']).all(): raise RuntimeError('Unknown content label')

    g=pd.read_excel(book,sheet_name='Table S8-Screen_all sgRNA')
    g['guide_source_row']=np.arange(2,len(g)+2)
    g['insulin_guide_source_symbol']=g.Gene
    g['gene']=symbol(g.Gene);g['source_guide_id']=g.sgrna
    clean=symbol(g.sgrna)
    g['replicate']=clean.str.extract(r'_(r[01])$',expand=False)
    g['guide_key']=clean.str.replace(r'_r[01]$','',regex=True)
    numeric=['LFC','FDR','control_count','treatment_count']
    for c in numeric:g[c]=pd.to_numeric(g[c],errors='coerce')
    g['finite']=np.isfinite(g[numeric]).all(axis=1)
    g['unique_key_replicate']=~g.duplicated(['guide_key','replicate'],keep=False)
    g['valid']=g.finite&g.unique_key_replicate&g.gene.notna()&g.replicate.notna()&g.guide_key.notna()
    g['count_sum']=g.control_count+g.treatment_count
    v=g[g.valid].copy()
    pairs=v.groupby('guide_key').agg(rows=('gene','size'),n_genes=('gene','nunique'),n_reps=('replicate','nunique'),min_sum=('count_sum','min'),gene=('gene','first'))
    covered_keys=pairs.index[(pairs.rows==2)&(pairs.n_genes==1)&(pairs.n_reps==2)&(pairs.min_sum>=20)]
    g['covered']=g.guide_key.isin(covered_keys)&g.valid
    covered=pairs.loc[covered_keys].groupby('gene').size().rename('covered_insulin_guides')
    v['pos_support']=(v.LFC>0)&(v.FDR<.1)
    v['neg_support']=(v.LFC<0)&(v.FDR<.1)
    agg=v.groupby(['gene','replicate']).agg(median_lfc=('LFC','median'),positive_support=('pos_support','sum'),negative_support=('neg_support','sum'),valid_guides=('guide_key','nunique'),source_rows=('guide_source_row',lambda x:';'.join(map(str,x))))
    a=agg.unstack('replicate');a.columns=[f'{x}_{y}' for x,y in a.columns]
    a=a.join(covered).join(g.groupby('gene').insulin_guide_source_symbol.first())
    a['covered_insulin_guides']=a.covered_insulin_guides.fillna(0).astype(int)
    a['insulin_covered']=(a.covered_insulin_guides>=3)&~a.index.isin(['EGFP','LacZ','luciferase'])

    e=pd.read_excel(book,sheet_name='Table S5 - RNASeq',usecols='A:C')
    e['expression_source_row']=np.arange(2,len(e)+2)
    e['expression_source_symbol']=e.GeneName;e['gene']=symbol(e.GeneName)
    e['baseMean']=pd.to_numeric(e.baseMean,errors='coerce')
    e['expression_ambiguous']=e.gene.duplicated(keep=False)&e.gene.notna()
    # All duplicate-name rows are excluded, not arbitrarily deduplicated.
    expression=e[e.gene.notna()&~e.expression_ambiguous].set_index('gene')[['baseMean','expression_source_row','expression_source_symbol','GeneID']]
    ambiguous=set(e.loc[e.expression_ambiguous,'gene'])
    f=s.set_index('gene').join(a,how='outer').join(h.set_index('gene'),how='outer').join(expression,how='left')
    f.index.name='gene'
    f['expression_ambiguous']=f.index.isin(ambiguous)
    f['published_content_label']=f.published_content_label.fillna('Not a published hit')
    for c in ['survival_eligible','insulin_covered']:f[c]=f[c].fillna(False).astype(bool)
    for c in ['positive_support_r0','positive_support_r1','negative_support_r0','negative_support_r1','covered_insulin_guides']:
        f[c]=f[c].fillna(0).astype(int)
    f['direction_audit_pass']=audit_mask(f,2)
    f['stronger_direction_audit_pass']=audit_mask(f,3)
    f['in_U0']=f.survival_eligible&f.insulin_covered
    f['in_U']=f.in_U0&~f.expression_ambiguous&np.isfinite(f.baseMean)&(f.baseMean>=10)&f.direction_audit_pass
    reasons=[]
    for gene,row in f.iterrows():
        r=[]
        if pd.isna(row.survival_source_symbol):r.append('no_exact_survival_symbol')
        elif not row.survival_eligible:
            if str(gene).startswith('hsa-'):r.append('hsa_prefix')
            if not row.survival_finite:r.append('nonfinite_survival')
            if row.survival_guides<3:r.append('survival_guides_lt3')
        if pd.isna(row.insulin_guide_source_symbol):r.append('no_exact_insulin_symbol')
        elif not row.insulin_covered:r.append('insulin_control' if gene in ['EGFP','LacZ','luciferase'] else 'covered_insulin_guides_lt3')
        if row.expression_ambiguous:r.append('ambiguous_expression_symbol')
        elif pd.isna(row.expression_source_symbol):r.append('no_exact_expression_symbol')
        elif not np.isfinite(row.baseMean):r.append('nonfinite_expression')
        elif row.baseMean<10:r.append('expression_lt10')
        if not row.direction_audit_pass:r.append('failed_published_direction_audit')
        reasons.append(';'.join(r) if r else 'eligible')
    f['eligibility_reason']=reasons
    audit={
      'survival_rows':len(s),'survival_eligible':int(s.survival_eligible.sum()),
      'survival_favoring':int((s.survival_eligible&(s.survival_beta>0)&(s.survival_fdr<.1)).sum()),
      'survival_harming':int((s.survival_eligible&(s.survival_beta<0)&(s.survival_fdr<.1)).sum()),
      'insulin_hit_rows':len(h),'insulin_decrease_hits':int((h.published_content_label=='Decrease').sum()),'insulin_increase_hits':int((h.published_content_label=='Increase').sum()),
      'insulin_guide_rows':len(g),'insulin_guide_replicates':{str(k):int(x) for k,x in g.replicate.value_counts().items()},
      'insulin_gene_labels':int(g.gene.nunique()),'covered_insulin_targets':int(a.insulin_covered.sum()),
      'U0':int(f.in_U0.sum()),'U':int(f.in_U.sum()),
      'primary_survival_favoring_margin':int((f.in_U&(f.survival_beta>0)&(f.survival_fdr<.1)).sum()),
      'primary_decrease_margin':int((f.in_U&(f.published_content_label=='Decrease')).sum()),
      'duplicated_expression_names':len(ambiguous),'invalid_guide_rows':int((~g.valid).sum()),
      'failed_direction_audits':f.index[~f.direction_audit_pass].tolist(),
      'SOCS1_survival_beta':float(f.loc['SOCS1','survival_beta']),
      'INS_replica_medians':[float(f.loc['INS',f'median_lfc_{rep}']) for rep in ['r0','r1']],
      'TOMM6_baseMean':float(f.loc['TOMM6','baseMean']),
      'TOMM6_positive_support':[int(f.loc['TOMM6',f'positive_support_{rep}']) for rep in ['r0','r1']],
      'survival_denominator_discrepancy':{'article':18703,'source_filter':18714,'difference':11,'disposition':'retain explicit deposited membership; author notebook contingency margins also total 18714'}
    }
    expected={'survival_rows':20708,'survival_eligible':18714,'survival_favoring':427,'survival_harming':440,'insulin_hit_rows':580,'insulin_decrease_hits':518,'insulin_increase_hits':62,'insulin_guide_rows':142180,'insulin_gene_labels':18056,'covered_insulin_targets':17597,'U0':17276,'U':12184,'primary_survival_favoring_margin':285,'primary_decrease_margin':479,'duplicated_expression_names':75}
    mismatch={k:{'expected':v,'actual':audit[k]} for k,v in expected.items() if audit[k]!=v}
    audit['marginal_check_mismatches']=mismatch
    write_json(EVID/'marginal_audit.json',audit)
    if mismatch:raise RuntimeError(f'Marginal source audit failed: {mismatch}')
    if audit['failed_direction_audits']!=['TOMM6'] or audit['SOCS1_survival_beta']>=0 or min(audit['INS_replica_medians'])<=0:raise RuntimeError('Sign audit failed')
    return f,g,e,audit

def audit_mask(f,n):
    dec=(f.positive_support_r0>=n)&(f.positive_support_r1>=n)
    inc=(f.negative_support_r0>=n)&(f.negative_support_r1>=n)
    return (~f.published_content_label.isin(['Decrease','Increase']))|((f.published_content_label=='Decrease')&dec)|((f.published_content_label=='Increase')&inc)

def stratum_labels(df,expression=True,bins=None):
    df=df.copy()
    if expression:
        z=np.log2(df.baseMean+1)
        if bins is None:
            q,bins=pd.qcut(z,5,labels=False,retbins=True,duplicates='drop')
        else:q=pd.cut(z,bins,labels=False,include_lowest=True,right=True)
        if q.isna().any():raise RuntimeError('Stratum assignment missing')
        df['expression_stratum']=q.astype(int)+1
        levels=list(range(1,len(bins)))
    else:
        bins=[];df['expression_stratum']=0;levels=[0]
    df['survival_guide_stratum']=np.where(df.survival_guides>=6,'6+','3-5')
    df['insulin_guide_stratum']=np.where(df.covered_insulin_guides>=4,'4+','3')
    df['stratum']=df.expression_stratum.astype(str)+'|'+df.survival_guide_stratum+'|'+df.insulin_guide_stratum
    combinations=[f'{q}|{s}|{i}' for q,s,i in itertools.product(levels,['3-5','6+'],['3','4+'])]
    return df, list(map(float,bins)),combinations

def summarize(name,df,combinations):
    strata=[]
    for label in combinations:
        z=df[df.stratum==label];n=len(z);m=int(z.A.sum());k=int(z.L.sum());a=int((z.A&z.L).sum())
        strata.append({'stratum':label,'n':n,'m':m,'k':k,'a':a,'A1_L0':m-a,'A0_L1':k-a,'A0_L0':n-m-k+a,'informative':max(0,m+k-n)<min(m,k)})
    rows=pd.DataFrame(strata);rows.to_csv(OUT/f'{name}_strata.csv',index=False)
    exact=exact_conditional(strata)
    n=len(df);m=int(df.A.sum());k=int(df.L.sum());a=int((df.A&df.L).sum())
    table=np.array([[a,m-a],[k-a,n-m-k+a]],dtype=int)
    fi=fisher_exact(table,alternative='greater');co=odds_ratio(table,kind='conditional');ci=co.confidence_interval(.95)
    def finite(x):return float(x) if np.isfinite(x) else ('Infinity' if x>0 else None)
    exact.update({'name':name,'n':n,'A':m,'L':k,'conflicts':a,'table_2x2':table.tolist(),
                  'fraction_A1':a/m if m else None,'fraction_A0':(k-a)/(n-m) if n>m else None,
                  'unstratified':{'sample_odds_ratio':finite(fi.statistic),'conditional_mle':finite(co.statistic),'ci_lower':finite(ci.low),'ci_upper':finite(ci.high),'p_enrichment':float(fi.pvalue)}})
    exact['conclusion']='supported' if exact['ci_lower']>1 else ('contradicted' if float(exact['ci_upper'])<=1 else 'unresolved')
    if exact['status']=='non_estimable':exact['conclusion']='non_estimable'
    df.to_csv(OUT/f'{name}_genes.csv',index=True)
    return exact

def run(f,g,e,audit):
    OUT.mkdir(exist_ok=True)
    # Verify the public analysis specification and scientific input files.
    inputs=json.loads((EVID/'analysis_inputs.json').read_text())
    for item in inputs['files']:
        if hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()!=item['sha256']:
            raise RuntimeError(f'Analysis input changed: {item["path"]}')
    f['A']=(f.survival_beta>0)&(f.survival_fdr<.1)&f.survival_eligible
    f['L']=(f.published_content_label=='Decrease')&f.direction_audit_pass
    f['survival_harming']=(f.survival_beta<0)&(f.survival_fdr<.1)&f.survival_eligible
    f.to_csv(OUT/'gene_eligibility_ledger.csv',index=True)
    g[['guide_source_row','source_guide_id','insulin_guide_source_symbol','gene','guide_key','replicate','LFC','FDR','control_count','treatment_count','valid','covered']].to_csv(OUT/'guide_evidence.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    e.to_csv(OUT/'expression_identifier_ledger.csv',index=False)
    frames={};bounds={};combos={}
    specs=[('primary',10,2,.1),('no_expression',None,2,.1),('expression_ge1',1,2,.1),('expression_ge100',100,2,.1),('survival_fdr005',10,2,.05),('three_directional_guides',10,3,.1)]
    for name,threshold,directional,fdr in specs:
        mask=f.in_U0&audit_mask(f,directional)
        if threshold is not None:mask&=~f.expression_ambiguous&np.isfinite(f.baseMean)&(f.baseMean>=threshold)
        z=f[mask].copy();z['A']=(z.survival_beta>0)&(z.survival_fdr<fdr)
        z['L']=(z.published_content_label=='Decrease')&(z.positive_support_r0>=directional)&(z.positive_support_r1>=directional)
        frames[name],bounds[name],combos[name]=stratum_labels(z,threshold is not None)
    for name in ['survival_harming','positive_vs_negative','three_guides_fixed_universe']:
        z=frames['primary'].copy()
        if name=='survival_harming':z['A']=z.survival_harming
        elif name=='positive_vs_negative':z=z[z.A|z.survival_harming].copy()
        else:z['L']=(z.published_content_label=='Decrease')&(z.positive_support_r0>=3)&(z.positive_support_r1>=3)
        frames[name]=z;bounds[name]=bounds['primary'];combos[name]=combos['primary']
    # Save boundaries before constructing overlaps, as specified in the protocol.
    write_json(OUT/'stratum_definitions.json',{'assignment':'right-closed, lowest included','log2_baseMean_plus1_boundaries':bounds,'all_strata':combos})
    results=[summarize(name,df,combos[name]) for name,df in frames.items()]
    write_json(OUT/'results.json',{'analyses':results,'source_margins':audit})
    pd.DataFrame([{k:v for k,v in x.items() if k not in ['unstratified','table_2x2']} for x in results]).to_csv(OUT/'analysis_summary.csv',index=False)

    candidates=f[f.A].copy()
    candidates['candidate_class']='Comparison unavailable or inconsistent'
    candidates.loc[candidates.in_U,'candidate_class']='Content effect unresolved'
    candidates.loc[candidates.in_U&(candidates.published_content_label=='Increase'),'candidate_class']='Survival-favoring with higher-content evidence'
    candidates.loc[candidates.in_U&candidates.L,'candidate_class']='Supported cross-screen conflict'
    candidates['stronger_evidence_flag']=candidates.published_content_label.isin(['Decrease','Increase'])&candidates.stronger_direction_audit_pass
    candidates['external_functional_evidence_status']='Unknown; not systematically assessed; this table supplies screen evidence only'
    order=['Supported cross-screen conflict','Survival-favoring with higher-content evidence','Content effect unresolved','Comparison unavailable or inconsistent']
    candidates['class_order']=candidates.candidate_class.map(dict(zip(order,range(4))))
    candidates=candidates.reset_index().sort_values(['class_order','survival_fdr','gene']).drop(columns='class_order')
    assert len(candidates)==427
    candidates.to_csv(OUT/'candidate_table_427.csv',index=False)
    write_json(OUT/'candidate_accounting.json',{'total':427,'primary_eligible':int(candidates.in_U.sum()),'classes':candidates.candidate_class.value_counts().to_dict(),'unavailable_reasons':candidates.loc[~candidates.in_U,'eligibility_reason'].value_counts().to_dict()})
    diagnostic=[]
    p=frames['primary']
    for label in combos['primary']:
        sub=p[p.stratum==label]
        for rep in ['r0','r1']:
            vals=[sub.loc[sub.A==v,f'median_lfc_{rep}'].to_numpy() for v in [False,True]]
            row={'stratum':label,'replicate':rep,'n_A0':len(vals[0]),'n_A1':len(vals[1])}
            for v,arr in enumerate(vals):
                for quant,q in [('q25',.25),('median',.5),('q75',.75)]:row[f'{quant}_A{v}']=float(np.quantile(arr,q)) if len(arr) else None
            row['median_difference_A1_minus_A0']=(row['median_A1']-row['median_A0']) if len(vals[0]) and len(vals[1]) else None
            diagnostic.append(row)
    pd.DataFrame(diagnostic).to_csv(OUT/'threshold_free_stratum_diagnostics.csv',index=False)
    plot(frames,results,candidates)
    print(json.dumps({'primary':results[0],'classes':candidates.candidate_class.value_counts().to_dict()},indent=2))

def plot(frames,results,candidates):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    folder=ROOT/'manuscript/figures';folder.mkdir(parents=True,exist_ok=True)
    p=frames['primary'];colors={'Content effect unresolved':'#365f9b','Survival-favoring with higher-content evidence':'#287b69','Supported cross-screen conflict':'#b33c35'}
    fig,axes=plt.subplots(1,2,figsize=(10,4),sharex=True,sharey=True)
    for ax,rep in zip(axes,['r0','r1']):
        ax.scatter(p.survival_beta,p[f'median_lfc_{rep}'],s=4,c='#b7bec5',alpha=.35,rasterized=True)
        for label,color in colors.items():
            z=candidates[candidates.candidate_class==label]
            ax.scatter(z.survival_beta,z[f'median_lfc_{rep}'],s=24 if label=='Supported cross-screen conflict' else 18,color=color,label=label,alpha=.85)
        ax.axhline(0,color='#777777',lw=.7);ax.axvline(0,color='#777777',lw=.7)
        ax.set_title(f'Insulin screen {rep}');ax.set_xlabel('Cytokine survival beta (MAGeCK)')
    axes[0].set_ylabel('Median insulin low/high guide LFC (log2)\nPositive values indicate lower stored insulin')
    fig.legend(*axes[1].get_legend_handles_labels(),loc='lower center',ncol=1,fontsize=8,bbox_to_anchor=(.5,-.05))
    fig.tight_layout(rect=[0,.16,1,1]);fig.savefig(folder/'signed_effects.png',bbox_inches='tight');fig.savefig(folder/'signed_effects.pdf',bbox_inches='tight');plt.close(fig)
    names=[x['name'] for x in results];labels=['Primary','No expression filter','Expression >=1','Expression >=100','Survival FDR <0.05','Three guides, exclusions','Survival-harming','Positive vs negative','Three guides, fixed U']
    fig,ax=plt.subplots(figsize=(7.8,4.5));y=np.arange(len(results))
    for i,x in enumerate(results):
        est=x['conditional_mle'];lo=x['ci_lower'];hi=x['ci_upper']
        ax.errorbar(est,i,xerr=[[est-lo],[hi-est]],fmt='o',color='#365f9b' if i==0 else '#555555',capsize=3)
    ax.axvline(1,ls='--',color='#aaaaaa');ax.set_xscale('log');ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlabel('Conditional odds ratio with exact 95% interval')
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    ax.xaxis.set_major_locator(FixedLocator([.3,.5,1,2,4]));ax.xaxis.set_major_formatter(FixedFormatter(['0.3','0.5','1','2','4']));ax.xaxis.set_minor_locator(NullLocator())
    fig.tight_layout();fig.savefig(folder/'sensitivity_forest.png');fig.savefig(folder/'sensitivity_forest.pdf');plt.close(fig)
    for rep in ['r0','r1']:
        fig,axes=plt.subplots(5,4,figsize=(8.5,10),sharex=True,sharey=True)
        xlo=float(p[f'median_lfc_{rep}'].min())-.15
        xhi=float(p[f'median_lfc_{rep}'].max())+.15
        for ax,label in zip(axes.flat, sorted(p.stratum.unique(),key=lambda s:(int(s.split('|')[0]),s.split('|')[1],s.split('|')[2]))):
            z=p[p.stratum==label]
            for exposure,color in [(False,'#777777'),(True,'#b33c35')]:
                arr=np.sort(z.loc[z.A==exposure,f'median_lfc_{rep}'])
                if len(arr):
                    ax.step(np.r_[xlo,arr,xhi],np.r_[0,np.arange(1,len(arr)+1)/len(arr),1],where='post',color=color,label=f'A={int(exposure)} n={len(arr)}')
            ax.set_title(label,fontsize=9);ax.legend(fontsize=8,loc='lower right');ax.axvline(0,color='#cccccc',lw=.5)
            from matplotlib.ticker import MaxNLocator
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4));ax.set_xlim(xlo,xhi)
        for ax in axes.flat[len(p.stratum.unique()):]:ax.set_visible(False)
        fig.supxlabel(f'Median insulin guide LFC ({rep}; log2 low/high)');fig.supylabel('Empirical cumulative fraction');fig.tight_layout()
        fig.savefig(folder/f'threshold_free_{rep}.png');fig.savefig(folder/f'threshold_free_{rep}.pdf');plt.close(fig)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--audit-only',action='store_true');args=ap.parse_args()
    f,g,e,audit=read_sources()
    if args.audit_only:print(json.dumps(audit,indent=2))
    else:run(f,g,e,audit)
