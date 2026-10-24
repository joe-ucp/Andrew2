"""Post-review published-catalog analysis; retain V1 arithmetic as historical evidence.

S8 stored zero is NOT interpreted as a validated q-value threshold here.
"""
import json, shutil
from pathlib import Path
import numpy as np
import pandas as pd
from analyze import stratum_labels, summarize, write_json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'

def run():
    current=json.loads((OUT/'results.json').read_text())
    if current.get('version')==2:
        raise RuntimeError('Already revised: run reproduce.py to regenerate V1 inputs before revising; historical outputs will not be overwritten.')
    historical=OUT/'legacy_protocol'; historical.mkdir(exist_ok=True)
    names=['primary','no_expression','expression_ge1','expression_ge100','survival_fdr005',
           'three_directional_guides','survival_harming','positive_vs_negative','three_guides_fixed_universe']
    for name in ['results.json','analysis_summary.csv','candidate_table_427.csv','candidate_accounting.json',
                 'gene_eligibility_ledger.csv','stratum_definitions.json']+[n+s for n in names for s in ['_genes.csv','_strata.csv']]:
        shutil.copy2(OUT/name,historical/name)
    old=json.loads((historical/'results.json').read_text())
    f=pd.read_csv(historical/'gene_eligibility_ledger.csv',keep_default_na=False,na_values=['']).set_index('gene')
    f['L']=f.published_content_label.eq('Decrease')
    f['in_U']=f.in_U0&~f.expression_ambiguous&np.isfinite(f.baseMean)&(f.baseMean>=10)
    f['eligibility_reason']=f.eligibility_reason.str.replace(';failed_published_direction_audit','',regex=False).str.replace('failed_published_direction_audit','',regex=False).replace('','eligible')
    f['primary_eligible']=f.in_U
    f['primary_eligibility_status']=np.where(f.in_U,'Eligible','Excluded')
    f['screen_evidence_status']=np.select([
        f.published_content_label.eq('Decrease'), f.published_content_label.eq('Increase'),
        f.insulin_guide_source_symbol.isna()],
        ['Published lower-content evidence','Published higher-content evidence','Screen measurement unavailable'],
        default='No published directional call')
    f['guide_significance_status']='S8 significance invalidated by Int64 coercion; selected S7 support recorded separately'
    f=f.rename(columns={**{f'{sign}_support_{rep}':f'stored_zero_{sign}_count_{rep}' for sign in ['positive','negative'] for rep in ['r0','r1']},
                        'direction_audit_pass':'legacy_stored_zero_rule_2_each',
                        'stronger_direction_audit_pass':'legacy_stored_zero_rule_3_each'})
    g=pd.read_csv(OUT/'guide_evidence.csv.gz',keep_default_na=False,na_values=[''])
    v=g[g.valid].copy()
    pairs=v.pivot(index=['gene','guide_key'],columns='replicate',values=['LFC','FDR'])
    pairrows=[]
    for (gene,guide),row in pairs.iterrows():
        a,b=row[('LFC','r0')],row[('LFC','r1')]
        zero=(row[('FDR','r0')]==0)&(row[('FDR','r1')]==0)
        pairrows.append({'gene':gene,'guide_key':guide,'lfc_r0':a,'lfc_r1':b,
          'stored_fdr_r0':row[('FDR','r0')],'stored_fdr_r1':row[('FDR','r1')],
          'same_positive_sign':a>0 and b>0,'opposite_signs':a*b<0,
          'same_positive_stored_zero':a>0 and b>0 and zero,
          'same_negative_stored_zero':a<0 and b<0 and zero})
    gp=pd.DataFrame(pairrows)
    gp=gp.sort_values(['gene','guide_key'])
    gp['guide_index_within_gene']=gp.groupby('gene').cumcount()+1
    counts=gp.groupby('gene')[['same_positive_sign','opposite_signs','same_positive_stored_zero','same_negative_stored_zero']].sum().add_suffix('_guides')
    f=f.join(counts)
    frames={};bounds={};combos={}
    for name,threshold,fdr in [('primary',10,.1),('no_expression',None,.1),('expression_ge1',1,.1),('expression_ge100',100,.1),('survival_fdr005',10,.05)]:
        mask=f.in_U0.copy()
        if threshold is not None:mask &= ~f.expression_ambiguous&np.isfinite(f.baseMean)&(f.baseMean>=threshold)
        z=f[mask].copy();z['A']=(z.survival_beta>0)&(z.survival_fdr<fdr)
        frames[name],bounds[name],combos[name]=stratum_labels(z,threshold is not None)
    for name in ['survival_harming','positive_vs_negative']:
        z=frames['primary'].copy()
        if name=='survival_harming':z['A']=z.survival_harming
        else:z=z[z.A|z.survival_harming].copy()
        frames[name]=z;bounds[name]=bounds['primary'];combos[name]=combos['primary']
    results=[summarize(name,df,combos[name]) for name,df in frames.items()]
    for r in results:
        r['overlap']=r.pop('conflicts')
        r['interpretation']='Published Table S1 catalog comparison; guide FDR significance not reconstructed'
        r['design_status']='Post-review source-semantic amendment; exact-symbol primary membership checked against original lock'
    assert set(frames['primary'].index)==set(pd.read_csv(historical/'primary_genes.csv').gene)
    original=old['analyses'][0]
    assert all(results[0][k]==original[k] for k in ['n','A','L','conditional_mle','ci_lower','ci_upper','p_enrichment'])
    assert results[0]['overlap']==14
    write_json(OUT/'results.json',{'version':2,'analyses':results,'semantic_status':'GUIDE_SIGNIFICANCE_UNRESOLVED',
       'primary_membership_and_estimates_unchanged':True,'legacy_results':'legacy_protocol/results.json',
       'source_margins':old['source_margins'], 'source_margins_note':'Historical margins include original stored-value audit; revised catalogs are described by analyses above.'})
    pd.DataFrame([{k:v for k,v in r.items() if k not in ['unstratified','table_2x2']} for r in results]).to_csv(OUT/'analysis_summary.csv',index=False)
    f.to_csv(OUT/'gene_eligibility_ledger.csv')
    write_json(OUT/'stratum_definitions.json',{'assignment':'right-closed, lowest included','log2_baseMean_plus1_boundaries':bounds,'all_strata':combos})
    c=f[f.A].copy().reset_index()
    order={'Published lower-content evidence':0,'Published higher-content evidence':1,'No published directional call':2,'Screen measurement unavailable':3}
    c['sort']=c.screen_evidence_status.map(order)
    c=c.sort_values(['sort','survival_fdr','gene']).drop(columns='sort')
    recovered=pd.read_csv(ROOT/'evidence/candidate_S7_selected_guide_evidence.csv',keep_default_na=False,na_values=[''])
    c=c.merge(recovered,on='gene',how='left',validate='one_to_one')
    c['external_functional_evidence_status']='See bounded candidate annotations; no matched survival/content/secretion validation performed'
    c.to_csv(OUT/'candidate_table_427.csv',index=False)
    gp[gp.gene.isin(c.gene)].to_csv(OUT/'candidate_paired_guide_effects.csv',index=False)
    write_json(OUT/'candidate_accounting.json',{'total':len(c),'common_coverage':int(c.in_U0.sum()),'primary_eligible':int(c.in_U.sum()),
       'coverage_losses':int((~c.in_U0).sum()),'additional_primary_losses':int((c.in_U0&~c.in_U).sum()),
       'screen_evidence_status_counts':c.screen_evidence_status.value_counts().to_dict(),
       'evidence_by_primary_eligibility':c.groupby(['screen_evidence_status','primary_eligibility_status']).size().reset_index(name='n').to_dict('records'),
       'primary_lower_content_candidates':c.loc[c.in_U&c.L,'gene'].tolist(),
       'excluded_with_published_lower_content':c.loc[~c.in_U&c.L,['gene','eligibility_reason']].to_dict('records'),
       'primary_survival_fdr_lt005_lower_content':int((c.in_U&c.L&(c.survival_fdr<.05)).sum()),
       'guide_statistics_semantics':'S8 q threshold invalidated by saved integer coercion; selected S7 continuous probabilities provide bounded support. Missing guide probabilities unknown; original zero rules are historical diagnostics.'})
    plots(c,gp,results,old['analyses'])
    print(json.dumps({'revision':2,'primary_unchanged':True,'all_candidates':len(c),'screen_status':c.screen_evidence_status.value_counts().to_dict()}))

def eligibility_flow_text(c):
    total=len(c)
    common=int(c.in_U0.sum())
    primary=int(c.in_U.sum())
    coverage_losses=int((~c.in_U0).sum())
    expression_losses=int((c.in_U0&~c.in_U).sum())
    lower=c.published_content_label.eq('Decrease')
    boxes=[(0,f'{total} source\ncandidates','Survival beta >0;\nsource FDR <0.10'),
           (3.55,f'{common} with common\ncoverage',f'{coverage_losses} excluded by coverage'),
           (7.1,f'{primary} primary\neligible',f'{expression_losses} further expression exclusions')]
    footer=f'Published lower-content labels: {int(lower.sum())} among all {total}; {int((c.in_U&lower).sum())} primary eligible; NUPL1 excluded only by expression mapping.'
    return boxes,footer


def plots(c,gp,results,old):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    folder=ROOT/'manuscript/figures';folder.mkdir(exist_ok=True)
    def save(fig,name):
        fig.savefig(folder/f'{name}.png',bbox_inches='tight');fig.savefig(folder/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,3.1));ax.set_xlim(-.2,10.1);ax.set_ylim(0,3);ax.axis('off')
    flow_boxes,flow_footer=eligibility_flow_text(c)
    for x,title,detail in flow_boxes:
        ax.add_patch(FancyBboxPatch((x,.85),2.75,1.35,boxstyle='round,pad=.10',facecolor='#ecf2f5',edgecolor='#48657c'))
        ax.text(x+1.375,1.75,title,ha='center',va='center',fontsize=12,weight='bold')
        ax.text(x+1.375,1.2,detail,ha='center',va='center',fontsize=8.6,wrap=True)
    for x in [2.9,6.45]:ax.annotate('',xy=(x+.5,1.5),xytext=(x,1.5),arrowprops={'arrowstyle':'->','lw':1.5})
    ax.text(5,.25,flow_footer,ha='center',fontsize=8.6)
    save(fig,'eligibility_flow')
    selected=c[c.in_U&c.L]
    fig,axes=plt.subplots(4,4,figsize=(9,8.7),sharex=True)
    for ax,(_,row) in zip(axes.flat,selected.iterrows()):
        sub=gp[gp.gene==row.gene].sort_values('guide_key')
        colors=plt.get_cmap('tab10')
        for i,(_,guide) in enumerate(sub.iterrows()):
            ax.plot([0,1],[guide.lfc_r0,guide.lfc_r1],'-o',color=colors(i),markersize=3,lw=1,alpha=.85)
        ax.axhline(0,color='#555555',lw=.6);ax.set_title(row.gene,fontsize=10,fontweight='bold')
        ax.set_xticks([0,1],['r0','r1']);ax.set_xlim(-.15,1.3)
        ax.tick_params(axis='x',labelbottom=True)
    for ax in axes.flat[len(selected):]:ax.axis('off')
    fig.supylabel('Insulin low/high guide LFC (log2); positive = lower content',x=.015)
    from matplotlib.lines import Line2D
    fig.legend([Line2D([0],[0],color=plt.get_cmap('tab10')(i),marker='o',lw=1,markersize=3) for i in range(4)],
               [f'Guide {i+1}' for i in range(4)],loc='lower center',ncol=4,fontsize=8,bbox_to_anchor=(.5,.018))
    fig.supxlabel('Within each gene, colors index the ordered guide keys in candidate_paired_guide_effects.csv',fontsize=8)
    fig.tight_layout(rect=[.025,.06,1,1]);save(fig,'paired_guide_effects')
    # All-range scatter retains outliers; neutral labels replace biological liability claims.
    p=pd.read_csv(OUT/'primary_genes.csv')
    fig,axes=plt.subplots(1,2,figsize=(9,3.7),sharex=True,sharey=True)
    for ax,rep in zip(axes,['r0','r1']):
        ax.scatter(p.survival_beta,p[f'median_lfc_{rep}'],s=4,c='#b7bec5',alpha=.35,rasterized=True)
        others=c[c.in_U&~c.published_content_label.isin(['Decrease','Increase'])]
        higher=c[c.in_U&c.published_content_label.eq('Increase')]
        for z,color,label in [(others,'#365f9b','No published directional call'),(higher,'#287b69','Published higher content'),(selected,'#b33c35','Published lower content')]:
            ax.scatter(z.survival_beta,z[f'median_lfc_{rep}'],s=18,color=color,label=label,alpha=.85)
        ax.axhline(0,color='#777777',lw=.6);ax.axvline(0,color='#777777',lw=.6)
        ax.set_title(f'Insulin replicate {rep}');ax.set_xlabel('Cytokine survival beta')
    axes[0].set_ylabel('Median guide LFC (log2)')
    fig.legend(*axes[1].get_legend_handles_labels(),loc='lower center',ncol=3,fontsize=8,bbox_to_anchor=(.5,-.02))
    fig.tight_layout(rect=[0,.08,1,1]);save(fig,'signed_effects')
    colors={'primary':'#203c50','no_expression':'#8c571d','expression_ge1':'#8c571d','expression_ge100':'#8c571d','survival_fdr005':'#436f86','survival_harming':'#73508b','positive_vs_negative':'#73508b'}
    labels={'primary':'Primary catalog','no_expression':'Universe: no expression filter','expression_ge1':'Universe: expression >=1','expression_ge100':'Universe: expression >=100','survival_fdr005':'Exposure threshold: FDR <0.05','survival_harming':'Comparator question: harming vs other','positive_vs_negative':'Comparator question: favoring vs harming'}
    fig,ax=plt.subplots(figsize=(9,4.6))
    for y,r in enumerate(results):
        ax.errorbar(r['conditional_mle'],y,xerr=[[r['conditional_mle']-r['ci_lower']],[r['ci_upper']-r['conditional_mle']]],fmt='o',color=colors[r['name']],capsize=3)
        ax.text(3.2,y,f"{r['conditional_mle']:.2f} ({r['ci_lower']:.2f}-{r['ci_upper']:.2f})",va='center',fontsize=9)
    ax.set_yticks(range(len(results)),[labels[r['name']] for r in results]);ax.invert_yaxis();ax.axvline(1,color='#777777',ls='--',lw=1)
    ax.set_xlim(.15,5.2);ax.set_xticks([.5,1,1.5,2,2.5,3]);ax.set_xlabel('Conditional odds ratio and exact 95% interval')
    fig.tight_layout();save(fig,'sensitivity_forest')

if __name__=='__main__':run()
