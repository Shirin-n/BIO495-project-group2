#!/usr/bin/env python3
"""Descriptive figures and complete reporting of frozen-v1 association results."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from prepare_microbiome_metrics import sha
from run_association_tests import PRIMARY

LABELS={PRIMARY:'Primary',
 'sensitivity_strict_all_mitotypes':'All mitotypes',
 'sensitivity_broad_mito_agree':'Broader health criteria',
 'sensitivity_broad_all_mitotypes':'Broader health + all mitotypes',
 'sensitivity_strict_mito_with_unassigned':'Include unassigned OTUs',
 'primary_depth_2000':'2,000 reads','primary_depth_5000':'5,000 reads','primary_depth_10000':'10,000 reads',
 'primary_min_2_colonies':'At least 2 recorded colonies','primary_min_3_colonies':'At least 3 recorded colonies',
 'primary_busco_ge_90':'BUSCO completeness >=90%',
 'primary_single_expedition':'Single expedition (E1)',
 'exposure_any_isoform_canonical':'Any-isoform canonical count',
 'exposure_nacht':'All NACHT loci','exposure_nb_arc':'NB-ARC loci',
 'exposure_canonical_per_10000_genes':'Canonical loci / 10,000 genes',
 'unrarefied_shannon':'Unrarefied Shannon'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    root=parser.parse_args().root;out=root/'results/association_v1';figdir=root/'results/figures'
    validation=json.loads((out/'validation.json').read_text());assert validation['status']=='complete'
    for name,digest in validation['products_sha256'].items():assert sha(root/name)==digest,name
    audit=json.loads((out/'independent_audit.json').read_text());assert audit['status']=='passed'
    assert audit['primary_results_sha256']==sha(out/'primary_results.tsv')
    results=pd.read_csv(out/'all_results.tsv',sep='\t')
    primary=results[results.scenario==PRIMARY].set_index('endpoint')
    species=pd.read_csv(out/'primary_species_data.tsv',sep='\t').sort_values('species')
    pairs=pd.read_csv(out/'primary_pairwise_data.tsv',sep='\t')
    loo=pd.read_csv(out/'leave_one_species_out.tsv',sep='\t')
    metrics=json.loads((out/'summary_metrics.json').read_text())
    colors={s:plt.get_cmap('tab10')(i) for i,s in enumerate(species.species)}
    short=lambda s:s.split()[0][0]+'. '+s.split()[1]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    fig,axes=plt.subplots(1,2,figsize=(12,5.8))
    colonies=pd.read_csv(root/'data/processed/microbiome_metrics_v1'/PRIMARY/'colony_alpha.tsv',sep='\t')
    for row in species.itertuples():
        values=colonies.loc[colonies.canonical_species.eq(row.species),'shannon_mean'].to_numpy()
        axes[0].scatter(np.full(len(values),row.canonical),values,s=19,color='#b4bac2',alpha=.8,zorder=1)
        if len(values)>1:axes[0].vlines(row.canonical,values.min(),values.max(),color='#b4bac2',lw=1,zorder=1)
        axes[0].scatter(row.canonical,row.shannon_mean,s=60,marker='D',color=colors[row.species],edgecolors='white',linewidth=.5,zorder=3)
        offset=(-7,9) if row.species=='Caulastraea furcata' else ((-7,-17) if row.species=='Porites cylindrica' else ((6,-15) if row.species=='Acropora hyacinthus' else (6,6)))
        axes[0].annotate(f'{short(row.species)} (n={row.n_samples})',(row.canonical,row.shannon_mean),
                         xytext=offset,textcoords='offset points',ha='right' if offset[0]<0 else 'left',fontsize=8.2)
    axes[0].set_xlabel('Canonical candidate gene loci');axes[0].set_ylabel('Species-mean Shannon diversity')
    axes[0].set_xlim(16,45);axes[0].set_ylim(-.2,7.8);axes[0].set_yticks(range(7))
    axes[0].set_title('A | Candidate count vs diversity',loc='left',fontweight='bold',fontsize=12)
    axes[1].scatter(pairs.gene_log1p_distance,pairs.bray_curtis,color='#3c638b',alpha=.75,s=32)
    axes[1].set_xlabel('Absolute difference in log(1 + candidate count)')
    axes[1].set_ylabel('Between-species Bray–Curtis distance')
    axes[1].set_title('B | Gene-count vs microbiome distances',loc='left',fontweight='bold',fontsize=12)
    axes[1].set_ylim(max(0,pairs.bray_curtis.min()-.08),1.2);axes[1].set_yticks(np.arange(.4,1.01,.1))
    for ax,endpoint in zip(axes,['A','B']):
        r=primary.loc[endpoint]
        ax.text(.025,.97,f"rho = {r.effect_rho:.3f}\nexact p = {r.p_raw:.4f}; Holm p = {r.p_holm_primary:.3f}",
                transform=ax.transAxes,va='top',fontsize=10,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
        ax.grid(alpha=.15);ax.set_axisbelow(True)
    fig.suptitle('Primary associations | 8 species; 40,320 label permutations',x=.055,ha='left',fontsize=15,fontweight='bold')
    fig.text(.055,.05,'A: diamonds = species means; gray points/ranges = recorded-colony values (not confidence intervals). n = sample count.',fontsize=8.5,color='#444444')
    fig.text(.055,.018,'B: 28 displayed pairs are dependent; inference permutes 8 species, not 28 pairwise observations. Both tests are two-sided.',fontsize=8.5,color='#444444')
    fig.subplots_adjust(left=.065,right=.98,bottom=.19,top=.83,wspace=.30)
    fig.savefig(figdir/'association_primary_v1.png');fig.savefig(figdir/'association_primary_v1.pdf');plt.close(fig)

    order=list(results.scenario.drop_duplicates());y=np.arange(len(order))
    fig,axes=plt.subplots(1,2,figsize=(12.8,9.2),sharey=True)
    labels=[f"{LABELS[name]}  (n={int(results.loc[results.scenario.eq(name),'n_species'].iloc[0])})" for name in order]
    for ax,endpoint in zip(axes,['A','B']):
        rows=results[results.endpoint==endpoint].set_index('scenario').loc[order]
        ax.axvline(0,color='#d1d5db',lw=1)
        ax.axvline(primary.loc[endpoint,'effect_rho'],color='#4B5563',lw=1,ls='--')
        for i,(name,row) in enumerate(rows.iterrows()):
            is_primary=name==PRIMARY;tested=row.status=='tested'
            ax.scatter(row.effect_rho,i,marker='D' if is_primary else 'o',s=46,
                       facecolors='black' if is_primary else ('#2457A7' if tested else 'white'),
                       edgecolors='black' if is_primary else ('#2457A7' if tested else '#9CA3AF'),zorder=3)
            p=row.p_holm_primary if is_primary else row.p_raw
            ax.text(1.16,i,'NA' if not np.isfinite(p) else f'{p:.4f}',ha='center',va='center',fontsize=8.5,color='#444444')
        ax.set_xlim(-1.05,1.4);ax.set_xticks([-1,-.5,0,.5,1]);ax.set_ylim(len(order)-.5,-1)
        ax.set_xlabel('Spearman rho' if endpoint=='A' else 'Spearman Mantel rho')
        ax.set_title(endpoint+' | '+('Diversity' if endpoint=='A' else 'Composition distances'),loc='left',fontweight='bold')
        ax.text(1.16,-.72,'p value',ha='center',fontsize=9,fontweight='bold')
        ax.grid(axis='y',alpha=.12);ax.set_axisbelow(True)
    axes[0].set_yticks(y,labels,fontsize=8.8);axes[1].tick_params(labelleft=False)
    fig.suptitle('All prespecified scenarios | no result selection',x=.035,ha='left',fontweight='bold',fontsize=15)
    fig.text(.035,.045,'Primary p values are Holm-adjusted; all sensitivity p values are unadjusted and exploratory. Dashed line: primary effect.',fontsize=9,color='#444444')
    fig.text(.035,.02,'Open markers: three-species descriptive effects only, with no p value. These are effect estimates, not confidence intervals.',fontsize=9,color='#444444')
    fig.subplots_adjust(left=.31,right=.98,bottom=.12,top=.90,wspace=.12)
    fig.savefig(figdir/'association_sensitivity_v1.png');fig.savefig(figdir/'association_sensitivity_v1.pdf');plt.close(fig)

    omitted=list(species.species);fig,axes=plt.subplots(1,2,figsize=(11.5,6),sharey=True)
    for ax,endpoint in zip(axes,['A','B']):
        frame=loo[loo.endpoint==endpoint].set_index('omitted_species').loc[omitted]
        ax.scatter(frame.effect_rho,np.arange(8),color='#C8862D',s=48,zorder=3)
        ax.axvline(primary.loc[endpoint,'effect_rho'],color='black',ls='--',lw=1,label='Full primary effect')
        ax.axvline(0,color='#c7cbd1',lw=.8)
        ax.set_xlim(-1,1);ax.set_xticks([-1,-.5,0,.5,1]);ax.set_ylim(7.5,-.5)
        ax.set_xlabel('Spearman rho' if endpoint=='A' else 'Spearman Mantel rho')
        bounds=metrics['leave_one_species_out_ranges'][endpoint]
        ax.set_title(f"{endpoint} | Range: {bounds['min']:.3f} to {bounds['max']:.3f}",loc='left',fontweight='bold',fontsize=11)
        ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    axes[0].set_yticks(np.arange(8),omitted);axes[1].tick_params(labelleft=False)
    axes[0].set_ylabel('Omitted species')
    fig.suptitle('Leave-one-species-out diagnostics | 7 species per estimate',x=.04,ha='left',fontweight='bold',fontsize=14)
    fig.text(.04,.035,'Dashed lines show the full eight-species effects. Omission ranges are not confidence intervals; no additional p values computed.',fontsize=8.8,color='#444444')
    fig.subplots_adjust(left=.26,right=.98,bottom=.15,top=.83,wspace=.25)
    fig.savefig(figdir/'association_leave_one_out_v1.png');fig.savefig(figdir/'association_leave_one_out_v1.pdf');plt.close(fig)

    tested=results[results.scenario.ne(PRIMARY)&results.status.eq('tested')]
    inference='Neither primary endpoint reaches the frozen 0.05 threshold after Holm correction.' if primary.p_holm_primary.ge(.05).all() else 'At least one primary endpoint reaches the frozen 0.05 threshold after Holm correction; interpretation remains exploratory.'
    lines=['CORAL NLR–MICROBIOME ASSOCIATION ANALYSIS — FROZEN V1','',
        inference,
        'The present data do not support the predicted positive associations. This is not evidence that no biological relationship exists.',
        '', 'PRIMARY RESULTS',primary[['n_species','n_samples','effect_rho','p_raw','p_holm_primary','permutations','extreme_count']].to_string(),
        '', 'SENSITIVITY RESULTS (ALL SCHEDULED; RAW P VALUES)',
        results[results.scenario.ne(PRIMARY)][['scenario','endpoint','n_species','n_samples','effect_rho','p_raw','status']].to_string(index=False),
        '', 'LEAVE-ONE-SPECIES-OUT EFFECTS',loo[['omitted_species','endpoint','n_species','effect_rho']].to_string(index=False),
        '', 'INTERPRETATION',
        f"{len(tested)} sensitivity endpoint tests were eligible; {int(tested.p_raw.lt(.05).sum())} had unadjusted p < 0.05.",
        'Sensitivity results are not corrected across scenarios and do not constitute independent replications of the primary tests.',
        'The any-isoform canonical and BUSCO>=90% scenarios coincide with primary inputs here; the unrarefied-Shannon change leaves B unchanged.',
        'The minimum-three-colony scenario contains only three species. Its coefficients are descriptive and its p values are undefined.',
        'The broad/all-mitotypes cohort is the prelisted combined scenario, not a one-factor perturbation.',
        'The leave-one-species-out ranges describe sensitivity to individual species, not confidence intervals.',
        '', 'METHODS',
        'A: Spearman rank correlation of canonical candidate gene count with species-mean standardized Shannon diversity.',
        'B: Spearman Mantel correlation of absolute log1p gene-count differences with Bray-Curtis distances between species profiles.',
        'Both tests use absolute-tail two-sided species-label permutation p values. Mantel permutes matrix rows and columns together.',
        'Eight-species tests enumerate all 40,320 labelled permutations, including identity. Tied counts retain average ranks and permutation multiplicities.',
        'Nine-species tests use 9,999 independent PCG64 label permutations, allowing repeats, with (extreme+1)/(B+1).',
        'Seed 20261117 uses SeedSequence children assigned by sorted scenario name. Stream specifications and permutation fingerprints are retained.',
        'Integer-scaled centered average ranks give exact comparisons of absolute permutation numerators, avoiding floating-point tie decisions.',
        'Primary A/B p values receive Holm correction with a fixed two-test family. Sensitivity p values remain explicitly unadjusted.',
        'Density exposures are recomputed from exact counts/annotation totals rather than rounded display values.',
        'Species are equally weighted; the 25 samples and 28 pairwise entries are not treated as independent inferential replicates.',
        '', 'LIMITATIONS',
        'Only eight species contribute to primary inference. The analysis is exploratory and cannot exclude biologically meaningful associations.',
        'Recorded colony IDs preserve dates and do not establish independent biological colonies. Five species have fewer than three recorded colony IDs.',
        'A 1,000-read Shannon standardization does not establish saturation. All prespecified higher-depth results are retained.',
        'Candidate counts are domain-architecture predictions, not experimentally validated immune gene repertoires; 22 canonical candidates have overlapping core-model assignments.',
        'Nominal permutation p values assume exchangeable species labels; shared ancestry and geography may violate this assumption.',
        'No suitable species-level tree is available. PGLS (lambda=1) and its lambda=0 contrast were not performed; ancestry is not adjusted.',
        'The single-expedition subset is a sensitivity analysis, not geographic adjustment. Site and species coverage remain uneven.',
        'No PERMANOVA or partial Mantel test was added.',
        '', 'VALIDATION',
        'Seven boundary/reference tests passed. Independent SciPy exact absolute-Spearman and raw-distance reranking checks reproduced the primary p values.',
        'Every saved permutation null distribution reproduces its reported extreme count and p value.',
        'All 65 frozen input representations and the completed diversity/NLR preparation hashes were verified.',
        '', 'REPRODUCTION',
        'conda activate bio-project',
        'python scripts/test_association_tests.py',
        'python scripts/run_association_tests.py',
        'python scripts/audit_association_results.py',
        'python scripts/summarize_association_tests.py',
        'The existing coral_nlr_microbiome_analysis.ipynb contains the complete executed workflow.',
        '', 'NEXT STAGE',
        'Synthesize the completed analysis into the project report and presentation, including negative results and all material limitations.',
        'Any new filtering rule, additional test or new tree-based analysis must be identified transparently relative to frozen v1.',
        '', 'METHOD REFERENCES',
        'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html',
        'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html',
        'https://scikit.bio/docs/dev/generated/skbio.stats.distance.mantel.html','']
    (root/'reports/association_v1_report.txt').write_text('\n'.join(lines))
    products=[root/'reports/association_v1_report.txt',out/'independent_audit.json']
    products += [figdir/(name+ext) for name in ['association_primary_v1','association_sensitivity_v1','association_leave_one_out_v1'] for ext in ['.png','.pdf']]
    completion=dict(status='complete',primary_conclusion=inference,pgls='not_performed_no_adequate_tree',
        source_script_sha256=sha(Path(__file__)),input_validation_sha256=sha(out/'validation.json'),
        products_sha256={str(p.relative_to(root)):sha(p) for p in products})
    (out/'completion.json').write_text(json.dumps(completion,indent=2)+'\n')
    print(inference)

if __name__=='__main__':main()
