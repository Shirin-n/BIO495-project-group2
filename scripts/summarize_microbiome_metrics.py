#!/usr/bin/env python3
"""Figures and descriptive report for frozen-v1 microbiome preparation."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from prepare_microbiome_metrics import sha,save

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    root=parser.parse_args().root
    out=root/'results/microbiome_metrics_v1';base=root/'data/processed/microbiome_metrics_v1'
    primary=base/'primary_strict_mito_agree';figdir=root/'results/figures'
    v=json.loads((out/'validation.json').read_text());assert v['status']=='complete'
    for name,digest in v['products_sha256'].items():assert sha(root/name)==digest,name
    sample=pd.read_csv(primary/'sample_alpha.tsv',sep='\t').set_index('sample_id')
    species=pd.read_csv(primary/'species_alpha.tsv',sep='\t').set_index('species').sort_index()
    scenarios=pd.read_csv(out/'scenario_coverage.tsv',sep='\t')
    diag=pd.read_csv(out/'rarefaction_diagnostics.tsv',sep='\t')
    pinfo=json.loads((out/'pcoa_diagnostics.json').read_text())
    species_order=list(species.index);palette=plt.get_cmap('tab10')
    colors={name:palette(i) for i,name in enumerate(species_order)}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    fig,ax=plt.subplots(figsize=(10.5,6.1))
    for y,name in enumerate(species_order):
        points=sample[sample.species==name].sort_index()
        jitter=np.linspace(-.16,.16,len(points)) if len(points)>1 else np.array([0])
        ax.scatter(points.shannon_mean,y+jitter,color=colors[name],s=35,alpha=.85,zorder=3)
        row=species.loc[name]
        ax.errorbar(row.shannon_mean,y,xerr=[[row.shannon_mean-row.shannon_mc_p05],[row.shannon_mc_p95-row.shannon_mean]],fmt='D',ms=6,color='black',capsize=3,zorder=4)
    labels=[f"{name}  (n={int(species.loc[name,'n_recorded_colonies'])})" for name in species_order]
    ax.set_yticks(range(8),labels);ax.invert_yaxis();ax.set_xlim(left=0)
    ax.set_xlabel('Shannon diversity (natural log; mean of 100 subsamples at 1,000 reads)')
    ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True)
    ax.set_title('Primary cohort | standardized Shannon diversity',loc='left',fontweight='bold',fontsize=14,pad=16)
    ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color='#6B7280',label='Sample mean'),Line2D([],[],marker='D',linestyle='',color='black',label='Species mean; bars: 5–95% Monte Carlo range')],loc='upper right',frameon=False,fontsize=9)
    fig.text(.04,.025,'n = recorded colony IDs, not independently verified biological replicates. Monte Carlo ranges are not biological confidence intervals.',fontsize=8.4,color='#444444')
    fig.subplots_adjust(left=.33,right=.98,bottom=.15,top=.88)
    fig.savefig(figdir/'microbiome_shannon_v1.png');fig.savefig(figdir/'microbiome_shannon_v1.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,6.8))
    species_short={s:s.split()[0][0]+'. '+s.split()[1] for s in species_order}
    expeditions=sorted(sample.expedition.unique());markers=dict(zip(expeditions,['o','^','s','P']))
    for ax,level,title in [(axes[0],'sample','Samples (descriptive EDA)'),(axes[1],'species','Equal-colony species profiles')]:
        coords=pd.read_csv(primary/f'{level}_pcoa_coordinates.tsv',sep='\t',index_col=0)
        inf=pinfo['primary_strict_mito_agree/'+level]
        for key,row in coords.iterrows():
            sp=sample.loc[key,'species'] if level=='sample' else key
            marker=markers[sample.loc[key,'expedition']] if level=='sample' else 'o'
            ax.scatter(row.PCoA1,row.PCoA2,color=colors[sp],marker=marker,s=48 if level=='sample' else 65,edgecolors='white',linewidth=.45,zorder=3)
            if level=='species':ax.annotate(species_short[sp],(row.PCoA1,row.PCoA2),xytext=(6,5),textcoords='offset points',fontsize=8.3)
        ax.margins(.24);ax.axhline(0,color='#e5e7eb',lw=.8);ax.axvline(0,color='#e5e7eb',lw=.8)
        ax.set_xlabel(f"PCoA1 ({inf['axis1_positive_fraction']:.1%} of positive inertia)")
        ax.set_ylabel(f"PCoA2 ({inf['axis2_positive_fraction']:.1%} of positive inertia)")
        ax.set_title(title,loc='left',fontweight='bold',fontsize=11)
        ax.text(.02,.02,f"Negative eigenvalues: {inf['negative_eigenvalues']}\nNegative |inertia| share: {abs(inf['negative_fraction_absolute_inertia'])*100:.3f}%",transform=ax.transAxes,fontsize=8,color='#444444',va='bottom')
    handles=[Line2D([],[],marker='o',linestyle='',color=colors[s],label=species_short[s]) for s in species_order]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=4,frameon=False,fontsize=9)
    fig.suptitle('Bray–Curtis PCoA | original distances retained',x=.06,ha='left',fontweight='bold',fontsize=15)
    fig.text(.06,.015,'Sample markers: circle = E1; triangle = E10; square = E7 (see metadata). No distance correction or group-significance test applied.',fontsize=8.4,color='#444444')
    # Derive marker labels from the actual lexicographic expedition order.
    fig.texts[-1].set_text('Sample markers: '+', '.join(f"{ {'o':'circle','^':'triangle','s':'square','P':'plus'}[markers[e]]} = {e.split('_')[0]}" for e in expeditions)+'. Original distances retained; no group-significance test.')
    fig.subplots_adjust(left=.075,right=.96,bottom=.24,top=.87,wspace=.27)
    fig.savefig(figdir/'microbiome_pcoa_v1.png');fig.savefig(figdir/'microbiome_pcoa_v1.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,6))
    for sid,rows in diag.groupby('sample_id'):
        rows=rows.sort_values('depth');color=colors[rows.species.iloc[0]]
        axes[0].plot(rows.depth,rows.shannon_mean,color=color,alpha=.65,lw=1.2)
        axes[1].plot(rows.depth,rows.observed_otus_mean,color=color,alpha=.65,lw=1.2)
    for ax in axes:
        ax.set_xscale('log');ax.set_xticks([100,250,500,1000,2000,5000,10000],[100,250,500,'1k','2k','5k','10k'])
        ax.set_xlabel('Reads per subsample');ax.axvline(1000,color='black',ls='--',lw=1)
        ax.grid(alpha=.15);ax.set_axisbelow(True)
    axes[0].set_ylabel('Mean Shannon diversity');axes[1].set_ylabel('Mean observed OTUs')
    axes[0].set_title('Shannon depth sensitivity',loc='left',fontweight='bold')
    axes[1].set_title('Observed richness is depth-dependent',loc='left',fontweight='bold')
    fig.suptitle('Rarefaction diagnostics | each line is one primary sample',x=.06,ha='left',fontweight='bold',fontsize=14)
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=4,frameon=False,fontsize=9)
    fig.text(.06,.012,'100 subsamples per point; curves stop at the available library size. Dashed line: frozen 1,000-read depth. No saturation claim.',fontsize=8.6,color='#444444')
    fig.subplots_adjust(left=.07,right=.98,bottom=.25,top=.86,wspace=.25)
    fig.savefig(figdir/'microbiome_rarefaction_v1.png');fig.savefig(figdir/'microbiome_rarefaction_v1.pdf');plt.close(fig)

    family=pd.read_csv(primary/'species_family_profiles.tsv',sep='\t',index_col=0).loc[species_order]
    top=family.mean(axis=0).sort_values(ascending=False).head(10).index.tolist()
    selected=family[top].copy();selected['Other families']=1-selected.sum(axis=1)
    assert np.allclose(selected.sum(axis=1),1)
    save(selected,out/'family_composition_plot_data.tsv',index=True)
    fig,ax=plt.subplots(figsize=(12,6.5));left=np.zeros(8);family_colors=plt.get_cmap('tab20')
    for i,col in enumerate(selected):
        ax.barh(np.arange(8),selected[col]*100,left=left*100,color=family_colors(i),label=col,height=.66)
        left+=selected[col].to_numpy()
    ax.set_yticks(range(8),species_order);ax.invert_yaxis();ax.set_xlim(0,100)
    ax.set_xlabel('Mean relative abundance (%)');ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
    ax.set_title('Primary species composition | source family annotations',loc='left',fontweight='bold',fontsize=13,pad=15)
    ax.legend(loc='upper left',bbox_to_anchor=(1.01,1.02),frameon=False,fontsize=8.5)
    fig.text(.04,.025,'Samples normalized before equal-colony averaging. Top 10 family labels by equal-species mean; source taxonomy retained.',fontsize=8.7,color='#444444')
    fig.subplots_adjust(left=.25,right=.68,bottom=.14,top=.87)
    fig.savefig(figdir/'microbiome_composition_v1.png');fig.savefig(figdir/'microbiome_composition_v1.pdf');plt.close(fig)

    wide=diag.pivot(index='sample_id',columns='depth',values='shannon_mean')
    delta=(wide[5000]-wide[1000]).dropna();recovery=diag[diag.depth==1000].fraction_source_observed_otus
    sample_pc=pinfo['primary_strict_mito_agree/sample'];species_pc=pinfo['primary_strict_mito_agree/species']
    summary=dict(primary_samples=len(sample),primary_species=len(species),recorded_colony_ids=int(species.n_recorded_colonies.sum()),
                 primary_depth=1000,subsamples_per_sample=100,shannon_log_base='e',
                 smallest_species_mean_shannon=float(species.shannon_mean.min()),largest_species_mean_shannon=float(species.shannon_mean.max()),
                 species_with_single_recorded_colony=species.index[species.n_recorded_colonies.eq(1)].tolist(),
                 species_with_at_least_3_recorded_colonies=int(species.n_recorded_colonies.ge(3).sum()),
                 scenarios_prepared=len(scenarios),shared_samples_at_1000_and_5000=len(delta),
                 matched_sample_shannon_5000_minus_1000_median=float(delta.median()),
                 matched_sample_shannon_5000_minus_1000_max=float(delta.max()),
                 median_fraction_source_observed_otus_at_1000=float(recovery.median()),
                 otu_fraction_interpretation='fraction of OTUs observed in the source library, not estimated population coverage',
                 sample_pcoa=sample_pc,species_pcoa=species_pc,
                 primary_species_analysis_ready=True,association_tests_run=False,
                 depth_interpretation='The frozen 1,000-read normalization is not evidence of saturation; retain all planned depth sensitivities.',
                 limitation='Recorded colony identifiers preserve source dates and do not establish independent biological colonies.')
    (out/'summary_metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
    report=['MICROBIOME DIVERSITY AND COMPOSITION PREPARATION — FROZEN V1','',
            'Completed the prespecified preprocessing stage. No NLR association, Mantel, PERMANOVA or phylogenetic test was performed.',
            'The frozen protocol, input data and NLR candidate exposures remain unchanged.','',
            'PRIMARY SPECIES METRICS',species[['n_samples','n_recorded_colonies','shannon_mean','observed_otus_mean','inverse_simpson_mean']].to_string(),'',
            'SCENARIO COVERAGE',scenarios.to_string(index=False),'','DIAGNOSTICS',json.dumps(summary,indent=2),'',
            'METHODS AND INTERPRETATION',
            'Shannon uses natural logs, averaging 100 hypergeometric subsamples without replacement at the specified read depth.',
            'Observed OTUs and inverse Simpson use the same draws; unrarefied Shannon is exported as a sensitivity endpoint.',
            'PCG64 and SeedSequence streams are assigned to the sorted union of frozen sample IDs, then sorted diagnostic depths.',
            'Positive OTUs are sorted as strings. Identical sample/count vectors reuse draws across cohorts. NumPy is pinned to 2.4.6.',
            'Metrics and profiles are first averaged within each species/recorded-colony ID, then equally across colonies within species.',
            'Recorded colony IDs retain dates. Do not strip dates, infer independent biological replication, or pool raw reads across species.',
            'Monte Carlo SDs/percentile ranges describe subsampling variability, not biological confidence intervals.',
            'Species with one recorded colony have undefined colony-level SD and no fabricated biological interval.',
            'Bray-Curtis uses unrarefied sample relative abundances, with equal-colony species averaging. No prevalence filter was added.',
            'PCoA uses the original distance matrices, retains all signed eigenvalues, and plots positive coordinates. Axis fractions use positive inertia.',
            'The source family taxonomy is unchanged. Top family labels are selected by equal-species mean solely for descriptive plotting.',
            'The broad-all-mitotypes cohort is the prelisted combined health/mitotype scenario, not a single-factor perturbation.',
            'Depth scenarios vary the primary cohort only; they are not crossed with all other sensitivity choices.',
            'The single-expedition scenario selects the expedition with the most eligible species (at least five), then lexicographic tie-break.',
            'The minimum-three-colony scenario retains only three species and is descriptive-only under the frozen five-species safeguard.',
            'Higher depth retains all eight primary species here but changes sample coverage and, for some species, colony coverage.',
            'Rarefaction curves do not establish saturation at 1,000 reads. Observed OTU counts are depth-specific, not complete biological richness.',
            'Retain v1 as the standardized exploratory analysis and report the planned depth sensitivities. No depth was chosen using association results.',
            '', 'DELIVERABLES',
            'data/processed/microbiome_metrics_v1/<scenario>/: sample/colony/species alpha, draw metrics, profiles, Bray-Curtis matrices, PCoA and integrated species data.',
            'results/microbiome_metrics_v1/: scenario membership, RNG streams, depth diagnostics, signed-eigenvalue diagnostics and validation hashes.',
            'coral_nlr_microbiome_analysis.ipynb: existing English notebook extended with these methods, tables and figures.',
            '', 'NEXT STAGE',
            'Run the frozen species-level A/B association tests, exact primary permutations and two-test Holm correction.',
            'Report all planned sensitivity effects and exclusions; do not select the most favorable depth, cohort or exposure.',
            'PGLS remains conditional on a suitable species-level host tree. The current source tree is not ready for this analysis.',
            '', 'REPRODUCTION',
            'conda activate bio-project',
            'python scripts/test_microbiome_metrics.py',
            'python scripts/prepare_microbiome_metrics.py',
            'python scripts/summarize_microbiome_metrics.py','',
            'METHOD REFERENCES',
            'https://numpy.org/doc/2.1/reference/random/generated/numpy.random.Generator.multivariate_hypergeometric.html',
            'https://numpy.org/doc/stable/reference/random/parallel.html',
            'https://scikit.bio/docs/latest/generated/skbio.stats.ordination.pcoa.html','']
    (root/'reports/microbiome_metrics_v1_report.txt').write_text('\n'.join(report))
    products=[out/'summary_metrics.json',out/'family_composition_plot_data.tsv',root/'reports/microbiome_metrics_v1_report.txt']
    products += [figdir/(stem+suffix) for stem in ['microbiome_shannon_v1','microbiome_pcoa_v1','microbiome_rarefaction_v1','microbiome_composition_v1'] for suffix in ['.png','.pdf']]
    (out/'completion.json').write_text(json.dumps(dict(status='complete',source_script_sha256=sha(Path(__file__)),
        input_validation_sha256=sha(out/'validation.json'),products_sha256={str(p.relative_to(root)):sha(p) for p in products}),indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
