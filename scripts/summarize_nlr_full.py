#!/usr/bin/env python3
"""Audit the full v1 run; export analysis-ready exposures and descriptive figures."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd
from run_nlr_pilot import sha, fasta, write_fasta


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    root=parser.parse_args().root
    source=root/'results/nlr_full_v1'
    processed=root/'data/processed/nlr_v1'
    figures=root/'results/figures'
    processed.mkdir(parents=True,exist_ok=True)
    (source/'completion.json').unlink(missing_ok=True)
    protocol=json.loads((root/'config/analysis_protocol_v1.json').read_text())
    validation=json.loads((source/'validation.json').read_text())
    assert validation['analysis_gate']=='pass'
    assert sha(root/'scripts/run_nlr_full.py')==validation['source_script_sha256']
    assert sha(root/'scripts/run_nlr_pilot.py')==validation['classification_module_sha256']
    for path,digest in validation['output_sha256'].items():assert sha(root/path)==digest,path
    read=lambda name:pd.read_csv(source/name,sep='\t')
    summary=read('species_summary.tsv')
    proteins=read('protein_classification.tsv')
    loci=read('locus_classification.tsv')
    domains=read('domain_evidence.tsv')
    pairs=read('blast_pairs.tsv')
    assert len(summary)==9 and summary.unresolved_genes.eq(0).all()
    assert not loci.duplicated(['genome_id','gene_id']).any()
    assert not proteins.duplicated(['genome_id','protein_id']).any()
    # Independently reconstruct representative selection from the exported sequences' QC.
    eligible=proteins[proteins.internal_stop_count.eq(0)]
    expected=eligible.sort_values(['genome_id','gene_id','length','protein_id'],ascending=[True,True,False,True]).drop_duplicates(['genome_id','gene_id'])
    representatives=proteins[proteins.is_representative].copy()
    keys=lambda frame:set(zip(frame.genome_id,frame.protein_id))
    assert keys(expected)==keys(representatives)
    assert len(representatives)==summary.annotated_genes.sum()
    accepted=domains[domains.retained].copy()
    assert accepted.model_coverage.ge(.5).all()
    assert accepted.sequence_score.ge(accepted.sequence_ga).all()
    assert accepted.domain_score.ge(accepted.domain_ga).all()
    # Independent architecture reconstruction, including negative proteins.
    nach=accepted[accepted.model.eq(protocol['nlr']['core_model'])]
    lrr=accepted[accepted.model.isin(protocol['nlr']['lrr_models'])]
    min_core_end=nach.groupby(['genome_id','protein_id']).ali_to.min()
    max_lrr_start=lrr.groupby(['genome_id','protein_id']).ali_from.max()
    joined=pd.concat([min_core_end.rename('core_end'),max_lrr_start.rename('lrr_start')],axis=1).dropna()
    expected_canonical=set(joined[joined.lrr_start>joined.core_end].index)
    assert expected_canonical==keys(proteins[proteins.canonical])
    canonical=representatives[representatives.canonical].copy()
    canonical['candidate_id']=canonical.genome_id+'__'+canonical.protein_id
    assert canonical.candidate_id.is_unique
    for row in summary.itertuples():
        pr=representatives[representatives.species.eq(row.species)]
        assert len(pr)==row.annotated_genes
        for col in ['canonical','nacht','nb_arc']:assert int(pr[col].sum())==getattr(row,col)
        any_loci=eligible[eligible.species.eq(row.species)&eligible.canonical].gene_id.nunique()
        assert any_loci==row.any_isoform_canonical

    # Scope flags are attached after sequence classification; outcomes are not inspected.
    summary['primary_eligible']=summary.species.isin(protocol['primary_species'])
    summary['scope']=np.where(summary.primary_eligible,'primary_and_sensitivity','sensitivity_only')
    assert summary.primary_eligible.sum()==8
    summary['any_isoform_delta']=summary.any_isoform_canonical-summary.canonical
    summary['busco_ge_90']=summary.busco_complete_percent.ge(90)
    for col in ['overlapping_core_models','mixed_core_models']:
        counts=canonical.groupby('species')[col].sum()
        summary['canonical_'+col]=summary.species.map(counts).fillna(0).astype(int)
    summary['qc_status']='verified_complete'
    summary.to_csv(processed/'species_exposures.tsv',sep='\t',index=False)
    summary[summary.primary_eligible].to_csv(processed/'primary_species_exposures.tsv',sep='\t',index=False)
    summary.to_csv(processed/'sensitivity_all_species_exposures.tsv',sep='\t',index=False)

    best=pairs.sort_values(['query_id','best_bitscore','reference_id'],ascending=[True,False,True]).drop_duplicates('query_id')
    canonical=canonical.merge(best,left_on='candidate_id',right_on='query_id',how='left',validate='one_to_one')
    canonical['has_reference_hit']=canonical.reference_id.notna()
    assert int(canonical.has_reference_hit.sum())==int(summary.canonical_with_reference_hit.sum())
    canonical['review_flag']=canonical.multiple_core_domains | canonical.missing_gtf_start_codon | canonical.missing_gtf_stop_codon | canonical.ambiguous_residues.gt(0)
    canonical.to_csv(processed/'canonical_candidates.tsv',sep='\t',index=False)
    canonical[canonical.review_flag].to_csv(source/'flagged_candidates.tsv',sep='\t',index=False)
    query_sequences=fasta(root/'data/interim/nlr_full_v1/blast_queries.faa')
    write_fasta(processed/'canonical_candidates.faa',{k:query_sequences[k] for k in canonical.candidate_id})
    assert len(fasta(processed/'canonical_candidates.faa'))==len(canonical)

    flagged_keys=canonical[canonical.multiple_core_domains][['genome_id','protein_id']]
    complex_domains=accepted[accepted.model.isin([protocol['nlr']['core_model'],protocol['nlr']['broad_core_model']])].merge(flagged_keys,on=['genome_id','protein_id'],validate='many_to_one')
    complex_domains.to_csv(source/'complex_core_domain_review.tsv',sep='\t',index=False)
    discordant=loci[loci.any_isoform_canonical & ~loci.canonical]
    alternatives=proteins.merge(discordant[['genome_id','gene_id']],on=['genome_id','gene_id'],validate='many_to_one')
    alternatives.to_csv(source/'isoform_sensitivity_review.tsv',sep='\t',index=False)
    assert len(discordant)==int(summary.any_isoform_delta.sum())

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    order=summary.sort_values('species').reset_index(drop=True)
    labels=[r.species+(' *' if not r.primary_eligible else '') for r in order.itertuples()]
    y=np.arange(len(order));colors=np.where(order.primary_eligible,'#2457A7','#CA8C38')
    fig,axes=plt.subplots(1,2,figsize=(11.8,6.6),sharey=True)
    for ax,col,title in [(axes[0],'canonical','Canonical candidate loci'),(axes[1],'canonical_per_10000_genes','Candidates per 10,000 annotated genes')]:
        bars=ax.barh(y,order[col],color=colors,height=.6)
        ax.bar_label(bars,fmt='%.0f' if col=='canonical' else '%.2f',padding=5,fontsize=10)
        ax.set_xlim(0,order[col].max()*1.25)
        ax.set_title(title,loc='left',fontweight='bold',fontsize=11,pad=14)
        ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True)
        ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    axes[0].set_yticks(y,labels);axes[0].invert_yaxis();axes[1].tick_params(labelleft=False)
    fig.suptitle('Nine-species NLR identification | frozen rule v1',x=.04,ha='left',fontsize=16,fontweight='bold')
    fig.legend(handles=[Patch(facecolor='#2457A7',label='Eight primary-analysis species'),Patch(facecolor='#CA8C38',label='* Sensitivity only')],loc='lower left',bbox_to_anchor=(.04,.025),ncol=2,frameon=False,fontsize=9)
    fig.text(.04,.016,'Unique GTF loci with NACHT + downstream LRR in the representative protein; architecture candidates, not validated function.',fontsize=8.7,color='#444444')
    fig.subplots_adjust(left=.245,right=.98,bottom=.16,top=.85,wspace=.28)
    fig.savefig(figures/'nlr_full_counts.png');fig.savefig(figures/'nlr_full_counts.pdf');plt.close(fig)

    # A count heatmap shows overlapping categories without suggesting an additive total.
    matrix=order[['nacht','canonical','nb_arc','any_isoform_canonical']].to_numpy()
    fig,ax=plt.subplots(figsize=(9.4,6.5))
    heat=ax.imshow(np.log1p(matrix),cmap='Blues',aspect='auto',vmin=0)
    for i in range(len(order)):
        for j in range(4):
            ax.text(j,i,str(matrix[i,j]),ha='center',va='center',fontsize=12,
                    color='white' if np.log1p(matrix[i,j])>np.log1p(matrix.max())*.66 else '#152E50')
    ax.set_yticks(y,labels);ax.set_xticks(range(4),['Any NACHT\n(representative)','NACHT + LRR\n(representative)','NB-ARC\n(representative)','NACHT + LRR\n(any eligible isoform)'])
    ax.xaxis.tick_top();ax.tick_params(length=0,pad=10)
    for spine in ax.spines.values():spine.set_visible(False)
    fig.suptitle('Architecture categories and isoform sensitivity',x=.04,ha='left',fontsize=15,fontweight='bold')
    fig.text(.04,.04,'Numbers are unique gene loci; colour uses log(1 + count). Categories overlap and must not be summed.\n* Echinopora horrida enters sensitivity cohorts only.',fontsize=9,color='#444444')
    fig.subplots_adjust(left=.30,right=.98,top=.80,bottom=.14)
    fig.savefig(figures/'nlr_full_architecture_counts.png');fig.savefig(figures/'nlr_full_architecture_counts.pdf');plt.close(fig)

    metrics=dict(
                 species=9,primary_species=8,protein_isoforms=int(summary.protein_isoforms.sum()),
                 annotated_gene_loci=int(summary.annotated_genes.sum()),canonical_candidates=len(canonical),
                 primary_canonical_candidates=int(summary.loc[summary.primary_eligible,'canonical'].sum()),
                 nacht_loci=int(summary.nacht.sum()),nb_arc_loci=int(summary.nb_arc.sum()),
                 any_isoform_canonical_loci=int(summary.any_isoform_canonical.sum()),
                 isoform_sensitive_loci=len(discordant),
                 canonical_overlapping_core_model_flags=int(canonical.overlapping_core_models.sum()),
                 canonical_multiple_core_flags=int(canonical.multiple_core_domains.sum()),
                 canonical_ambiguous_residue_flags=int(canonical.ambiguous_residues.gt(0).sum()),
                 canonical_partial_gtf_flags=int((canonical.missing_gtf_start_codon | canonical.missing_gtf_stop_codon).sum()),
                 canonical_with_reference_hit=int(canonical.has_reference_hit.sum()),
                 mapping_percent=100.0,unresolved_loci=int(summary.unresolved_genes.sum()),
                 source_busco_complete_min=float(summary.busco_complete_percent.min()),
                 pilot_regression='all five evidence tables identical for the three pilot species',
                 independent_representative_and_architecture_checks='passed',association_tests_run=False)
    (source/'summary_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    brief=summary[['species','scope','annotated_genes','canonical','nacht','nb_arc','any_isoform_canonical','canonical_per_10000_genes','canonical_overlapping_core_models','busco_complete_percent']]
    lines=['NINE-SPECIES NLR IDENTIFICATION — V1 COMPLETION REPORT','',
           'The frozen v1 architecture rule was applied to all nine matched Metz v1.0.0 proteomes.',
           'Eight species enter the primary analysis; Echinopora horrida is sensitivity-only.',
           'No microbiome diversity metric, association test or phylogenetic model was computed in this stage.',
           '', 'COUNTS',brief.to_string(index=False),'', 'VALIDATION',json.dumps(metrics,indent=2),
           '', 'INTERPRETATION AND REVIEW',
           'Counts represent unique GTF gene loci, not isoforms or validated immune functions.',
           'Primary representatives are the longest eligible isoform after removal of a terminal stop, with lexical ID tie-break.',
           'Primary architecture requires accepted NACHT and at least one downstream accepted LRR in the same representative.',
           'Model-specific GA thresholds, >=50% HMM coverage and repeat-overlap rules are unchanged from the pilot.',
           'All-isoform counts retain the same sequence-QC rule. Changed loci are listed in isoform_sensitivity_review.tsv.',
           'Overlapping NACHT/NB-ARC assignments are flagged, retained under v1 and counted only once per category.',
           'Multiple core-model matches need not imply multiple physically distinct nucleotide-binding domains.',
           'Missing GTF start/stop features and ambiguous residues are review flags, not automatic exclusions under v1.',
           'NB-ARC and NACHT categories overlap; neither their sum nor a reference BLAST hit establishes total NLR repertoire.',
           'BLAST hits to the five-sequence reference panel support homology review; no identity or human-match filter is applied.',
           'The candidate table includes best-hit identity and separate union query/subject coverage. Partial matches remain partial evidence.',
           'All five reference controls and their BLAST self-matches passed. No validated coral positive control is available.',
           'BUSCO values are the original source-reported, OCR-checked annotation scores, not newly run assessments.',
           'Genes lacking a resolved primary representative would be unresolved, not zero; none were encountered in this run.',
           '', 'PRIMARY DELIVERABLES',
           'data/processed/nlr_v1/species_exposures.tsv: nine-species exposure table with scope and review flags.',
           'data/processed/nlr_v1/primary_species_exposures.tsv: eight-species primary exposure table.',
           'data/processed/nlr_v1/sensitivity_all_species_exposures.tsv: nine-species sensitivity exposure table.',
           'data/processed/nlr_v1/canonical_candidates.tsv and .faa: per-candidate evidence and sequences.',
           'results/nlr_full_v1/: full protein, locus, domain and BLAST evidence, review tables and validation records.',
           'coral_nlr_microbiome_analysis.ipynb: the existing English notebook, extended with full-run methods and outputs.',
           '', 'REPRODUCTION',
           'conda activate bio-project',
           'python scripts/test_nlr_pilot.py',
           'python scripts/run_nlr_full.py',
           'python scripts/summarize_nlr_full.py',
           'HMMER outputs are reused only when input, command, tool version and output hashes agree.',
           '', 'NEXT STAGE',
           'Compute frozen repeated-rarefaction Shannon diversity and species-level composition summaries.',
           'Retain colony identifiers and equal colony/species weighting as defined by v1.',
           'Check compatibility with these frozen NLR exposures before proceeding to the planned association tests.',
           'The host phylogeny and coral-specific functional reference limitations remain unchanged.','']
    (root/'reports/nlr_full_v1_report.txt').write_text('\n'.join(lines))
    products=list(processed.glob('*'))+[source/name for name in ['flagged_candidates.tsv','complex_core_domain_review.tsv','isoform_sensitivity_review.tsv','summary_metrics.json']]
    products += [figures/name for name in ['nlr_full_counts.png','nlr_full_counts.pdf','nlr_full_architecture_counts.png','nlr_full_architecture_counts.pdf']]
    products += [root/'reports/nlr_full_v1_report.txt']
    completion=dict(status='complete',metrics=metrics,source_script_sha256=sha(Path(__file__)),
                    full_validation_sha256=sha(source/'validation.json'),
                    products_sha256={str(f.relative_to(root)):sha(f) for f in products})
    (source/'completion.json').write_text(json.dumps(completion,indent=2)+'\n')
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':main()
