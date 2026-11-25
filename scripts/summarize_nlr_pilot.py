#!/usr/bin/env python3
"""Generate figures and a technical review from validated pilot outputs."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
root=parser.parse_args().root
out=root/'results/nlr_pilot_v1'
figdir=root/'results/figures'
figdir.mkdir(exist_ok=True)
v=json.loads((out/'validation.json').read_text())
assert v['pilot_gate']=='pass'
for name,digest in v['output_sha256'].items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
s=pd.read_csv(out/'species_summary.tsv',sep='\t')
p=pd.read_csv(out/'protein_classification.tsv',sep='\t')
d=pd.read_csv(out/'domain_evidence.tsv',sep='\t')
c=pd.read_csv(out/'control_classification.tsv',sep='\t')
cd=pd.read_csv(out/'control_domains.tsv',sep='\t')
canonical=p[p.is_representative & p.canonical].copy()
complex_candidates=canonical[canonical.multiple_core_domains].copy()
review=[]
for r in complex_candidates.itertuples():
    hh=d[(d.genome_id==r.genome_id)&(d.protein_id==r.protein_id)&d.retained]
    for h in hh[hh.model.isin(['PF05729.20','PF00931.30'])].itertuples():
        review.append(dict(species=r.species,gene_id=r.gene_id,protein_id=r.protein_id,
                           model=h.model,domain=h.domain_name,start=h.ali_from,end=h.ali_to,
                           bit_score=h.domain_score,model_coverage=h.model_coverage,
                           interpretation='Overlapping NACHT/NB-ARC model evidence; retained under frozen v1, not two independent cores.'))
pd.DataFrame(review).to_csv(out/'overlapping_core_review.tsv',sep='\t',index=False)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
labels=[x.replace(' ','\n',1) for x in s.species]
fig,axes=plt.subplots(1,2,figsize=(11.7,5),gridspec_kw={'width_ratios':[1.45,1]})
x=np.arange(len(s));width=.24
for offset,col,label,color in [(-1,'nacht','Any NACHT (representative)','#72A6AD'),(0,'canonical','NACHT + downstream LRR','#2457A7'),(1,'nb_arc','NB-ARC (separate category)','#D48B46')]:
    bars=axes[0].bar(x+offset*width,s[col],width,label=label,color=color)
    axes[0].bar_label(bars,padding=3,fontsize=10)
axes[0].set_xticks(x,labels);axes[0].set_ylabel('Unique gene loci');axes[0].set_ylim(0,125)
axes[0].set_title('Candidate counts',loc='left',fontweight='bold')
axes[0].legend(frameon=False,fontsize=8.3,loc='upper left')
bars=axes[1].bar(x,s.canonical_per_10000_genes,color='#2457A7',width=.57)
axes[1].bar_label(bars,fmt='%.2f',padding=4)
axes[1].set_xticks(x,labels);axes[1].set_ylim(0,16)
axes[1].set_ylabel('Canonical candidate loci / 10,000 annotated genes')
axes[1].set_title('Annotation-size normalization',loc='left',fontweight='bold')
for ax in axes:
    ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
fig.suptitle('Three-species NLR pilot | frozen architecture rule v1',x=.06,ha='left',fontweight='bold',fontsize=15)
fig.text(.06,.035,'Counts describe sequence architecture, not validated immune function. Categories overlap and must not be added.',fontsize=9,color='#444444')
fig.subplots_adjust(left=.07,right=.98,bottom=.20,top=.82,wspace=.34)
fig.savefig(figdir/'nlr_pilot_counts.png');fig.savefig(figdir/'nlr_pilot_counts.pdf');plt.close(fig)

colors={'NACHT':'#2457A7','NB-ARC':'#D48B46','LRR':'#469A67','Effector':'#9956A5','Other repeat':'#A0A5AC'}
def category(name):
    if name in ['NACHT','NB-ARC']:return name
    if name.startswith('LRR'):return 'LRR'
    if name in ['CARD','Death','DED']:return 'Effector'
    return 'Other repeat'
selection=[]
for species in s.species:
    candidates=canonical[(canonical.species==species)&~canonical.multiple_core_domains].copy()
    candidates['median_distance']=(candidates.length-candidates.length.median()).abs()
    selection.append(candidates.sort_values(['median_distance','protein_id']).iloc[0])
selected=pd.DataFrame(selection).drop(columns='median_distance')
selected.to_csv(out/'architecture_examples.tsv',sep='\t',index=False)
fig,axes=plt.subplots(2,1,figsize=(12,7.7),gridspec_kw={'height_ratios':[5,3.2]})
names={'Q9Y239':'Human NOD1','Q9HC29':'Human NOD2','Q96P20':'Human NLRP3','ADU79234.1':'Hydra NLR type 1','O14727':'Human APAF1'}
for ax,records,domain_table,control in [(axes[0],c,cd,True),(axes[1],selected,d,False)]:
    titles=[]
    for y,r in enumerate(records.itertuples()):
        hits=domain_table[(domain_table.protein_id==r.protein_id)&domain_table.retained]
        if not control:hits=hits[hits.genome_id==r.genome_id]
        ax.plot([1,r.length],[y,y],color='#5D6571',lw=2,zorder=1)
        for h in hits.itertuples():
            ax.add_patch(Rectangle((h.ali_from,y-.19),h.ali_to-h.ali_from+1,.38,
                                   facecolor=colors[category(h.domain_name)],edgecolor='white',lw=.35,zorder=2))
        ax.text(r.length+15,y,f'{r.length:,} aa',va='center',fontsize=8.5,color='#444444')
        title=f'{names[r.protein_id]}  |  {r.protein_id}' if control else f'{r.species}\n{r.protein_id}'
        titles.append(title)
    ax.set_yticks(np.arange(len(records)),titles,fontsize=9)
    ax.set_ylim(len(records)-.4,-.6)
    ax.set_xlim(0,max(records.length)*1.15)
    ax.set_xlabel('Amino-acid position (alignment coordinates, 1-based)')
    ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0,pad=10)
    ax.grid(axis='x',alpha=.15)
axes[0].set_title('Reference controls: three canonical, Hydra noncanonical, APAF1 NB-ARC',loc='left',fontweight='bold',fontsize=11,pad=12)
axes[1].set_title('One illustrative canonical candidate per coral species',loc='left',fontweight='bold',fontsize=11,pad=12)
fig.legend(handles=[Patch(facecolor=color,label=name) for name,color in colors.items()],loc='upper center',bbox_to_anchor=(.62,.985),ncol=5,frameon=False,fontsize=9)
fig.text(.03,.025,'Coral examples are nearest the median length among candidates without overlapping core-model flags. Only retained model hits are drawn.',fontsize=8.5,color='#444444')
fig.subplots_adjust(left=.29,right=.98,top=.90,bottom=.10,hspace=.64)
fig.savefig(figdir/'nlr_pilot_architectures.png');fig.savefig(figdir/'nlr_pilot_architectures.pdf');plt.close(fig)

lines=['NLR PILOT V1 — TECHNICAL COMPLETION REPORT','',
       'Scope: Pocillopora damicornis, Stylophora pistillata, Isopora palifera.',
       'Technical gate: PASS. This establishes reproducible classification behavior, not functional validation.',
       'Frozen v1 protocol and all 65 frozen input representations remained unchanged.',
       'No microbial diversity normalization, association test, phylogenetic model or full nine-species run was performed.',
       '', 'RESULTS',s.to_string(index=False),'','VALIDATION',
       f"All {int(s.protein_isoforms.sum()):,} protein isoforms map explicitly to {int(s.annotated_genes.sum()):,} GTF gene loci (100%).",
       'No internal-stop exclusions or unresolved loci in these three proteomes.',
       'Longest eligible isoform and any-eligible-isoform canonical counts agree in each species.',
       'All three human NLR controls satisfy the frozen architecture rule.',
       'Hydra ADU79234.1 yields DEATH-NACHT without LRR, as reported by Lange et al. Figure 4; it is a noncanonical control.',
       'APAF1 yields NB-ARC and WD40 without a canonical NACHT-LRR call. All five BLAST self-controls pass.',
       f'All {len(canonical)} canonical candidates have at least one BLAST hit to the five-sequence panel at E <= 1e-5.',
       'BLAST was supporting evidence only; identities and separate union query/subject coverage are in blast_pairs.tsv.',
       'BLAST matches may be partial and do not establish orthology, NLR function or complete domain architecture.',
       'GA, HMM-coverage, coordinate, repeat-overlap and output-hash checks passed.',
       '', 'REVIEW FLAGS',
       'Five canonical candidates (2/2/1 across the species above) also match NB-ARC over the same NACHT-region coordinates.',
       'These are overlapping model assignments, not evidence for two independent nucleotide-binding domains.',
       'They remain counted once under frozen v1; no score threshold or inclusion rule was changed after viewing counts.',
       'See overlapping_core_review.tsv for accession-level coordinates, model coverage and scores.',
       'All 85 canonical representatives have annotated start/stop codons and no ambiguous residues; this does not prove correct gene models.',
       'The panel contains no experimentally validated coral positive reference; sensitivity and specificity in coral remain unestimated.',
       'Noncanonical NACHT candidates and NB-ARC candidates are preserved separately. Categories overlap and cannot be summed.',
       'BUSCO completeness is source-reported (OCR-checked), not newly computed: 99.4%, 99.2%, 94.0%.',
       '', 'IMPLEMENTATION DETAILS',
       'HMMER 3.4: hmmsearch --cut_ga, 15 frozen models, all isoforms scanned, >=50% HMM coverage.',
       'Primary count: unique GTF locus with NACHT followed by a nonoverlapping downstream LRR in its longest eligible isoform.',
       'Representative length excludes one terminal stop. Internal-stop isoforms are flagged and excluded from counts.',
       'Repeat-model competition: greedy descending domain score, accession tie-break, >50% overlap of shorter alignment.',
       'Partial-model flags mean absent GTF start/stop feature; absence is not itself proof of biological truncation.',
       'BLAST+ 2.17.0: five references, E <= 1e-5; no identity filter and no requirement for a human match.',
       'HMMER resumes only after input, command, version and output hash checks. Raw model outputs and commands are retained.',
       '', 'REPRODUCTION',
       'conda activate bio-project',
       'python scripts/test_nlr_pilot.py',
       'python scripts/run_nlr_pilot.py',
       'python scripts/summarize_nlr_pilot.py',
       'Or execute coral_nlr_microbiome_analysis.ipynb with the Python (bio-project) kernel.',
       '', 'NEXT STAGE',
       'Apply the unchanged v1 rule to all nine matched species; retain eight-species primary and nine-species sensitivity scopes.',
       'Carry the overlapping-core flags forward for interpretation. Any new exclusion rule requires a versioned pre-association amendment.',
       'Then perform the already planned microbiome diversity preparation and species-level aggregation.',
       '', 'METHOD SOURCES',
       'HMMER user guide: https://eddylab.org/software/hmmer/CURRENT/Userguide.pdf',
       'NCBI BLAST output fields: https://www.ncbi.nlm.nih.gov/books/NBK279684/',
       'Hydra architecture: https://doi.org/10.1093/molbev/msq349',
       '']
(root/'reports/nlr_pilot_v1_report.txt').write_text('\n'.join(lines))
review_manifest={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest()
                 for f in [out/'overlapping_core_review.tsv',out/'architecture_examples.tsv',
                            figdir/'nlr_pilot_counts.png',figdir/'nlr_pilot_architectures.png',root/'reports/nlr_pilot_v1_report.txt']}
(out/'review_manifest.json').write_text(json.dumps(review_manifest,indent=2)+'\n')
print('Pilot figures, complex-core review and completion report written.')
