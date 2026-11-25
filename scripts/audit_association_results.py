#!/usr/bin/env python3
"""Independent floating-rank/SciPy audit of the integer-rank primary tests."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import permutation_test,rankdata,pearsonr,spearmanr
from prepare_microbiome_metrics import sha

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    root=parser.parse_args().root;out=root/'results/association_v1'
    p=pd.read_csv(out/'primary_species_data.tsv',sep='\t').sort_values('species')
    results=pd.read_csv(out/'primary_results.tsv',sep='\t').set_index('endpoint')
    x=p.canonical.to_numpy();y=p.shannon_mean.to_numpy();yr=rankdata(y)
    def abs_spearman(values,axis=-1):
        ranked=rankdata(values,method='average',axis=axis)
        centered=ranked-ranked.mean(axis=axis,keepdims=True);target=yr-yr.mean()
        return np.abs(np.sum(centered*target,axis=axis)/np.sqrt(np.sum(centered**2,axis=axis)*np.sum(target**2)))
    # The frozen two-sided definition is absolute-tail, not SciPy's doubled smaller tail.
    reference_a=permutation_test((x,),abs_spearman,permutation_type='pairings',vectorized=True,
                                  n_resamples=np.inf,alternative='greater',batch=4096)
    assert np.isclose(reference_a.pvalue,results.loc['A','p_raw'],rtol=0,atol=1e-12)
    folder=root/'data/processed/microbiome_metrics_v1/primary_strict_mito_agree'
    bray=pd.read_csv(folder/'species_bray_curtis.tsv',sep='\t',index_col=0).loc[p.species,p.species].to_numpy()
    gene=np.abs(np.log1p(x[:,None])-np.log1p(x[None,:]));ii,jj=np.triu_indices(len(x),1)
    perm=np.array(list(itertools.permutations(range(len(x)))))
    # Independently permute raw distances and rerank every permutation.
    ranked=rankdata(gene[perm[:,ii],perm[:,jj]],axis=1,method='average')
    target=rankdata(bray[ii,jj]);reference_b=pearsonr(ranked,np.broadcast_to(target,ranked.shape),axis=1).statistic
    observed_b=spearmanr(gene[ii,jj],bray[ii,jj]).statistic
    extreme_b=int((np.abs(reference_b)>=abs(observed_b)-1e-12).sum())
    assert extreme_b==int(results.loc['B','extreme_count'])
    assert np.isclose(extreme_b/len(perm),results.loc['B','p_raw'],rtol=0,atol=1e-12)
    assert np.isclose(observed_b,results.loc['B','effect_rho'],rtol=0,atol=1e-12)
    all_results=pd.read_csv(out/'all_results.tsv',sep='\t')
    for name,group in all_results.groupby('scenario',sort=False):
        eligible=group[group.status=='tested']
        if not len(eligible):
            assert group.p_raw.isna().all();continue
        null=np.load(out/'permutation_nulls'/f'{name}.npz')
        for row in eligible.itertuples():
            extreme=int((np.abs(null[row.endpoint])>=abs(row.observed_rank_dot)).sum())
            assert extreme==row.extreme_count
            expected=extreme/row.permutations if row.permutation_method=='exact' else (extreme+1)/(row.permutations+1)
            assert np.isclose(row.p_raw,expected,rtol=0,atol=1e-12)
    audit=dict(status='passed',primary_A_independent_scipy_abs_statistic_p=float(reference_a.pvalue),
               primary_B_independent_raw_matrix_reranking_p=extreme_b/len(perm),
               primary_B_extreme_count=extreme_b,primary_permutations=len(perm),
               all_saved_nulls_reproduce_reported_pvalues=True,undefined_pvalues_not_imputed=True,
               primary_results_sha256=sha(out/'primary_results.tsv'),source_script_sha256=sha(Path(__file__)))
    (out/'independent_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))

if __name__=='__main__':main()
