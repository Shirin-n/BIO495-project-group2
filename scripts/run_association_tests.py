#!/usr/bin/env python3
"""Frozen-v1 species-level Spearman and Spearman Mantel permutation analyses.

Two-sided p values count absolute statistics (not doubled smaller tails).
Integer-scaled average ranks allow exact extreme-count comparisons without
floating-point tolerance decisions. Sensitivity p values are unadjusted.
"""
import argparse
import csv
import gzip
import itertools
import json
import math
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata,spearmanr
from prepare_microbiome_metrics import sha,save

PRIMARY='primary_strict_mito_agree'

def integer_centered_ranks(x):
    x=np.asarray(x,dtype=float)
    assert x.ndim==1 and np.isfinite(x).all()
    doubled=np.rint(2*rankdata(x,method='average')).astype(np.int64)
    return len(x)*doubled-doubled.sum()

def rank_effect(x,y):
    a=integer_centered_ranks(x);b=integer_centered_ranks(y)
    denom=np.sqrt(float(a@a)*float(b@b))
    return (float(a@b)/denom if denom else np.nan),a,b,denom

@lru_cache(maxsize=None)
def exact_permutations(n):
    assert n<=8
    return np.array(list(itertools.permutations(range(n))),dtype=np.int16)

def permutations(n,seed,random_replicates=9999):
    if n<=8:return exact_permutations(n),'exact'
    rng=np.random.Generator(np.random.PCG64(seed))
    return np.array([rng.permutation(n) for _ in range(random_replicates)],dtype=np.int16),'monte_carlo'

def test_endpoints(x,y,bray,perms=None,mode='descriptive',minimum_species=5):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float);bray=np.asarray(bray,dtype=float)
    n=len(x);assert len(y)==n and bray.shape==(n,n)
    assert np.isfinite(x).all() and (x>=0).all() and np.isfinite(bray).all()
    assert np.allclose(bray,bray.T) and np.allclose(np.diag(bray),0)
    triu=np.triu_indices(n,1);gene=np.abs(np.log1p(x[:,None])-np.log1p(x[None,:]))
    a_rho,a,b,a_denom=rank_effect(x,y)
    b_rho,g,m,b_denom=rank_effect(gene[triu],bray[triu])
    assert not np.isfinite(a_rho) or np.isclose(a_rho,spearmanr(x,y).statistic)
    assert not np.isfinite(b_rho) or np.isclose(b_rho,spearmanr(gene[triu],bray[triu]).statistic)
    gene_rank=np.zeros((n,n),dtype=np.int64);gene_rank[triu]=g;gene_rank[(triu[1],triu[0])]=g
    results=[];nulls={}
    for endpoint,effect,left,right,denom in [('A',a_rho,a,b,a_denom),('B',b_rho,g,m,b_denom)]:
        numerator=int(left@right)
        row=dict(endpoint=endpoint,n_species=n,n_pairs=n*(n-1)//2 if endpoint=='B' else None,
                 effect_rho=effect,status='tested',reason='',permutation_method=mode,
                 permutations=0,extreme_count=None,p_raw=np.nan,p_holm_primary=np.nan,
                 random_p_mc_se=np.nan,observed_rank_dot=numerator,rank_norm_product=denom)
        if not denom:
            row.update(status='undefined',reason='constant exposure, endpoint or distance vector',permutation_method='none')
        elif n<minimum_species:
            row.update(status='descriptive_only',reason='fewer than five species',permutation_method='none')
        elif perms is None:
            row.update(status='effect_only',reason='leave-one-species-out effect diagnostic',permutation_method='none')
        else:
            if endpoint=='A':dots=(a[perms]@b).astype(np.int64)
            else:dots=(gene_rank[perms[:,triu[0]],perms[:,triu[1]]]@m).astype(np.int64)
            extreme=int((np.abs(dots)>=abs(numerator)).sum())
            p=extreme/len(perms) if mode=='exact' else (extreme+1)/(len(perms)+1)
            if mode=='exact':assert np.array_equal(perms[0],np.arange(n)) and dots[0]==numerator and extreme>=1
            row.update(permutations=len(perms),extreme_count=extreme,p_raw=p,
                       random_p_mc_se=np.sqrt(p*(1-p)/(len(perms)+1)) if mode=='monte_carlo' else np.nan)
            nulls[endpoint]=dots
        results.append(row)
    return results,nulls,gene

def holm_two(pvalues):
    assert len(pvalues)==2
    p=np.asarray(pvalues,dtype=float);defined=np.isfinite(p)
    filled=np.where(defined,p,1.0);order=np.argsort(filled,kind='stable')
    adjusted=np.maximum.accumulate(filled[order]*np.array([2,1]))
    result=np.empty(2);result[order]=np.minimum(adjusted,1)
    result[~defined]=np.nan
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    root=parser.parse_args().root.resolve()
    protocol=json.loads((root/'config/analysis_protocol_v1.json').read_text())
    assert protocol['protocol_id']=='coral-nlr-v1.0'
    assert protocol['statistics']['exact_permutations_up_to_species']==8
    assert np.__version__=='2.4.6'
    frozen=list(csv.DictReader((root/'config/analysis_inputs_v1.tsv').open(),delimiter='\t'))
    import hashlib
    for r in frozen:
        data=(root/r['path']).read_bytes()
        if r['hash_mode']=='gzip_decompressed_bytes':data=gzip.decompress(data)
        assert hashlib.sha256(data).hexdigest()==r['sha256'],r['path']
    inputs=root/'data/processed/microbiome_metrics_v1'
    prep=root/'results/microbiome_metrics_v1'
    prep_validation=json.loads((prep/'validation.json').read_text())
    assert prep_validation['status']=='complete'
    for name,digest in prep_validation['products_sha256'].items():assert sha(root/name)==digest,name
    assert prep_validation['candidate_exposures_sha256']==sha(root/'data/processed/nlr_v1/species_exposures.tsv')
    coverage=pd.read_csv(prep/'scenario_coverage.tsv',sep='\t').fillna('')
    out=root/'results/association_v1';out.mkdir(exist_ok=True)
    null_dir=out/'permutation_nulls';null_dir.mkdir(exist_ok=True)
    (out/'validation.json').unlink(missing_ok=True)
    scenarios=[dict(scenario=r.scenario,data_scenario=r.scenario,exposure='canonical',alpha='shannon_mean',factor=r.factor)
               for r in coverage.itertuples()]
    for exposure in ['any_isoform_canonical','nacht','nb_arc','canonical_per_10000_genes']:
        scenarios.append(dict(scenario='exposure_'+exposure,data_scenario=PRIMARY,exposure=exposure,alpha='shannon_mean',factor='exposure'))
    scenarios.append(dict(scenario='unrarefied_shannon',data_scenario=PRIMARY,exposure='canonical',alpha='shannon_unrarefied',factor='alpha_normalization'))
    scenario_names=sorted(r['scenario'] for r in scenarios)
    seeds=dict(zip(scenario_names,np.random.SeedSequence(protocol['statistics']['seed']).spawn(len(scenario_names))))
    specification=dict(protocol_id=protocol['protocol_id'],scenarios=scenarios,seed=protocol['statistics']['seed'],
        random_streams={name:list(seeds[name].spawn_key) for name in scenario_names},
        rng='PCG64; SeedSequence children in sorted scenario-name order; random label permutations with replacement',
        two_sided='count abs(permuted rank dot product) >= abs(observed rank dot product)',
        integer_rank_scaling='z = n * (2 * average_rank) - sum(2 * average_rank); constant positive scaling cancels in rho',
        ties='average ranks; all labelled permutations retained, including duplicate statistics',
        holm='two primary endpoints only; undefined endpoint reserves one of the two hypotheses',
        sensitivity_pvalues='unadjusted exploratory; never a selection criterion for primary conclusions',
        density_exposure='10000 * canonical / annotated_genes, recomputed without reporting-rounding',
        matrix_permutation='same label permutation on rows and columns of gene-distance rank matrix',
        leave_one_species_out='effects only; ranges are diagnostics, not confidence intervals',
        numpy=np.__version__,scipy=scipy.__version__,source_script_sha256=sha(Path(__file__)))
    # Save execution details before evaluating any association statistic.
    (out/'execution_specification.json').write_text(json.dumps(specification,indent=2)+'\n')
    rows=[];species_rows=[];null_manifest=[];primary_data=None
    for scenario in scenarios:
        name=scenario['scenario'];folder=inputs/scenario['data_scenario']
        table=pd.read_csv(folder/'analysis_ready_species.tsv',sep='\t').sort_values('species').reset_index(drop=True)
        ids=list(table.species)
        dist=pd.read_csv(folder/'species_bray_curtis.tsv',sep='\t',index_col=0).loc[ids,ids].to_numpy()
        x=(10000*table.canonical/table.annotated_genes).to_numpy() if scenario['exposure']=='canonical_per_10000_genes' else table[scenario['exposure']].to_numpy()
        y=table[scenario['alpha']].to_numpy()
        perm,mode=permutations(len(ids),seeds[name],protocol['statistics']['monte_carlo_permutations_above_8_species']) if len(ids)>=5 else (None,'none')
        result,nulls,gene=test_endpoints(x,y,dist,perm,mode)
        metadata=dict(scenario=name,data_scenario=scenario['data_scenario'],factor=scenario['factor'],
                      exposure=scenario['exposure'],alpha_metric=scenario['alpha'],n_samples=int(table.n_samples.sum()),
                      n_recorded_colonies=int(table.n_recorded_colonies.sum()),
                      excluded_primary_species=';'.join(sorted(set(protocol['primary_species'])-set(ids))),
                      additional_species=';'.join(sorted(set(ids)-set(protocol['primary_species']))),
                      pvalue_family='primary_two_test_Holm' if name==PRIMARY else 'unadjusted_sensitivity',
                      note='B unchanged from primary; alpha endpoint only changes' if name=='unrarefied_shannon' else '')
        if name==PRIMARY:
            assert set(ids)==set(protocol['primary_species']) and len(ids)==8
            adjusted=holm_two([r['p_raw'] for r in result])
            for r,p in zip(result,adjusted):r['p_holm_primary']=float(p)
            primary_data=(table,x,y,dist,gene)
        for r in result:rows.append({**metadata,**r})
        for i,species in enumerate(ids):
            species_rows.append(dict(scenario=name,species=species,exposure_value=float(x[i]),alpha_value=float(y[i]),
                                     n_samples=int(table.n_samples.iloc[i]),n_recorded_colonies=int(table.n_recorded_colonies.iloc[i])))
        if nulls:
            path=null_dir/(name+'.npz');np.savez_compressed(path,**nulls)
            null_manifest.append(dict(scenario=name,path=str(path.relative_to(root)),sha256=sha(path),
                                      permutations=len(perm),method=mode,permutation_order_sha256=hashlib.sha256(perm.tobytes()).hexdigest()))
        print('TESTED',name,[(r['endpoint'],r['effect_rho'],r['p_raw'],r['status']) for r in result],flush=True)
    result_frame=pd.DataFrame(rows)
    save(result_frame,out/'all_results.tsv')
    save(result_frame[result_frame.scenario.eq(PRIMARY)],out/'primary_results.tsv')
    save(result_frame[~result_frame.scenario.eq(PRIMARY)],out/'sensitivity_results.tsv')
    save(pd.DataFrame(species_rows),out/'scenario_species_values.tsv')
    save(pd.DataFrame(null_manifest),out/'null_manifest.tsv')
    table,x,y,dist,gene=primary_data
    save(table,out/'primary_species_data.tsv')
    pair_rows=[]
    for i,j in zip(*np.triu_indices(len(table),1)):
        pair_rows.append(dict(species_i=table.species.iloc[i],species_j=table.species.iloc[j],gene_log1p_distance=gene[i,j],bray_curtis=dist[i,j]))
    save(pd.DataFrame(pair_rows),out/'primary_pairwise_data.tsv')
    leave_one=[]
    for i,species in enumerate(table.species):
        keep=np.arange(len(table))!=i
        results,_,_=test_endpoints(x[keep],y[keep],dist[np.ix_(keep,keep)])
        for row in results:leave_one.append(dict(omitted_species=species,**row))
    save(pd.DataFrame(leave_one),out/'leave_one_species_out.tsv')
    md=pd.read_csv(inputs/PRIMARY/'sample_metadata.tsv',sep='\t')
    save(pd.crosstab(md.canonical_species,md.site),out/'species_site_sample_counts.tsv',index=True)
    save(pd.crosstab(md.canonical_species,md.expedition),out/'species_expedition_sample_counts.tsv',index=True)
    tree_status=dict(status='not_performed',reason='Frozen source tree not accepted as species-level PGLS input; no validated replacement present.',
                     frozen_tree=str(root/'data/raw/phylogeny/pollock2018/host_tree_from_step_11.newick'),
                     mapping_sha256=sha(root/'config/species_tree_map.tsv'),ancestry_adjusted=False,
                     lambda_1_run=False,lambda_0_contrast_run=False,
                     geographic_adjustment='none; preplanned single-expedition sensitivity reported separately')
    (out/'phylogeny_status.json').write_text(json.dumps(tree_status,indent=2)+'\n')
    diagnostics=dict(primary_results=result_frame[result_frame.scenario.eq(PRIMARY)].replace({np.nan:None}).to_dict('records'),
                     leave_one_species_out_ranges={endpoint:dict(min=float(pd.DataFrame(leave_one).query('endpoint == @endpoint').effect_rho.min()),
                         max=float(pd.DataFrame(leave_one).query('endpoint == @endpoint').effect_rho.max())) for endpoint in ['A','B']},
                     pgls_status=tree_status,scenarios=len(scenarios),sensitivity_tests=int((result_frame.scenario.ne(PRIMARY)&result_frame.status.eq('tested')).sum()),
                     sensitivity_nominal_p_below_0_05=int((result_frame.scenario.ne(PRIMARY)&result_frame.p_raw.lt(.05)).sum()),
                     interpretation='Exploratory nominal p-values; species exchangeability, ancestry, geography and annotation uncertainty limit inference.')
    (out/'summary_metrics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    products=[out/name for name in ['all_results.tsv','primary_results.tsv','sensitivity_results.tsv','scenario_species_values.tsv','null_manifest.tsv',
               'primary_species_data.tsv','primary_pairwise_data.tsv','leave_one_species_out.tsv','species_site_sample_counts.tsv',
               'species_expedition_sample_counts.tsv','phylogeny_status.json','summary_metrics.json','execution_specification.json']]
    products += [root/r['path'] for r in null_manifest]
    validation=dict(status='complete',protocol_id=protocol['protocol_id'],frozen_input_representations_verified=len(frozen),
                    scenario_count=len(scenarios),result_rows=len(rows),primary_permutations=math.factorial(8),
                    source_script_sha256=sha(Path(__file__)),microbiome_preparation_validation_sha256=sha(prep/'validation.json'),
                    products_sha256={str(p.relative_to(root)):sha(p) for p in products})
    (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print('ASSOCIATION ANALYSIS COMPLETE',flush=True)

if __name__=='__main__':main()
