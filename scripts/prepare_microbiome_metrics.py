#!/usr/bin/env python3
"""Frozen-v1 alpha diversity, compositions and analysis-ready species summaries.
No association tests. The notebook calls this reproducible batch implementation.
"""
import argparse
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import pdist, squareform

METRICS=['shannon','observed_otus','inverse_simpson']
DIAGNOSTIC_DEPTHS=[100,250,500,1000,2000,5000,10000]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(frame,path,index=False):
    compression={'method':'gzip','mtime':0} if str(path).endswith('.gz') else None
    frame.to_csv(path,sep='\t',index=index,compression=compression,float_format='%.12g')

def alpha(counts):
    a=np.atleast_2d(np.asarray(counts,dtype=np.float64))
    assert np.isfinite(a).all() and (a>=0).all() and (a.sum(axis=1)>0).all()
    p=a/a.sum(axis=1,keepdims=True)
    logs=np.zeros_like(p);np.log(p,out=logs,where=p>0)
    return np.column_stack((-(p*logs).sum(axis=1),(a>0).sum(axis=1),1/(p*p).sum(axis=1)))

def rarefy(counts,depth,seed_sequence,repeats=100):
    counts=np.asarray(counts,dtype=np.int64)
    assert counts.ndim==1 and (counts>=0).all() and counts.sum()>=depth and counts.sum()<10**9
    rng=np.random.Generator(np.random.PCG64(seed_sequence))
    draws=rng.multivariate_hypergeometric(counts,depth,size=repeats,method='marginals')
    assert (draws>=0).all() and (draws<=counts).all() and (draws.sum(axis=1)==depth).all()
    return alpha(draws)

def aggregate(values,metadata):
    """Rows indexed by sample; explicit species/recorded-colony equal weighting."""
    assert values.index.is_unique and set(values.index)==set(metadata.index)
    info=metadata.loc[values.index,['canonical_species','colony_id']]
    assert info.notna().all().all() and info.ne('').all().all()
    colonies=values.groupby([info.canonical_species,info.colony_id],sort=True).mean()
    species=colonies.groupby(level=0,sort=True).mean()
    return colonies,species

def pcoa(distance):
    d=np.asarray(distance,dtype=float)
    assert np.allclose(d,d.T) and np.allclose(np.diag(d),0) and (d>=0).all()
    j=np.eye(len(d))-np.ones_like(d)/len(d)
    b=-.5*j@(d*d)@j
    eigenvalues,vectors=np.linalg.eigh(b)
    order=np.argsort(eigenvalues)[::-1];eigenvalues=eigenvalues[order];vectors=vectors[:,order]
    for k in range(len(d)):
        if vectors[np.argmax(np.abs(vectors[:,k])),k]<0:vectors[:,k]*=-1
    tolerance=max(1,float(np.max(np.abs(eigenvalues))))*1e-12
    positive=eigenvalues>tolerance;negative=eigenvalues < -tolerance
    # Retain the complete signed spectrum; plot real coordinates of positive axes only.
    assert np.allclose((vectors*eigenvalues)@vectors.T,b,atol=1e-10)
    reconstructed=np.diag(b)[:,None]+np.diag(b)[None,:]-2*b
    assert np.allclose(reconstructed,d*d,atol=1e-10)
    coords=vectors[:,positive]*np.sqrt(eigenvalues[positive])
    total=float(eigenvalues[positive].sum());neg=float(-eigenvalues[negative].sum())
    proportions=np.where(positive,eigenvalues/total if total else 0,0)
    info=dict(negative_eigenvalues=int(negative.sum()),negative_inertia=neg,positive_inertia=total,
              negative_fraction_absolute_inertia=neg/(total+neg) if total+neg else 0,
              largest_negative_to_largest_positive=float(-eigenvalues[negative].min()/eigenvalues[0]) if negative.any() else 0,
              numerical_zero_tolerance=tolerance,axis1_positive_fraction=float(proportions[0]),
              axis2_positive_fraction=float(proportions[1]) if len(proportions)>1 else 0,
              distance_correction='none',axis_percent_denominator='sum of positive eigenvalues')
    return eigenvalues,proportions,coords,info

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    root=parser.parse_args().root.resolve()
    config=json.loads((root/'config/analysis_protocol_v1.json').read_text())
    rule=config['microbiome']['alpha'];depth=int(rule['depth']);repeats=int(rule['replicates'])
    assert depth==1000 and repeats==100 and np.__version__=='2.4.6'
    frozen=list(csv.DictReader((root/'config/analysis_inputs_v1.tsv').open(),delimiter='\t'))
    for r in frozen:
        data=(root/r['path']).read_bytes()
        if r['hash_mode']=='gzip_decompressed_bytes':data=gzip.decompress(data)
        assert hashlib.sha256(data).hexdigest()==r['sha256'],r['path']
    nlr_completion=json.loads((root/'results/nlr_full_v1/completion.json').read_text())
    assert nlr_completion['status']=='complete'
    for name,digest in nlr_completion['products_sha256'].items():assert sha(root/name)==digest,name
    nlr=pd.read_csv(root/'data/processed/nlr_v1/species_exposures.tsv',sep='\t')
    out=root/'results/microbiome_metrics_v1';out.mkdir(exist_ok=True)
    dest=root/'data/processed/microbiome_metrics_v1';dest.mkdir(exist_ok=True)
    (out/'validation.json').unlink(missing_ok=True)
    # All sample and OTU IDs are lexicographically sorted, including numeric-looking OTUs.
    cohorts={}
    for name,definition in config['cohorts'].items():
        folder=root/definition['path']
        md=pd.read_csv(folder/'sample_metadata.tsv',sep='\t',dtype={'sample_id':str,'colony_id':str}).set_index('sample_id').sort_index()
        assert list(md.index)==sorted(definition['sample_ids'])
        counts=pd.read_csv(folder/'counts.tsv.gz',sep='\t',index_col=0)
        counts.index=counts.index.map(str);counts=counts.sort_index().loc[:,md.index]
        assert counts.index.is_unique and counts.columns.is_unique
        assert np.isfinite(counts.to_numpy()).all() and (counts.to_numpy()>=0).all()
        assert np.equal(counts.to_numpy(),counts.to_numpy().astype(np.int64)).all()
        assert counts.sum().ge(depth).all()
        recorded_depth=md['reads_including_unassigned' if name.endswith('with_unassigned') else 'microbial_reads'].astype(int)
        assert np.array_equal(counts.sum().to_numpy(),recorded_depth.to_numpy())
        taxonomy=pd.read_csv(folder/'taxonomy.tsv',sep='\t',dtype={'feature_id':str}).set_index('feature_id')
        assert set(counts.index)==set(taxonomy.index)
        cohorts[name]=(counts,md,taxonomy)
    sample_ids=sorted(set().union(*(set(v[1].index) for v in cohorts.values())))
    streams=np.random.SeedSequence(rule['seed']).spawn(len(sample_ids))
    seeds={}
    seed_rows=[]
    for sample_index,(sid,stream) in enumerate(zip(sample_ids,streams)):
        for d,child in zip(DIAGNOSTIC_DEPTHS,stream.spawn(len(DIAGNOSTIC_DEPTHS))):
            seeds[sid,d]=child
            seed_rows.append(dict(sample_id=sid,sample_sort_index=sample_index,depth=d,spawn_key='.'.join(map(str,child.spawn_key))))
    save(pd.DataFrame(seed_rows),out/'rng_streams.tsv')
    cache={};draw_tables={};diagnostics=[]
    def metrics_for(counts,sid,d):
        positive=counts[sid][counts[sid]>0]
        signature=hashlib.sha256(('\n'.join(positive.index)+'\n').encode()+positive.to_numpy(dtype='<i8').tobytes()).hexdigest()
        key=(sid,d,signature)
        if key not in cache:cache[key]=rarefy(positive.to_numpy(),d,seeds[sid,d],repeats)
        return cache[key]
    primary_name='primary_strict_mito_agree'
    primary_counts,primary_md,primary_tax=cohorts[primary_name]
    for sid in primary_md.index:
        observed=int((primary_counts[sid]>0).sum())
        for d in DIAGNOSTIC_DEPTHS:
            if primary_counts[sid].sum()<d:continue
            vals=metrics_for(primary_counts,sid,d)
            row=dict(sample_id=sid,species=primary_md.loc[sid,'canonical_species'],depth=d,source_reads=int(primary_counts[sid].sum()),source_observed_otus=observed)
            for k,metric in enumerate(METRICS):
                row[metric+'_mean']=float(vals[:,k].mean());row[metric+'_mc_sd']=float(vals[:,k].std(ddof=1))
            row['fraction_source_observed_otus']=float(vals[:,1].mean()/observed)
            diagnostics.append(row)
    save(pd.DataFrame(diagnostics),out/'rarefaction_diagnostics.tsv')
    variants=[dict(name=primary_name,cohort=primary_name,depth=depth,filter='none',factor='primary')]
    for name in config['sensitivities']['cohorts']:
        variants.append(dict(name=name,cohort=name,depth=depth,filter='none',factor='combined_health_and_mitotype' if name=='sensitivity_broad_all_mitotypes' else name))
    for d in config['sensitivities']['depths']:
        variants.append(dict(name=f'primary_depth_{d}',cohort=primary_name,depth=d,filter='none',factor='depth'))
    for minimum in config['sensitivities']['minimum_recorded_colonies']:
        variants.append(dict(name=f'primary_min_{minimum}_colonies',cohort=primary_name,depth=depth,filter=f'colonies:{minimum}',factor='colony_count'))
    variants.append(dict(name='primary_busco_ge_90',cohort=primary_name,depth=depth,filter='busco',factor='annotation_quality'))
    expedition_counts=primary_md.groupby('expedition').canonical_species.nunique()
    eligible_expeditions=expedition_counts[expedition_counts>=5]
    selected_expedition=sorted(eligible_expeditions.index,key=lambda e:(-eligible_expeditions[e],e))[0] if len(eligible_expeditions) else None
    if selected_expedition:
        variants.append(dict(name='primary_single_expedition',cohort=primary_name,depth=depth,filter='expedition:'+selected_expedition,factor='expedition'))
    scenario_rows=[];pcoa_infos={};membership=[]
    for variant in variants:
        name=variant['name'];counts,md,tax=cohorts[variant['cohort']]
        md=md[counts.sum().ge(variant['depth'])].copy()
        if variant['filter'].startswith('colonies:'):
            n=int(variant['filter'].split(':')[1]);numbers=md.groupby('canonical_species').colony_id.nunique()
            md=md[md.canonical_species.isin(numbers[numbers>=n].index)]
        elif variant['filter']=='busco':md=md[md.canonical_species.isin(nlr.loc[nlr.busco_complete_percent>=90,'species'])]
        elif variant['filter'].startswith('expedition:'):md=md[md.expedition==selected_expedition]
        folder=dest/name;folder.mkdir(exist_ok=True)
        save(md,folder/'sample_metadata.tsv',index=True)
        base_md=cohorts[variant['cohort']][1]
        for sid,row in base_md.iterrows():
            membership.append(dict(scenario=name,sample_id=sid,species=row.canonical_species,included=sid in md.index,
                                   reason='included' if sid in md.index else ('below_depth' if counts[sid].sum()<variant['depth'] else variant['filter'])))
        assert len(md)>0
        sample_values=[];rep_rows=[]
        for sid,row in md.iterrows():
            vals=metrics_for(counts,sid,variant['depth'])
            raw=alpha(counts[sid].to_numpy())[0]
            r=dict(sample_id=sid,species=row.canonical_species,colony_id=row.colony_id,expedition=row.expedition,site=row.site,
                   source_reads=int(counts[sid].sum()),rarefaction_depth=variant['depth'])
            for k,metric in enumerate(METRICS):
                r[metric+'_mean']=float(vals[:,k].mean());r[metric+'_mc_sd']=float(vals[:,k].std(ddof=1))
                r[metric+'_mc_mean_se']=float(vals[:,k].std(ddof=1)/np.sqrt(repeats))
                r[metric+'_mc_p05']=float(np.quantile(vals[:,k],.05));r[metric+'_mc_p95']=float(np.quantile(vals[:,k],.95))
                r[metric+'_unrarefied']=float(raw[k])
            sample_values.append(r)
            for replicate,val in enumerate(vals):rep_rows.append(dict(sample_id=sid,replicate=replicate,**dict(zip(METRICS,val))))
        sample=pd.DataFrame(sample_values).set_index('sample_id');rep=pd.DataFrame(rep_rows)
        save(sample,folder/'sample_alpha.tsv',index=True);save(rep,folder/'alpha_draw_metrics.tsv.gz')
        columns=[m+'_mean' for m in METRICS]+[m+'_unrarefied' for m in METRICS]
        colony_means,species_means=aggregate(sample[columns],md)
        save(colony_means,folder/'colony_alpha.tsv',index=True)
        species_means.index.name='species'
        species_means['n_samples']=md.groupby('canonical_species').size()
        species_means['n_recorded_colonies']=md.groupby('canonical_species').colony_id.nunique()
        species_means['n_expeditions']=md.groupby('canonical_species').expedition.nunique()
        species_means['n_sites']=md.groupby('canonical_species').site.nunique()
        species_means['single_colony_flag']=species_means.n_recorded_colonies.eq(1)
        species_means['fewer_than_3_colonies']=species_means.n_recorded_colonies.lt(3)
        species_means['colony_shannon_sd']=colony_means.groupby(level=0).shannon_mean.std(ddof=1)
        # Replicate-wise averages give Monte Carlo variability, never biological confidence intervals.
        for metric in METRICS:
            matrix=rep.pivot(index='sample_id',columns='replicate',values=metric).loc[md.index]
            _,sp_draws=aggregate(matrix,md)
            assert np.allclose(sp_draws.mean(axis=1),species_means[metric+'_mean'])
            species_means[metric+'_mc_sd']=sp_draws.std(axis=1,ddof=1)
            species_means[metric+'_mc_p05']=sp_draws.quantile(.05,axis=1)
            species_means[metric+'_mc_p95']=sp_draws.quantile(.95,axis=1)
        save(species_means,folder/'species_alpha.tsv',index=True)
        # Normalize each sample first. Do not pool reads or weight species by sample counts.
        profiles=counts.loc[:,md.index].T.astype(float)
        profiles=profiles.loc[:,profiles.sum(axis=0)>0]
        profiles=profiles.div(profiles.sum(axis=1),axis=0)
        colony_profiles,species_profiles=aggregate(profiles,md)
        assert all(np.allclose(table.sum(axis=1),1) for table in [profiles,colony_profiles,species_profiles])
        save(profiles,folder/'sample_profiles.tsv.gz',index=True)
        save(colony_profiles,folder/'colony_profiles.tsv.gz',index=True)
        species_profiles.index.name='species';save(species_profiles,folder/'species_profiles.tsv.gz',index=True)
        for level,profile in [('sample',profiles),('species',species_profiles)]:
            distance=squareform(pdist(profile.to_numpy(),metric='braycurtis'))
            assert np.isfinite(distance).all() and distance.max()<=1+1e-12
            save(pd.DataFrame(distance,index=profile.index,columns=profile.index),folder/f'{level}_bray_curtis.tsv',index=True)
            eig,prop,coords,info=pcoa(distance)
            save(pd.DataFrame(dict(axis=np.arange(1,len(eig)+1),eigenvalue=eig,positive_inertia_fraction=prop)),folder/f'{level}_pcoa_eigenvalues.tsv')
            save(pd.DataFrame(coords,index=profile.index,columns=[f'PCoA{i+1}' for i in range(coords.shape[1])]),folder/f'{level}_pcoa_coordinates.tsv',index=True)
            pcoa_infos[name+'/'+level]=info
        joined=species_means.reset_index().merge(nlr,on='species',how='left',validate='one_to_one',indicator=True)
        assert joined['_merge'].eq('both').all();joined=joined.drop(columns='_merge')
        save(joined,folder/'analysis_ready_species.tsv')
        # Source taxonomy labels are retained; this is descriptive composition aggregation.
        def family_label(fid):
            levels=str(tax.loc[fid,'taxonomy']).split(';');families=[x for x in levels if x.startswith('f__') and len(x)>3]
            if families:return families[0][3:]
            parent=next((x[3:] for prefix in ['p__','k__'] for x in levels if x.startswith(prefix) and len(x)>3),'unassigned')
            return 'Unclassified family ('+parent+')'
        family=species_profiles.T.groupby([family_label(fid) for fid in species_profiles.columns]).sum().T
        assert np.allclose(family.sum(axis=1),1)
        save(family,folder/'species_family_profiles.tsv',index=True)
        species_set=set(md.canonical_species)
        scenario_rows.append(dict(scenario=name,source_cohort=variant['cohort'],factor=variant['factor'],depth=variant['depth'],
                                  samples=len(md),recorded_colonies=len(colony_means),species=len(species_set),
                                  lost_primary_species=';'.join(sorted(set(config['primary_species'])-species_set)),
                                  gained_species_vs_primary=';'.join(sorted(species_set-set(config['primary_species']))),
                                  fewer_than_5_species=len(species_set)<5,
                                  selected_expedition=selected_expedition if variant['filter'].startswith('expedition:') else ''))
        print('PREPARED',name,len(md),'samples',len(species_set),'species',flush=True)
    save(pd.DataFrame(scenario_rows),out/'scenario_coverage.tsv')
    save(pd.DataFrame(membership),out/'scenario_membership.tsv')
    (out/'pcoa_diagnostics.json').write_text(json.dumps(pcoa_infos,indent=2)+'\n')
    settings=dict(protocol_id=config['protocol_id'],numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,
                  seed=rule['seed'],bit_generator='PCG64',sampler='multivariate_hypergeometric(method=marginals)',
                  replicates=repeats,sample_order=sample_ids,depth_order=DIAGNOSTIC_DEPTHS,
                  stream_policy='SeedSequence(seed).spawn(sorted union sample IDs), then spawn(sorted diagnostic depths). Identical sample/count vectors reuse draws across cohorts; positive OTUs sorted as strings.',
                  uncertainty='Monte Carlo resampling variability, not biological confidence intervals.',
                  aggregation='mean within species/recorded-colony, then equal mean across recorded colonies',
                  recorded_colonies='IDs retained verbatim; distinct labels do not establish biological independence',
                  expedition_selection='most eligible species, at least 5; lexicographic tie break',
                  selected_expedition=selected_expedition,unrarefied_shannon_role='sensitivity only',
                  combined_cohort_note='broad_all_mitotypes is the prelisted combined health/mitotype scenario, not a single-factor change',
                  source_script_sha256=sha(Path(__file__)))
    (out/'implementation.json').write_text(json.dumps(settings,indent=2)+'\n')
    products=[p for p in dest.rglob('*') if p.is_file()]+[out/name for name in ['rng_streams.tsv','rarefaction_diagnostics.tsv','scenario_coverage.tsv','scenario_membership.tsv','pcoa_diagnostics.json','implementation.json']]
    validation=dict(status='complete',frozen_representations_verified=len(frozen),scenario_count=len(variants),
                    unique_samples=len(sample_ids),unique_rarefaction_jobs=len(cache),draws_per_job=repeats,
                    primary_samples=25,primary_species=8,association_tests_run=False,
                    candidate_exposures_sha256=sha(root/'data/processed/nlr_v1/species_exposures.tsv'),
                    source_script_sha256=sha(Path(__file__)),
                    products_sha256={str(p.relative_to(root)):sha(p) for p in products})
    (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print('MICROBIOME PREPARATION COMPLETE',len(variants),'scenarios',flush=True)

if __name__=='__main__':main()
