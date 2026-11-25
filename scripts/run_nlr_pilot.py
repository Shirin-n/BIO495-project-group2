#!/usr/bin/env python3
"""Run the frozen v1 three-species technical pilot, without microbiome outcomes.

Usage: python scripts/run_nlr_pilot.py --root /path/to/bio-project
Raw data and frozen protocol files are read-only. Derived outputs are regenerated.
HMMER jobs resume only when their command, inputs and tool version match.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def fasta(path):
    result = {}
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        identifier, sequence = None, []
        for line in handle:
            if line.startswith('>'):
                if identifier is not None:
                    assert identifier not in result, identifier
                    result[identifier] = ''.join(sequence).upper()
                identifier, sequence = line[1:].split()[0], []
            else:
                sequence.append(line.strip())
        if identifier is not None:
            assert identifier not in result, identifier
            result[identifier] = ''.join(sequence).upper()
    return result


def write_fasta(path, records):
    with path.open('w') as handle:
        for identifier, sequence in records.items():
            handle.write(f'>{identifier}\n')
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i:i+80] + '\n')


def parse_gtf(path):
    mapping, genes, transcript_rows = {}, set(), set()
    features = defaultdict(set)
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            assert len(fields) == 9, 'Malformed GTF record'
            kind, attributes = fields[2], fields[8]
            attrs = dict(re.findall(r'(\w+)\s+"([^"]+)"', attributes))
            if kind == 'gene':
                genes.add(attrs.get('gene_id', attributes.strip().rstrip(';')))
            if kind == 'transcript':
                transcript_rows.add(attrs.get('transcript_id', attributes.strip().rstrip(';')))
            # Mapping is explicitly supported by GTF gene_id/transcript_id pairs.
            if 'transcript_id' in attrs and 'gene_id' in attrs:
                tid, gid = attrs['transcript_id'], attrs['gene_id']
                assert tid not in mapping or mapping[tid] == gid, f'Conflicting mapping: {tid}'
                mapping[tid] = gid
                features[tid].add(kind)
    assert set(mapping.values()) <= genes, 'Transcript refers to missing GTF gene'
    assert set(mapping) == transcript_rows, 'GTF transcript rows and explicit mappings differ'
    return mapping, genes, features


def overlap_fraction(a, b):
    overlap = max(0, min(a['ali_to'], b['ali_to']) - max(a['ali_from'], b['ali_from']) + 1)
    return overlap / min(a['ali_to']-a['ali_from']+1, b['ali_to']-b['ali_from']+1)


def parse_domains(path, models, sequences, settings):
    rows, by_protein = [], defaultdict(list)
    repeats = set(settings['lrr_models']) | {'PF00023.38', 'PF00400.40', 'PF00515.36'}
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        f = line.split(maxsplit=22)
        accession = f[4]
        m = models[accession]
        sequence_ga, domain_ga = map(float, m['gathering_threshold'].strip(';').split())
        r = dict(protein_id=f[0], model=accession, domain_name=f[3], protein_length=int(f[2]),
                 model_length=int(f[5]), sequence_evalue=float(f[6]), sequence_score=float(f[7]),
                 domain_index=int(f[9]), reported_domain_count=int(f[10]),
                 conditional_evalue=float(f[11]), independent_evalue=float(f[12]),
                 domain_score=float(f[13]), hmm_from=int(f[15]), hmm_to=int(f[16]),
                 ali_from=int(f[17]), ali_to=int(f[18]), env_from=int(f[19]), env_to=int(f[20]),
                 sequence_ga=sequence_ga, domain_ga=domain_ga)
        assert r['model_length'] == int(m['length'])
        assert r['protein_length'] == len(sequences[r['protein_id']])
        assert 1 <= r['hmm_from'] <= r['hmm_to'] <= r['model_length']
        assert 1 <= r['env_from'] <= r['ali_from'] <= r['ali_to'] <= r['env_to'] <= r['protein_length']
        r['model_coverage'] = (r['hmm_to']-r['hmm_from']+1) / r['model_length']
        r['passes_ga'] = r['sequence_score'] >= sequence_ga and r['domain_score'] >= domain_ga
        r['passes_coverage'] = r['model_coverage'] >= settings['minimum_model_coverage']
        r['retained'] = r['passes_ga'] and r['passes_coverage']
        r['rejection_reason'] = '' if r['retained'] else ('below_GA' if not r['passes_ga'] else 'model_coverage_below_0.5')
        by_protein[r['protein_id']].append(r)
        rows.append(r)
    for hits in by_protein.values():
        competing = sorted([h for h in hits if h['retained'] and h['model'] in repeats],
                           key=lambda h: (-h['domain_score'], h['model'], h['ali_from'], h['ali_to']))
        kept = []
        for h in competing:
            if any(overlap_fraction(h, k) > 0.5 for k in kept):
                h['retained'], h['rejection_reason'] = False, 'competing_repeat_overlap_gt_0.5'
            else:
                kept.append(h)
    return rows, by_protein


def architecture(hits, settings):
    kept = [h for h in hits if h['retained']]
    core = [h for h in kept if h['model'] == settings['core_model']]
    nb = [h for h in kept if h['model'] == settings['broad_core_model']]
    lrr = [h for h in kept if h['model'] in settings['lrr_models']]
    return dict(
                    canonical=any(l['ali_from'] > n['ali_to'] for n in core for l in lrr),
                    nacht=bool(core), nb_arc=bool(nb), nacht_domains=len(core), nb_arc_domains=len(nb),
                    lrr_domains=len(lrr), multiple_core_domains=len(core)+len(nb) > 1,
                    mixed_core_models=bool(core and nb),
                    overlapping_core_models=any(overlap_fraction(n,b)>0.5 for n in core for b in nb),
                    domains=';'.join(f"{h['domain_name']}:{h['ali_from']}-{h['ali_to']}" for h in sorted(kept, key=lambda x:x['ali_from'])))


def union_length(intervals):
    end, total = 0, 0
    for left, right in sorted(intervals):
        total += max(0, right-max(left, end+1)+1)
        end = max(end, right)
    return total


@contextmanager
def blast_database_path(database):
    # BLAST interprets whitespace in database paths as a separator.
    if not any(char.isspace() for char in str(database)):
        yield database
        return
    with tempfile.TemporaryDirectory(prefix='bio495-blast-', dir='/tmp') as directory:
        alias = Path(directory) / 'interim'
        alias.symlink_to(database.parent)
        yield alias / database.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    root = args.root.resolve()
    config = json.loads((root/'config/analysis_protocol_v1.json').read_text())
    settings = config['nlr']
    assert config['protocol_id'] == 'coral-nlr-v1.0'
    # Check every frozen representation before creating any results.
    for row in read_tsv(root/'config/analysis_inputs_v1.tsv'):
        data = (root/row['path']).read_bytes()
        if row['hash_mode'] == 'gzip_decompressed_bytes':
            data = gzip.decompress(data)
        assert hashlib.sha256(data).hexdigest() == row['sha256'], row['path']
    out = root/'results/nlr_pilot_v1'
    interim = root/'data/interim/nlr_pilot_v1'
    out.mkdir(parents=True, exist_ok=True)
    interim.mkdir(parents=True, exist_ok=True)
    # An interrupted/failed rerun cannot leave a stale success marker.
    (out/'validation.json').unlink(missing_ok=True)
    tools = {name: str(Path(sys.executable).parent/name) for name in ['hmmsearch','blastp','makeblastdb']}
    versions = {name: subprocess.check_output([exe, '-h' if name == 'hmmsearch' else '-version'], text=True).splitlines()[:3] for name,exe in tools.items()}
    commands = []
    def run(command, logfile):
        print('RUN', ' '.join(map(str, command)), flush=True)
        start = time.monotonic()
        with logfile.open('w') as handle:
            subprocess.run(list(map(str,command)), stdout=handle, stderr=subprocess.STDOUT, check=True)
        commands.append(dict(command=list(map(str,command)), seconds=round(time.monotonic()-start,2), log=str(logfile.relative_to(root))))
        (out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    models = {r['accession']:r for r in read_tsv(root/'config/hmm_model_manifest.tsv')}
    assert len(models) == 15
    panel = interim/'architecture_models.hmm'
    panel.write_bytes(b'\n'.join(gzip.decompress((root/m['path']).read_bytes()) for m in models.values()))
    def scan(label, sequences):
        faa, dom, stamp = (interim/f'{label}{suffix}' for suffix in ['.faa','.domtblout','.scan.json'])
        write_fasta(faa, {k:v.replace('*','X') for k,v in sequences.items()})
        command = [tools['hmmsearch'],'--cut_ga','--cpu',str(args.threads),'--seed','42','--noali','--domtblout',str(dom),str(panel),str(faa)]
        fingerprint = dict(fasta_sha256=sha(faa),panel_sha256=sha(panel),command=command,version=versions['hmmsearch'])
        prior = json.loads(stamp.read_text()) if stamp.exists() else {}
        if prior.get('fingerprint') == fingerprint and dom.exists() and prior.get('domtblout_sha256') == sha(dom):
            print('REUSE verified HMMER output:', label, flush=True)
            commands.append(dict(command=command,cached=True,domtblout_sha256=sha(dom)))
        else:
            run(command,interim/f'{label}.hmmsearch.log')
            assert dom.read_text().rstrip().endswith('# [ok]'), 'HMMER did not complete'
            stamp.write_text(json.dumps(dict(fingerprint=fingerprint,domtblout_sha256=sha(dom)),indent=2)+'\n')
        return parse_domains(dom,models,sequences,settings)

    controls = fasta(root/'data/reference/nlr_sequences/nlr_architecture_references.faa')
    controls.update(fasta(root/'data/reference/nlr_sequences/non_nlr_specificity_comparator.faa'))
    control_hits, control_by = scan('controls',controls)
    control_rows = [dict(protein_id=k,length=len(v),**architecture(control_by[k],settings)) for k,v in controls.items()]
    write_tsv(out/'control_domains.tsv',control_hits)
    write_tsv(out/'control_classification.tsv',control_rows)
    control_index = {r['protein_id']:r for r in control_rows}
    positive_ok = all(control_index[k]['canonical'] for k in ['Q9Y239','Q9HC29','Q96P20'])
    specificity_ok = control_index['O14727']['nb_arc'] and not control_index['O14727']['canonical']
    hydra_ok = (control_index['ADU79234.1']['nacht'] and not control_index['ADU79234.1']['canonical']
                and control_index['ADU79234.1']['lrr_domains'] == 0
                and 'Death:' in control_index['ADU79234.1']['domains'])
    assert positive_ok, 'A human positive control failed the frozen architecture rule; inspect before proceeding'
    assert specificity_ok, 'APAF1 specificity comparator failed; inspect before proceeding'
    assert hydra_ok, 'Hydra DEATH-NACHT control differs from the published architecture; inspect before proceeding'
    print('CONTROL CHECK:',json.dumps(control_rows),flush=True)

    qc_index = {r['species']:r for r in read_tsv(root/'results/qc/genome_annotation_quality.tsv')}
    summaries, proteins_all, loci_all, domains_all, candidate_sequences = [], [], [], [], {}
    mapping_qc = []
    for species in settings['pilot_species']:
        quality = qc_index[species]
        gid = quality['genome_id']
        folder = root/'data/raw/genomes/metz_v1'/f"{gid}_{species.replace(' ','_')}"
        seqs_raw = fasta(folder/'braker.aa.gz')
        seqs = {k:v[:-1] if v.endswith('*') else v for k,v in seqs_raw.items()}
        mapping, genes, features = parse_gtf(folder/'braker.gtf.gz')
        assert set(seqs) == set(mapping), f'GTF/FASTA transcript mismatch: {species}'
        assert len(genes) == int(quality['source_gene_count']), f'Gene count mismatch: {species}'
        assert set(mapping.values()) == genes, 'GTF genes without protein records'
        assert all(seqs.values()), 'Empty protein'
        assert all(set(s) <= set('ACDEFGHIKLMNPQRSTVWYBJOUXZ*') for s in seqs.values()), 'Invalid protein alphabet'
        by_gene = defaultdict(list)
        for protein_id,gene in mapping.items(): by_gene[gene].append(protein_id)
        representatives = {}
        for gene, identifiers in by_gene.items():
            eligible = [p for p in identifiers if '*' not in seqs[p]]
            representatives[gene] = sorted(eligible,key=lambda p:(-len(seqs[p]),p))[0] if eligible else None
        write_fasta(interim/f'{gid}.representatives.faa', {p:seqs[p] for p in representatives.values() if p is not None})
        hits, by_protein = scan(gid,seqs)
        for h in hits: domains_all.append(dict(species=species,genome_id=gid,gene_id=mapping[h['protein_id']],**h))
        pr = {}
        for pid,seq in seqs.items():
            row = dict(species=species,genome_id=gid,gene_id=mapping[pid],protein_id=pid,length=len(seq),
                       terminal_stop_removed=seqs_raw[pid].endswith('*'),internal_stop_count=seq.count('*'),
                       ambiguous_residues=sum(c not in 'ACDEFGHIKLMNPQRSTVWY*' for c in seq),
                       missing_gtf_start_codon='start_codon' not in features[pid],
                       missing_gtf_stop_codon='stop_codon' not in features[pid],
                       is_representative=representatives[mapping[pid]] == pid,eligible_for_count='*' not in seq,
                       **architecture(by_protein[pid],settings))
            pr[pid] = row
            proteins_all.append(row)
            if row['is_representative'] and (row['nacht'] or row['nb_arc']):
                candidate_sequences[f'{gid}__{pid}'] = seq
        lr = []
        for gene, ids in by_gene.items():
            rep = representatives[gene]
            eligible = [pr[p] for p in ids if pr[p]['eligible_for_count']]
            row = dict(species=species,genome_id=gid,gene_id=gene,isoform_count=len(ids),
                       representative_protein=rep or '',status='resolved' if rep else 'unresolved_all_internal_stops',
                       canonical=pr[rep]['canonical'] if rep else '',
                       nacht=pr[rep]['nacht'] if rep else '',nb_arc=pr[rep]['nb_arc'] if rep else '',
                       any_isoform_canonical=any(p['canonical'] for p in eligible) if eligible else '',
                       any_isoform_nacht=any(p['nacht'] for p in eligible) if eligible else '',
                       any_isoform_nb_arc=any(p['nb_arc'] for p in eligible) if eligible else '')
            lr.append(row)
        loci_all.extend(lr)
        unresolved = sum(r['status'] != 'resolved' for r in lr)
        sums = {key:sum(r[key] is True for r in lr) for key in ['canonical','nacht','nb_arc','any_isoform_canonical','any_isoform_nacht','any_isoform_nb_arc']}
        summaries.append(dict(species=species,genome_id=gid,annotated_genes=len(genes),protein_isoforms=len(seqs),
                              mapped_proteins=len(mapping),mapping_percent=100.0,resolved_genes=len(genes)-unresolved,
                              unresolved_genes=unresolved,internal_stop_proteins=sum('*' in s for s in seqs.values()),
                              **sums,canonical_per_10000_genes=round(sums['canonical']/len(genes)*10000,4),
                              representative_core_candidates=sums['nacht']+sums['nb_arc']-sum(r['nacht'] is True and r['nb_arc'] is True for r in lr),
                              busco_complete_percent=float(quality['busco_C']),
                              canonical_partial_gtf_flags=sum(r['is_representative'] and r['canonical'] and (r['missing_gtf_start_codon'] or r['missing_gtf_stop_codon']) for r in pr.values()),
                              canonical_ambiguous_residue_flags=sum(r['is_representative'] and r['canonical'] and r['ambiguous_residues']>0 for r in pr.values()),
                              canonical_multiple_core_flags=sum(r['is_representative'] and r['canonical'] and r['multiple_core_domains'] for r in pr.values())))
        mapping_qc.append(dict(species=species,protein_records=len(seqs),gtf_transcripts=len(mapping),gtf_genes=len(genes),
                               source_gene_count=int(quality['source_gene_count']),unmapped_proteins=0,missing_protein_transcripts=0,
                               genes_without_proteins=0,representatives=len(genes)-unresolved,
                               note='Mapping verified using explicit GTF gene_id/transcript_id pairs; no suffix inference.'))
        print('SPECIES SUMMARY:',json.dumps(summaries[-1]),flush=True)
    write_tsv(out/'species_summary.tsv',summaries)
    write_tsv(out/'mapping_qc.tsv',mapping_qc)
    write_tsv(out/'protein_classification.tsv',proteins_all)
    write_tsv(out/'locus_classification.tsv',loci_all)
    write_tsv(out/'domain_evidence.tsv',domains_all)
    write_tsv(out/'canonical_candidates.tsv',[p for p in proteins_all if p['is_representative'] and p['canonical']],fields=list(proteins_all[0]))

    # BLAST supports sequence review only. Counts above never depend on BLAST hits.
    ref = interim/'blast_reference_panel.faa'
    write_fasta(ref,controls)
    db = interim/'reference_panel_db'
    queries = {**{f'control__{p}':s for p,s in controls.items()},**candidate_sequences}
    query = interim/'blast_queries.faa'
    write_fasta(query,queries)
    blastout = interim/'reference_hits.tsv'
    fields = 'qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp'.split()
    with blast_database_path(db) as blast_db:
        run([tools['makeblastdb'],'-in',blast_db.parent/ref.name,'-dbtype','prot','-parse_seqids','-out',blast_db],interim/'makeblastdb.log')
        run([tools['blastp'],'-query',query,'-db',blast_db,'-evalue',str(settings['blast_reporting']['evalue_max']),
             '-num_threads',args.threads,'-max_target_seqs',len(controls),'-outfmt','6 '+' '.join(fields),'-out',blastout],interim/'blastp.log')
    with blastout.open() as handle:
        hsps = list(csv.DictReader(handle,fieldnames=fields,delimiter='\t'))
    grouped = defaultdict(list)
    for h in hsps:
        # makeblastdb recognizes accessions and may emit sp|ACC| / gb|ACC|.
        # Resolve exactly one known source ID; preserve original output in interim.
        matching = {s for s in h['sseqid'].split('|') if s in controls}
        assert len(matching) == 1, f"Unknown or ambiguous BLAST subject: {h['sseqid']}"
        h['sseqid'] = matching.pop()
        grouped[h['qseqid'],h['sseqid']].append(h)
    pairs = []
    for (qid,sid),hh in grouped.items():
        assert qid in queries and sid in controls
        qlen,slen = len(queries[qid]),len(controls[sid])
        for h in hh:
            assert int(h['qlen'])==qlen and int(h['slen'])==slen
            assert 1 <= int(h['qstart']) <= int(h['qend']) <= qlen
            assert 1 <= int(h['sstart']) <= int(h['send']) <= slen
            assert float(h['evalue']) <= settings['blast_reporting']['evalue_max']
        best = max(hh,key=lambda h:float(h['bitscore']))
        pairs.append(dict(query_id=qid,reference_id=sid,hsp_count=len(hh),best_hsp_pident=float(best['pident']),
                          best_bitscore=float(best['bitscore']),minimum_evalue=min(float(h['evalue']) for h in hh),
                          query_coverage_union=union_length([(int(h['qstart']),int(h['qend'])) for h in hh])/qlen,
                          subject_coverage_union=union_length([(int(h['sstart']),int(h['send'])) for h in hh])/slen))
    write_tsv(out/'blast_hsps.tsv',hsps,fields=fields)
    write_tsv(out/'blast_pairs.tsv',pairs)
    assert all(any(p['query_id']=='control__'+c and p['reference_id']==c and p['query_coverage_union']>.99 for p in pairs) for c in controls), 'BLAST self-control failed'
    hit_queries = {p['query_id'] for p in pairs}
    for r in summaries:
        prefix = r['genome_id']+'__'
        r['core_candidates_with_reference_hit'] = sum(q.startswith(prefix) and q in hit_queries for q in candidate_sequences)
        canonical_ids = {prefix+p['protein_id'] for p in proteins_all if p['species']==r['species'] and p['is_representative'] and p['canonical']}
        r['canonical_with_reference_hit'] = len(canonical_ids & hit_queries)
    write_tsv(out/'species_summary.tsv',summaries)
    (out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    validation = dict(protocol_id=config['protocol_id'],stage='three-species NLR technical pilot',
                      completed_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),tools=versions,
                      frozen_input_representations_verified=len(read_tsv(root/'config/analysis_inputs_v1.tsv')),
                      human_canonical_controls_pass=positive_ok,apaf1_specificity_pass=specificity_ok,
                      hydra_noncanonical_control_pass=hydra_ok,
                      hydra_expected_architecture_source='https://doi.org/10.1093/molbev/msq349 (Figure 4; DEATH-NACHT without LRRs)',
                      hydra_canonical=control_index['ADU79234.1']['canonical'],
                      blast_self_controls_pass=True,protein_gene_mapping_pass=True,domain_coordinates_pass=True,
                      count_unit='unique GTF gene locus',primary_representative='longest eligible isoform; lexical tie-break',
                      internal_stop_scan_policy='Replace internal * by X for HMMER only; retain flags and exclude these isoforms from counts.',
                      partial_model_flag='Missing annotated start_codon or stop_codon; review flag, not proof of biological truncation.',
                      blast_coverage='Union of all reported HSP intervals separately on query and subject; not sum of overlaps.',
                      outcome_blind=True,association_tests_run=False,
                      pilot_gate='pass' if positive_ok and specificity_ok and hydra_ok and all(s['unresolved_genes']==0 for s in summaries) else 'review_required',
                      source_script_sha256=sha(Path(__file__)),
                      output_sha256={str(f.relative_to(root)):sha(f) for f in sorted(out.glob('*.tsv'))})
    (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print('PILOT COMPLETE:',validation['pilot_gate'],flush=True)


if __name__ == '__main__':
    main()
