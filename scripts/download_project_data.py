import concurrent.futures
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / 'config/source_metadata'
RESERVE = 3 * 1024**3

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def main():
    jobs = []
    species = []
    for group, repo in [('genome', 'sebametz/coral-genomes-annotation'), ('microbiome', 'zaneveld/GCMP_Australia_Coevolution')]:
        tree = json.loads((META / (group + '_tree.json')).read_text())['tree']
        commit = json.loads((META / (group + '_commit.json')).read_text())['sha']
        genome_dirs = {str(Path(f['path']).parent) for f in tree if f['path'].endswith('/braker.aa.gz')}
        for f in tree:
            if f['type'] != 'blob':
                continue
            p = Path(f['path'])
            dest = None
            if group == 'genome':
                if str(p.parent) in genome_dirs and (p.name in ['README.md', 'braker.aa.gz', 'braker.gtf.gz'] or p.name.endswith(('.stats', '_busco.jpeg', '.emapper.annotations.gz'))):
                    dest = Path('data/raw/genomes/metz_v1') / p.relative_to('results')
                elif str(p) in ['README.md', 'LICENSE', 'methods/README.md', 'methods/other.md']:
                    dest = Path('data/raw/genomes/metz_v1/source_docs') / p
            else:
                base = '1_canonical_starting_files/input/'
                if str(p) in [base + x for x in ['gcmp16S_map_r25.txt', 'gcmp16S_map_r25_with_mitochondrial_data.txt', 'otu_table_mc2_wtax_no_pynast_failures_no_organelles.biom']]:
                    dest = Path('data/raw/microbiome/pollock2018') / p.name
                elif str(p) in [base + 'host_tree_from_step_11.newick', base + 'gg_constrained_rep_set_fastttree.tre']:
                    folder = 'phylogeny' if p.suffix == '.newick' else 'microbiome'
                    dest = Path('data/raw') / folder / 'pollock2018' / p.name
                elif str(p) == '0_otu_table_generation/input/GCMP_16S_samples_sequenced.txt':
                    dest = Path('data/raw/microbiome/pollock2018') / p.name
                elif str(p) == 'README.md' or (str(p).startswith('0_otu_table_generation/procedure/') and p.suffix in ['.sh', '.r']):
                    dest = Path('data/raw/microbiome/pollock2018/source_docs') / p.name
            if dest:
                jobs.append(dict(dataset=group, source_commit=commit, source_path=str(p), local_path=str(dest), expected_bytes=f['size'], git_blob_sha1=f['sha'], url='https://raw.githubusercontent.com/' + repo + '/' + commit + '/' + str(p)))
        if group == 'genome':
            for dirname in sorted(genome_dirs):
                tolid, name = Path(dirname).name.split('_', 1)
                species.append(dict(genome_id=tolid, source_species=name.replace('_',' '), annotation_version='v1.0.0', source_commit=commit, protein_path=str(Path('data/raw/genomes/metz_v1') / Path(dirname).name / 'braker.aa.gz'), annotation_path=str(Path('data/raw/genomes/metz_v1') / Path(dirname).name / 'braker.gtf.gz'), taxonomy_status='source name; harmonization pending'))
    write_tsv(ROOT/'config/species_manifest.tsv', species, list(species[0]))
    fields = list(jobs[0])
    write_tsv(ROOT/'config/planned_downloads.tsv', jobs, fields)
    total = sum(j['expected_bytes'] for j in jobs)
    missing = sum(j['expected_bytes'] for j in jobs if not (ROOT / j['local_path']).exists())
    free = shutil.disk_usage(ROOT).free
    print(json.dumps(dict(files=len(jobs), genome_sets=len(species), planned_bytes=total, missing_bytes=missing, free_bytes=free, reserve_bytes=RESERVE)), flush=True)
    if missing and missing + RESERVE > free:
        raise RuntimeError('Insufficient disk space for planned download plus reserve')

    def fetch(job):
        row = dict(job)
        dest = ROOT / job['local_path']
        dest.parent.mkdir(parents=True, exist_ok=True)
        part = dest.with_name(dest.name + '.part')
        try:
            if not dest.exists():
                proc = subprocess.run(['curl', '-L', '--fail', '--retry', '2', '--retry-delay', '2', '--connect-timeout', '20', '--max-time', '240', '-sS', job['url'], '-o', str(part)], capture_output=True, text=True)
                if proc.returncode:
                    raise RuntimeError(proc.stderr.strip())
                candidate = part
            else:
                candidate = dest
            size = candidate.stat().st_size
            gitsha = hashlib.sha1(('blob ' + str(size) + '\0').encode())
            sha256 = hashlib.sha256()
            with candidate.open('rb') as source:
                for block in iter(lambda: source.read(1024*1024), b''):
                    gitsha.update(block)
                    sha256.update(block)
            if size != job['expected_bytes'] or gitsha.hexdigest() != job['git_blob_sha1']:
                raise RuntimeError('Size or Git blob checksum mismatch')
            if dest.suffix == '.gz':
                with gzip.open(candidate, 'rb') as source:
                    for block in iter(lambda: source.read(1024*1024), b''):
                        pass
            if candidate == part:
                part.rename(dest)
            row.update(status='verified', actual_bytes=size, sha256=sha256.hexdigest(), checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), error='')
        except Exception as exc:
            row.update(status='failed', error=str(exc))
        return row

    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for row in pool.map(fetch, jobs):
            rows.append(row)
            if len(rows) % 20 == 0 or row['status'] != 'verified':
                print(json.dumps(dict(completed=len(rows), total=len(jobs), last=row['local_path'], status=row['status'], error=row.get('error',''))), flush=True)
            write_tsv(ROOT/'config/download_manifest.tsv', rows, fields + ['status','actual_bytes','sha256','checked_utc','error'])
    summary = dict(files=len(rows), verified=sum(r['status']=='verified' for r in rows), failed=sum(r['status']!='verified' for r in rows), downloaded_bytes=sum(r.get('actual_bytes',0) for r in rows), free_bytes=shutil.disk_usage(ROOT).free)
    print(json.dumps(summary), flush=True)
    if summary['failed']:
        raise SystemExit(1)

if __name__ == '__main__':
    main()
