from pathlib import Path
import csv, gzip, hashlib, json, shutil

root=Path(__file__).resolve().parents[1]
manifest=list(csv.DictReader((root/'config/download_manifest.tsv').open(),delimiter='\t'))
assert len(manifest)==255
assert all(r['status']=='verified' for r in manifest)
genomes=root/'data/raw/genomes/metz_v1'
proteins=list(genomes.glob('*/braker.aa.gz'))
gtfs=list(genomes.glob('*/braker.gtf.gz'))
assert len(proteins)==len(gtfs)==40
for p in proteins:
    with gzip.open(p,'rt') as f:
        assert f.readline().startswith('>'),str(p)
micro=root/'data/raw/microbiome/pollock2018'
biomfile=micro/'otu_table_mc2_wtax_no_pynast_failures_no_organelles.biom'
with biomfile.open('rb') as f:
    signature=f.read(8)
sample_ids=None
if signature == b'\x89HDF\r\n\x1a\n':
    biom_summary=dict(format='HDF5 BIOM',validation='Git source checksum and HDF5 signature verified',internal_validation='full matrix validation is performed by the analysis notebook')
else:
    with biomfile.open() as f:
        biom=json.load(f)
    assert biom['shape']==[len(biom['rows']),len(biom['columns'])]
    assert biom['matrix_type']=='sparse'
    assert all(0<=r<biom['shape'][0] and 0<=c<biom['shape'][1] and v>=0 for r,c,v in biom['data'])
    sample_ids={r['id'] for r in biom['columns']}
    biom_summary=dict(format=biom['format'],shape=biom['shape'],nonzero_entries=len(biom['data']))
metadata={}
for p in micro.glob('gcmp16S_map*.txt'):
    with p.open() as f:
        reader=csv.DictReader(f,delimiter='\t'); records=list(reader)
        first=reader.fieldnames[0]
        ids={r[first] for r in records}
        metadata[p.name]=dict(rows=len(records),fields=reader.fieldnames,unique_sample_ids=len(ids),biom_samples_matched=len(ids & sample_ids) if sample_ids is not None else None,biom_samples_missing=len(sample_ids-ids) if sample_ids is not None else None)
host_tree=root/'data/raw/phylogeny/pollock2018/host_tree_from_step_11.newick'
assert host_tree.read_text().strip().endswith(';')
summary=dict(genome_protein_sets=len(proteins),genome_gtf_sets=len(gtfs),main_files_verified=len(manifest),main_download_bytes=sum(int(r['actual_bytes']) for r in manifest),biom=biom_summary,metadata=metadata,free_bytes=shutil.disk_usage(root).free,scope='Download integrity and basic format checks only; BIOM internal parsing, species overlap, health filtering and biological validation pending.')
(root/'results/qc').mkdir(parents=True, exist_ok=True)
(root/'results/qc/download_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
