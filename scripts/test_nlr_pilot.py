"""Boundary checks for classification, competing repeats and BLAST coverage."""
import json
import tempfile
import unittest
from pathlib import Path
from run_nlr_pilot import architecture, overlap_fraction, parse_domains, union_length, parse_gtf
import gzip

S = dict(core_model='PF05729.20',broad_core_model='PF00931.30',lrr_models=['PF13516.14'],minimum_model_coverage=0.5)


def hit(model,start,end,retained=True):
    return dict(model=model,ali_from=start,ali_to=end,retained=retained,domain_name=model)


class PilotTests(unittest.TestCase):
    def test_downstream_must_be_strict_and_in_same_protein(self):
        n=hit(S['core_model'],100,200)
        self.assertFalse(architecture([n,hit('PF13516.14',190,210)],S)['canonical'])
        self.assertFalse(architecture([n,hit('PF13516.14',200,220)],S)['canonical'])
        self.assertTrue(architecture([n,hit('PF13516.14',201,220)],S)['canonical'])
        self.assertFalse(architecture([n,hit('PF13516.14',201,220,False)],S)['canonical'])
        self.assertFalse(architecture([hit('PF13516.14',201,220)],S)['canonical'])

    def test_nb_arc_is_not_nacht(self):
        a=architecture([hit(S['broad_core_model'],100,200),hit('PF13516.14',210,230)],S)
        self.assertTrue(a['nb_arc']); self.assertFalse(a['canonical'])

    def test_overlapping_core_models_are_flagged(self):
        a=architecture([hit(S['core_model'],100,200),hit(S['broad_core_model'],110,190)],S)
        self.assertTrue(a['overlapping_core_models']); self.assertTrue(a['mixed_core_models'])

    def test_coverage_union_does_not_double_count(self):
        self.assertEqual(union_length([(1,100),(90,150),(200,210)]),161)
        self.assertEqual(union_length([(1,100),(5,20),(1,100)]),100)

    def test_ga_coverage_and_repeat_overlap_boundaries(self):
        m={'PF13516.14':dict(length='24',gathering_threshold='23.4 8.7;')}
        lines=[]
        # Format follows hmmsearch domtblout, with independent sequence/domain scores.
        for number,(score,frm,to,left,right) in enumerate([(9,1,12,10,21),(10,1,12,12,23),(11,1,11,40,50),(8.6,1,24,60,83)],1):
            lines.append(f'p - 100 LRR_6 PF13516.14 24 1e-8 30 0 {number} 4 1e-6 1e-6 {score} 0 {frm} {to} {left} {right} {left} {right} 0.9 -')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'test.dom';p.write_text('\n'.join(lines))
            rows,_=parse_domains(p,m,{'p':'A'*100},S)
        self.assertEqual([r['retained'] for r in rows],[False,True,False,False])
        self.assertEqual(rows[0]['rejection_reason'],'competing_repeat_overlap_gt_0.5')
        self.assertEqual(rows[1]['model_coverage'],0.5)
        self.assertEqual(rows[2]['rejection_reason'],'model_coverage_below_0.5')
        self.assertEqual(rows[3]['rejection_reason'],'below_GA')

    def test_mapping_uses_explicit_gtf_attributes_not_id_suffix(self):
        lines=['ctg\tA\tgene\t1\t90\t.\t+\t.\tlocus42',
               'ctg\tA\ttranscript\t1\t90\t.\t+\t.\tproteinX',
               'ctg\tA\tCDS\t1\t90\t.\t+\t0\ttranscript_id "proteinX"; gene_id "locus42";']
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'test.gtf.gz'
            with gzip.open(p,'wt') as f:f.write('\n'.join(lines))
            mapping,genes,_=parse_gtf(p)
        self.assertEqual(mapping,{'proteinX':'locus42'});self.assertEqual(genes,{'locus42'})


if __name__=='__main__':unittest.main(verbosity=2)
