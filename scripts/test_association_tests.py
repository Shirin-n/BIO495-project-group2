import unittest
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist,squareform
from run_association_tests import exact_permutations,permutations,test_endpoints,holm_two,rank_effect

class AssociationTests(unittest.TestCase):
    def fixture(self,n=5):
        rng=np.random.default_rng(12)
        profiles=rng.dirichlet([1,2,3],size=n)
        return squareform(pdist(profiles,metric='braycurtis'))
    def test_ranks_with_ties_match_scipy(self):
        x=[2,2,4,8,9];y=[.4,.6,.1,.7,.2]
        self.assertAlmostEqual(rank_effect(x,y)[0],spearmanr(x,y).statistic)
    def test_exact_both_endpoints_match_independent_brute_force(self):
        x=np.array([1,1,3,5,7]);y=np.array([.8,.1,.4,.3,.9]);d=self.fixture()
        perm=exact_permutations(5);rows,nulls,gene=test_endpoints(x,y,d,perm,'exact')
        triu=np.triu_indices(5,1)
        expected_a=np.array([spearmanr(x[p],y).statistic for p in perm])
        expected_b=np.array([spearmanr(np.abs(np.log1p(x[p,None])-np.log1p(x[p][None,:]))[triu],d[triu]).statistic for p in perm])
        for r,expected in zip(rows,[expected_a,expected_b]):
            np.testing.assert_allclose(nulls[r['endpoint']]/r['rank_norm_product'],expected,atol=1e-14)
            extreme=int((abs(expected)>=abs(expected[0])-1e-12).sum())
            self.assertEqual(r['extreme_count'],extreme)
            self.assertEqual(r['p_raw'],extreme/120)
    def test_perfect_rank_association_exact_includes_identity(self):
        rows,_,_=test_endpoints(np.arange(5),np.arange(5),self.fixture(),exact_permutations(5),'exact')
        self.assertEqual(rows[0]['extreme_count'],2)
        self.assertAlmostEqual(rows[0]['p_raw'],2/120)
    def test_undefined_is_not_zero_or_evidence_for_null(self):
        rows,_,_=test_endpoints(np.ones(5),np.arange(5),self.fixture(),exact_permutations(5),'exact')
        for r in rows:
            self.assertEqual(r['status'],'undefined');self.assertTrue(np.isnan(r['p_raw']));self.assertTrue(np.isnan(r['effect_rho']))
    def test_low_species_is_descriptive_only(self):
        rows,_,_=test_endpoints(np.arange(3),np.arange(3),self.fixture(3))
        self.assertEqual(rows[0]['effect_rho'],1)
        self.assertEqual(rows[0]['status'],'descriptive_only');self.assertTrue(np.isnan(rows[0]['p_raw']))
    def test_holm_reserves_two_endpoints(self):
        np.testing.assert_allclose(holm_two([.03,.04]),[.06,.06])
        np.testing.assert_allclose(holm_two([np.nan,.03]),[np.nan,.06],equal_nan=True)
        self.assertTrue(np.isnan(holm_two([np.nan,np.nan])).all())
    def test_monte_carlo_seed_and_plus_one(self):
        p,mode=permutations(9,np.random.SeedSequence(17),99)
        q,_=permutations(9,np.random.SeedSequence(17),99)
        np.testing.assert_array_equal(p,q)
        rows,_,_=test_endpoints(np.arange(9),np.arange(9),self.fixture(9),p,mode)
        for r in rows:self.assertEqual(r['p_raw'],(r['extreme_count']+1)/100)

if __name__=='__main__':unittest.main(verbosity=2)
