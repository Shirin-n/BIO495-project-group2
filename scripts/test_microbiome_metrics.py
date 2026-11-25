import unittest
import numpy as np
import pandas as pd
from prepare_microbiome_metrics import alpha,rarefy,aggregate,pcoa
from scipy.spatial.distance import pdist,squareform

class DiversityTests(unittest.TestCase):
    def test_alpha_known_distributions(self):
        np.testing.assert_allclose(alpha([5,5,0])[0],[np.log(2),2,2])
        np.testing.assert_allclose(alpha([0,10])[0],[0,1,1])
    def test_without_replacement_at_full_depth(self):
        observed=rarefy([1,2,7],10,np.random.SeedSequence(7),100)
        np.testing.assert_allclose(observed,np.tile(alpha([1,2,7]),(100,1)))
    def test_rarefaction_deterministic_and_bounded(self):
        a=rarefy([2,5,20],8,np.random.SeedSequence(19),100)
        b=rarefy([2,5,20],8,np.random.SeedSequence(19),100)
        np.testing.assert_array_equal(a,b)
        self.assertTrue((a[:,1]<=3).all())
        with self.assertRaises(AssertionError):rarefy([2,5,20],28,np.random.SeedSequence(19))
    def test_colonies_have_equal_weight_not_sample_count(self):
        values=pd.DataFrame({'x':[0.,0.,1.]},index=['s1','s2','s3'])
        md=pd.DataFrame({'canonical_species':['A']*3,'colony_id':['c1','c1','c2']},index=values.index)
        colony,species=aggregate(values,md)
        self.assertEqual(species.loc['A','x'],.5)
        self.assertNotEqual(species.loc['A','x'],values.x.mean())
        _,reordered=aggregate(values.iloc[::-1],md)
        pd.testing.assert_frame_equal(species,reordered)
    def test_euclidean_pcoa_reconstructs_distances(self):
        d=squareform(pdist([[0,0],[1,0],[0,2],[2,3]]))
        eig,prop,coords,info=pcoa(d)
        np.testing.assert_allclose(squareform(pdist(coords)),d,atol=1e-10)
        self.assertEqual(info['negative_eigenvalues'],0)
        self.assertAlmostEqual(prop.sum(),1)
    def test_non_euclidean_pcoa_reports_negative_spectrum(self):
        d=np.array([[0,.2,.2],[.2,0,1],[.2,1,0]])
        original=d.copy();eig,prop,coords,info=pcoa(d)
        self.assertGreater(info['negative_eigenvalues'],0)
        self.assertGreater(info['negative_fraction_absolute_inertia'],0)
        np.testing.assert_array_equal(d,original)
        self.assertEqual(info['distance_correction'],'none')

if __name__=='__main__':unittest.main(verbosity=2)
