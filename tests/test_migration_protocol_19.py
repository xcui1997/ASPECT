import unittest
import numpy as np
from aspect.migration_test import MigrationTester
class ProtocolTests(unittest.TestCase):
 def test_exact_subsamples_reproducible(self):
  t=np.arange(103,dtype=float);c=np.c_[t,t*.1,np.zeros(len(t))];m=MigrationTester({'migration_test_params':{}})
  a=m._perturbations(c,t,np.random.default_rng(321));b=m._perturbations(c,t,np.random.default_rng(321))
  self.assertEqual(len(a),18)
  for x,y in zip(a,b):
   self.assertEqual(len(x[0]),92);self.assertEqual(x[2:],(70,10));np.testing.assert_array_equal(x[0],y[0]);self.assertEqual(len(np.unique(x[1])),92);self.assertTrue(np.all(np.diff(x[1])>0))
 def test_reject_legacy_mixed_vote(self):
  for p in [{'n_subsample':14},{'origin_n_first':[8,12]},{'percentiles':[75,85]},{'drop_rate':.2}]:
   with self.assertRaises(ValueError):MigrationTester({'migration_test_params':{'perturbation':p}})
if __name__=='__main__':unittest.main()
