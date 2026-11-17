"""Classification is invariant to time translation and input row order."""
import unittest
import numpy as np
import pandas as pd
from aspect.classifier import ClusterClassifier


class TimeOriginTests(unittest.TestCase):
    def setUp(self):
        self.classifier = ClusterClassifier({'classification_params':{
            'use_features':['mean_t','skew','mag_diff'],
            'features_greater_than':['mean_t'], 'features_less_than':['skew','mag_diff'],
            'mean_t':.32,'skew':6,'mag_diff':.5}})
        t=np.arange(20,dtype=float)**1.4
        self.df=pd.DataFrame({'stime':t*86400,'time_day':t,'mag':np.r_[3.,np.repeat(1.,19)]})

    def test_parent_time_offset_and_row_order_do_not_change_features(self):
        baseline=self.classifier.calculate_metrics(self.df)
        shifted=self.df.copy();shifted['stime']+=1e9;shifted['time_day']+=1000
        changed=self.classifier.calculate_metrics(shifted.sample(frac=1,random_state=9))
        for key in ('mean_t','skew','mean_t_index','skew_index','mag_diff'):
            self.assertAlmostEqual(baseline[key],changed[key],places=9)

    def test_time_day_only_is_rezeroed(self):
        a=self.classifier.calculate_metrics(self.df.drop(columns='stime'))
        b=self.df.drop(columns='stime').copy();b['time_day']+=1000
        b=self.classifier.calculate_metrics(b)
        self.assertAlmostEqual(a['mean_t'],b['mean_t'],places=12)

if __name__=='__main__':unittest.main()
