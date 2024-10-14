import numpy as np
import pandas as pd

from  numpy.random import MT19937
from numpy.random import RandomState, SeedSequence

class ProteinProfileGenerator():
    def __init__(self, peptide_profiles):
        self._peptide_profiles = peptide_profiles
        
        self.protein_profile_dataframe = None
        self._generate_protein_profile_dataframe()

    def _generate_protein_profile_dataframe(self):
        collected_profiles = [x.peptide_profile_vector for x in self._peptide_profiles]
        protnames_for_index = [x.protein_name for x in self._peptide_profiles]
        pepnames_for_index = [f'{idx}' for idx in range(len(self._peptide_profiles))]
        self.protein_profile_dataframe = pd.DataFrame(collected_profiles,index=[protnames_for_index, pepnames_for_index])
        self.protein_profile_dataframe.index.names = ['protein', 'ion']
        self.protein_profile_dataframe = np.log2(self.protein_profile_dataframe.replace(0, np.nan))



class PeptideProfile():
    def __init__(self, protein_name, fraction_zeros_in_profile, systematic_peptide_shift, add_noise, num_samples = 20, min_intensity = 1e6, max_intensity = 1e10):


        self._fraction_zeros_in_profile = fraction_zeros_in_profile
        self._systematic_peptide_shift = systematic_peptide_shift
        self._add_noise = add_noise
        self._min_intensity = min_intensity
        self._max_intensity = max_intensity
        self._num_samples = num_samples

        self.protein_name = protein_name
        self.peptide_profile_vector = []
        self._define_peptide_profile_vector()

    def _define_peptide_profile_vector(self):
        self.peptide_profile_vector = self._get_single_peptide_profile_template()
        self._scale_profile_vector()
        if self._add_noise:
            self._apply_poisson_noise_to_profilevector()
        self._add_zeros_to_profilevector()

    def _get_single_peptide_profile_template(self):
        rs = RandomState(MT19937(SeedSequence(42312)))
        return rs.randint(low=self._min_intensity, high=self._max_intensity,size=self._num_samples)

    def _scale_profile_vector(self):
        self.peptide_profile_vector = self.peptide_profile_vector*self._systematic_peptide_shift

    def _apply_poisson_noise_to_profilevector(self):
        self.peptide_profile_vector = np.random.poisson(lam=self.peptide_profile_vector, size=len(self.peptide_profile_vector))

    def _add_zeros_to_profilevector(self):
        num_elements_to_set_zero = int(self._num_samples*self._fraction_zeros_in_profile)
        idxs_to_set_zero = np.random.choice(self._num_samples,size=num_elements_to_set_zero, replace=False)
        self.peptide_profile_vector[idxs_to_set_zero] = 0
        




class RatioChecker():
    def __init__(self, formatted_df : pd.DataFrame, organism2expectedfc : dict, organism2CI95 : dict):
        """the ratio checker takes in a dataframe with columns 'log2fc' and 'organism' and checks if the ratios are consistent with the expected median fold changes and confidence intervals
        
        Args:
            formatted_df (pd.DataFrame): dataframe with columns 'log2fc' and 'organism' from mixed species experiment
            organism2expectedfc (dict): the expected log2fc for this organism in the mixed species experiment
            organism2CI95 (dict): the expected confidence interval for this organism in the mixed species experiment (i.e. 95% of protein ratios should be within this interval)
        """
        self._formatted_df = formatted_df
        self._organism2expectedfc = organism2expectedfc
        self._organism2deviation_threshold = organism2CI95
        self._organism2fcs = self._get_organism2fcs()
        self._check_ratio_consistency_per_organism()

    def _get_organism2fcs(self):
        return self._formatted_df.groupby('organism')['log2fc'].apply(list).to_dict()
    
    def _check_ratio_consistency_per_organism(self):
        for organism in self._organism2expectedfc.keys():
            expected_fc = self._organism2expectedfc[organism]
            deviation_threshod = self._organism2deviation_threshold[organism]
            fcs = self._organism2fcs[organism]
            RatioConsistencyChecker(fcs = fcs, expected_fc = expected_fc, deviation_threshold = deviation_threshod)			

class RatioConsistencyChecker():
    def __init__(self, fcs, expected_fc, deviation_threshold):
        self._fcs = np.array([fc for fc in fcs if np.isfinite(fc)])
        self._expected_fc = expected_fc
        self._fc_deviations = self._calc_fc_deviations()
        self._deviation_threshold = deviation_threshold
        self._fc_deviation_center = self._calc_deviation_center() #should be 0 if no bias
        self._fraction_consistent = self._get_fraction_consistent()

        self._assert_no_bias()

    def _calc_fc_deviations(self):
        return abs(self._fcs - self._expected_fc)
    
    def _calc_deviation_center(self):
        return np.nanmedian(self._fc_deviations)
    
    def _get_fraction_consistent(self):
        consistent_fcs= sum(self._fc_deviations <= self._deviation_threshold)
        total_fcs = len(self._fcs)
        return consistent_fcs/total_fcs
    
    def _assert_no_bias(self):
        assert self._fc_deviation_center < self._deviation_threshold, f"Deviation from center: {self._fc_deviation_center:.2f}"
        assert self._fraction_consistent >0.95, f"Fraction consistent below 95%: {self._fraction_consistent:.2f}"
        print("Checked ratios, no bias detected")