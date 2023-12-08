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
        