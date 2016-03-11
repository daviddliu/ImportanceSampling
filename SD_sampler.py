#!/bin/python

import numpy as np
import random as rand
import ipdb

class Sampler(object):
    """
    Takes a dataset of DNA samples and samples from their possible phylogenies.

    in:
     + sequence_array: list of lists
        type: np.array
         where the rows are allele k, and the columns are mutation site m

     + multiplicity_vector: vector of multiplicities
        type: np.array
         where the index is allele k, and the value n_k is the count of allele k
    """
    def __init__(self, sequence_array, multiplicity_vector):
        self.sequence_array = np.array(sequence_array)
        self.multiplicity_vector = np.array(multiplicity_vector)


    def contains_allele(self, allele):
        """
        Allele is a numpy array with mutation sites.
        :return: -1 if no, index if yes
        """
        for i, row in enumerate(self.sequence_array):
            if np.array_equal(allele, row):
                return i
        return -1

    def delete_site_from_allele(self, allele):
        for i, num in enumerate(allele):
            if num == 1:
                new_allele = allele
                allele[i] = 0
                return allele


    def singleton_alleles(self):
        """
        This is the "M" in the paper. Set of row indices with alleles bearing at least one singleton mutation (a mutation borne by one allele).
        """
        singleton_alleles = set()
        for column in self.sequence_array.T:
            if sum(column) == 1:
                # The first nonzero element of that column: this is the allele we want
                singleton_alleles.add(np.amin(np.nonzero(column)))
        return singleton_alleles

        
    def is_unique_mutation(self, allele, mutation_index):
        """
        Returns True if this mutation is the only one in the column
        """
        return sum(self.sequence_array.T[mutation_index]) == 1


    def generate_distribution_SD(self):
        """
        :param allele: the allele k, a number
        :return: An array of probabilities for each allele, according to the Stephen-Donnelly proposal distribution
        """
        total_mutation_sites = sum(self.multiplicity_vector)
        distribution_vector = np.empty((len(self.multiplicity_vector)))
        M = self.singleton_alleles()
        for allele in range(len(self.multiplicity_vector)):
            n_k = self.multiplicity_vector[allele]
            # Multiplicity >= 2 or the allele is still in the sample
            if (allele in M) or (n_k >= 2):
                distribution_vector[allele] = float(n_k)/float(total_mutation_sites)
            else: # n_k = 1 and k not in M
                distribution_vector[allele] = 0

        return distribution_vector

    

    def update_sequences_and_counts(self, allele):
        """
        Given an allele action, update the sequence matrix and multiplicity vector.
        """
        M = self.singleton_alleles()
        n_k = self.multiplicity_vector[allele]
        # Multiplicity >= 2 or the allele is still in the sample
        if (allele in M) or (n_k >= 2):
            # Update the array and counts
            if n_k >= 2: # Multiple copies
                # Subtract 1 from n_k
                self.multiplicity_vector[allele] -= 1
            else: # Allele is in M
                # Delete the first column with a nonzero entry (but it must be a unique mutation)
                mutation_indices = np.nonzero(self.sequence_array[allele])[0]
                deleted_site = None
                for mutation in mutation_indices:
                    if self.is_unique_mutation(allele, mutation):
                        # Delete this column
                        deleted_site = mutation
                        self.sequence_array[:, deleted_site] = 0
                        print "Deleted a column"
                        break

                new_allele_sequence = self.sequence_array[allele]
                new_allele_sequence[deleted_site] = 0
                coalesce_partner = self.contains_allele(new_allele_sequence)
#                ipdb.set_trace()
                # If s_k^m does occur in the sample then merge k with whatever else has this
                if coalesce_partner != -1:
                    self.sequence_array[allele, :] = 0
                    self.multiplicity_vector[allele] = 0
                    self.multiplicity_vector[coalesce_partner] += 1
                # If s_k^m does not occur in the sample then don't change the multiplicity
                else:
                    # Already deleted the column
                    pass

        else: # n_k = 1 and k not in M
            pass

        return None

    def random_choice(self, array):
        rand_num = rand.uniform(0, 1)
        tally = 0
        for i, num in enumerate(array):
            tally += num
            if tally > rand_num:
                return i
        return len(array) - 1


    def generate_sample(self):
        """
        Starts from a random allele k and generates a sample based on the SD distribution.
        :return: An evolutionary history, with its weight.
        """
        path = []
        while True:
            print self.multiplicity_vector
            distribution = self.generate_distribution_SD()
            # Choose an event with this distribution
            # mutation_site = np.random.choice(len(np.nonzero(self.multiplicity_vector)), 1, p=distribution)[0]
            mutation_site = self.random_choice(distribution)
            path.append(self.sequence_array[mutation_site].tolist())
            # Update the allele_array and multiplicity_vector
            self.update_sequences_and_counts(mutation_site)
            # Recur until the sequence matrix is all zeroes
            if np.all(self.sequence_array == 0) and (np.count_nonzero(self.multiplicity_vector) == len(self.multiplicity_vector) - 1):
                return None
        

toy_sequence_array = np.array([
    [1, 0, 0, 0, 0],
    [0, 1, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 1, 1]
])
toy_multiplicity_vector = np.array([3, 1, 1, 1])
toy_sampler = Sampler(toy_sequence_array, toy_multiplicity_vector)

print toy_sampler.generate_sample()
