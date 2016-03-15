#!/bin/python

import numpy as np
import random as rand
import ipdb
import copy

class Sample(object):
    """
    Takes a dataset of DNA samples and generates a sample from their possible phylogenies.

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
        distribution_vector = np.empty((len(self.multiplicity_vector)))
        M = self.singleton_alleles()
        for allele in range(len(self.multiplicity_vector)):
            n_k = self.multiplicity_vector[allele]
            # Multiplicity >= 2 or the allele is still in the sample
            if (allele in M) or (n_k >= 2):
                distribution_vector[allele] = float(n_k)
            else: # n_k = 1 and k not in M
                distribution_vector[allele] = 0
        total_mutation_sites = sum(distribution_vector)
        return map(lambda count: count/float(total_mutation_sites), distribution_vector)

    def generate_distribution_Hobolth(self):
            """
            :param: allele: the allele, k, a number
            :return: An array of probabilities for each allele, according to the Hobolth proposal distribution
            """
            distribution_vector = np.empty((len(self.multiplicity_vector)))
            M = self.singleton_alleles()
            theta = self.watterson_est()

            #Row vector times sequence matrix
            d_m = np.dot(self.multiplicity_vector,self.sequence_array)
            
            u_km = copy.copy(self.sequence_array) #temporarily copy sequence_array to preserve dimensions
            for row in range(np.shape(self.sequence_array)[0]):
                for column in range(np.shape(self.sequence_array)[1]):
                    if(int(sequence_array[row][column]) == 1):
                        u[row][column] = p_theta(theta, d_m[column]) * self.multiplicity_vector[row] / float(d_m[column])
                    else: #if sequence_array[row][column] == 0 WILL NEED TO CHANGE FOR 2 LOCI CASE
                        u[row][column] = ((1.0 - p_theta(theta, d_m[column]))*self.multiplicity_vector[row])/(float(sum(self.multiplicity_vector)) - float(d_m[column]))
            
            for allele in range(len(self.multiplicity_vector)):
                ### Generate the probability for each allele according to Hobolth.
                ### distribution_vector[allele] = probability
                n_k = self.multiplicity_vector[allele]
                if (allele in M) or (n_k >= 2):
                    distribution_vector[allele] = float(sum([u_km[allele][column] for column in range(np.shape(u_km)[1])]))
                else:
                    distribution_vector[allele] = 0
                
            total_mutation_sites = sum(distribution_vector)
            return map(lambda count: count/float(total_mutation_sites), distribution_vector)
    
    def watterson_est(self):
        M = np.shape(self.sequence_array)[1] #number of columns        
        return float(M)/sum([1/float(i) for i in range(1,sum(self.multiplicity_vector)-1)])
        
    def p_theta(self, theta, d):
        n = sum(self.multiplicity_vector)
        p = 0
        for k in range(2, n - d + 1):
            f1 = (float(d) - 1)/(n - k)
            f2 = 1.0/float(k-1+theta)
            f3 = stat.comb(n-d-1,k-2)
            f4 = stat.comb(n-1, k-1)
            f5 = 0
            for i in range(2,n-d+1):
                t1 = 1.0/ float(i-1+theta)
                t2 = stat.comb(n-d-1,i-2)
                t3 = stat.comb(n-1,i-1)
                f5 += t1 * t2 / float(t3)
                
            p += f1 * f2 * f3/(f4 * f5)
            
            return p
        

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
                        break
                new_allele_sequence = copy.copy(self.sequence_array[allele])
                new_allele_sequence[deleted_site] = 0
                coalesce_partner = self.contains_allele(new_allele_sequence)
                # Actually delete the column now that changes have been made.
                # If s_k^m does occur in the sample then merge k with whatever else has this
                if coalesce_partner != -1:
                    # Don't coalesce if it's all 0's
                    if np.all(new_allele_sequence == 0) and not np.all(self.sequence_array == 0):
                        self.multiplicity_vector[allele] -= 1
                    else:
                        # Delete the duplicate row and adjust the multiplicities
                        self.multiplicity_vector[coalesce_partner] += copy.copy(self.multiplicity_vector[allele])
                        self.multiplicity_vector[allele] = 0
                        self.sequence_array[allele, :] = 0
                # If s_k^m does not occur in the sample then don't change the multiplicity
                else:
                    # Already deleted the column
                    pass
                self.sequence_array[:, deleted_site] = 0

        else: # n_k = 1 and k not in M
            pass

        return None

    def generate_sample(self, proposal="Hobolth"):
        """
        Starts from a random allele k and generates a sample based on the SD distribution.
        :return: An evolutionary history, with its weight.
        """
        path = [ ]
        while True:
            # Recur until the sequence matrix is all zeroes and there is only one non-zero entry in the multiplicity vector
            if np.all(self.sequence_array == 0) and sum(self.multiplicity_vector) <= 1:
                path.append((self.sequence_array.tolist(), self.multiplicity_vector.tolist(), 1))
                for state in path:
                    print state[1]
                    print_matrix(state[0])
                    print "Weight %f" % state[2]
                    print "---"
                return path
            distribution = None
            if proposal == "Hobolth":
                distribution = self.generate_distribution_SD()
            elif proposal == "SD":
                distribution = self.generate_distribution_Hobolth()
            else:
                raise Exception("No such distribution %s exists" % proposal)
            # Choose an event with this distribution
            mutation_site = np.random.choice(len(self.multiplicity_vector), 1, p=distribution)[0]
            path.append((self.sequence_array.tolist(), self.multiplicity_vector.tolist(), distribution[mutation_site]))
            # Update the allele_array and multiplicity_vector
            self.update_sequences_and_counts(mutation_site)


class Sampler(object):
    """
    Generates N samples, using the Sample object.
    """

    def __init__(self, sequence_array, multiplicity_vector):
        self.sequence_array = np.array(sequence_array)
        self.multiplicity_vector = np.array(multiplicity_vector)

    def generate_N_samples(self, n, proposal_dist="Hobolth"):
        """
        Generates n samples according to the desired proposal distribution. Hobolth default.
        :param: n, number of samples
        :kwargs: proposal, {Hobolth, SD}
        :return: N samples in a list where EACH sample is 
        [
            (state tuple
                sequence_array,
                multiplicity_vector,
                weight
            ),
            ...
        ]
        """
        samples = [ ]
        for i in range(n):
           sample = Sample(copy.copy(self.sequence_array), copy.copy(self.multiplicity_vector)).generate_sample()
           samples.append(sample)
        return samples


    def calculate_likelihood(self, samples):
        """
        Calculates the likelihood of a set of states generated by the Sampler.
        """
        from operator import mul
        total_weight = 0
        for sample in samples:
            # Take the product over the weights, which is the third entry of each tuple in the list of states.
            total_weight += reduce(mul, map(lambda state_tuple: state_tuple[2], sample), 1)
        return total_weight/len(samples)

def print_matrix(matrix):
    for row in matrix:
        print row

toy_sequence_array = np.array([
    [1, 0, 0, 0, 0],
    [0, 1, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 1, 1]
])
toy_multiplicity_vector = np.array([3, 1, 1, 1])
toy_sampler = Sampler(toy_sequence_array, toy_multiplicity_vector)

"""
states = toy_sampler.generate_sample()
for state in states:
    print state[1]
    print_matrix(state[0])
    print "Weight %f" % state[2]
    print "---"
"""

samples = toy_sampler.generate_N_samples(100, proposal_dist="Hobolth")
print toy_sampler.calculate_likelihood(samples)
