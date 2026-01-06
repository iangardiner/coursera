import unittest
from typing import Dict

import motif


class MotifTest(unittest.TestCase):


    def test_find_motif_enumerations(self):

        for i in range(1, 7):
            with open(f"MotifEnumeration/inputs/input_{i}.txt") as input_file:
                params = input_file.readline().strip().split(" ")
                k = int(params[0])
                d = int(params[1])
                genomes = input_file.readline().strip().split(" ")

            result = motif.find_motif_enumerations(genomes, k, d)

            with open(f"MotifEnumeration/outputs/output_{i}.txt") as output_file:
                line = output_file.readline().strip()
                expected = set()
                if line:
                    expected.update(line.split(" "))

            print(f"Found {result} motif enumerations present in all {len(genomes)} genomes")
            self.assertEqual(expected, result)


    def test_big_motif_enumerations(self):

        with open("dataset_30302_8.txt") as input_file:
            params = input_file.readline().strip().split(" ")
            k = int(params[0])
            d = int(params[1])
            genomes = input_file.readline().strip().split(" ")

        result = motif.find_motif_enumerations(genomes, k, d)

        print(" ".join(result))


    def test_compute_entropy(self):

        motifs = [
            "AAAA",
            "CCCC",
            "GGGG",
            "TTTT"
        ]
        entropy = motif.compute_entropy(motifs)
        self.assertEqual(entropy, 8)

        motifs = [
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA"
        ]
        entropy = motif.compute_entropy(motifs)
        self.assertEqual(entropy, 0)

        motifs = [
            "AAAA",
            "AAAA",
            "TTTT",
            "TTTT"
        ]
        entropy = motif.compute_entropy(motifs)
        self.assertEqual(entropy, 4)

    def test_compute_motif_entropy(self):
        motifs = [
            "TCGGGGGTTTTT",
            "CCGGTGACTTAC",
            "ACGGGGATTTTC",
            "TTGGGGACTTTT",
            "AAGGGGACTTCC",
            "TTGGGGACTTCC",
            "TCGGGGATTCAT",
            "TCGGGGATTCCT",
            "TAGGGGAACTAC",
            "TCGGGTATAACC"
        ]
        entropy = motif.compute_entropy(motifs)
        print(f"motif entropy is {entropy}")


    def test_find_distances(self):

        for i in range(1, 5):
            with open(f"DistanceBetweenPatternAndStrings/inputs/input_{i}.txt") as input_file:
                pattern = input_file.readline().strip()
                genomes = input_file.readline().strip().split(" ")

            result = motif.find_distances_between_pattern_and_strings(pattern, genomes)

            with open(f"DistanceBetweenPatternAndStrings/outputs/output_{i}.txt") as output_file:
                expected = int(output_file.readline().strip())

            self.assertEqual(expected, result)


    def test_find_big_distances(self):

        with open("dataset_30312_1.txt") as input_file:
            pattern = input_file.readline().strip()
            genomes = input_file.readline().strip().split(" ")

        result = motif.find_distances_between_pattern_and_strings(pattern, genomes)

        print(result)


    def test_find_median_string(self):

        for i in range(1, 5):
            with open(f"MedianString/inputs/input_{i}.txt") as input_file:
                k = int(input_file.readline().strip())
                genomes = input_file.readline().strip().split(" ")

            result = motif.find_median_string(genomes, k)

            with open(f"MedianString/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip()

            print(f"median strings found are {result} compared to expected {expected}")
            self.assertIn(expected, result)


    def test_median_string_test(self):

        genomes = [
            "CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
            "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
            "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG",
        ]

        result = motif.find_median_string(genomes, 7)
        print(result)


    def test_big_median_string(self):

        with open("dataset_30304_9.txt") as input_file:
            k = int(input_file.readline().strip())
            genomes = input_file.readline().strip().split(" ")

        result = motif.find_median_string(genomes, k)

        print(result)


    def read_profile(self, input_file) -> Dict[str, list[float]]:
        result = dict()
        for base in ('A', 'C', 'G', 'T'):
            result[base] = list(map(float, input_file.readline().strip().split(" ")))
        return result

    def test_most_probable_kmers(self):

        for i in range(1, 4):
            with open(f"ProfileMostProbableKmer/inputs/input_{i}.txt") as input_file:
                genome = input_file.readline().strip()
                k = int(input_file.readline().strip())
                profile = self.read_profile(input_file)


            result = motif.find_profile_most_probable_kmer(genome, profile, k)

            with open(f"ProfileMostProbableKmer/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip()

            print(f"profile most probable kmers are {result} compared to expected {expected}")
            self.assertEqual(expected, result[0])


    def test_big_most_probable_kmers(self):

        with open("dataset_30305_3.txt") as input_file:
            genome = input_file.readline().strip()
            k = int(input_file.readline().strip())
            profile = self.read_profile(input_file)

        result = motif.find_profile_most_probable_kmer(genome, profile, k)

        print(f"profile most probable kmers are {result}")


    def test_greedy_motif_search(self):

        for i in range(5, 6):
            with open(f"GreedyMotifSearch/inputs/input_{i}.txt") as input_file:
                args = input_file.readline().strip().split(" ")
                k = int(args[0])
                genomes = input_file.readline().strip().split(" ")

            result = motif.greedy_motif_search(genomes, k)

            with open(f"GreedyMotifSearch/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip().split(" ")

            print(f"profile most probable kmers are {result} compared to expected {expected}")
            self.assertEqual(expected, result)


    def test_big_greedy_motif_search(self):

        with open("dataset_30305_5.txt") as input_file:
            args = input_file.readline().strip().split(" ")
            k = int(args[0])
            genomes = input_file.readline().strip().split(" ")

        result = motif.greedy_motif_search(genomes, k)

        string = " ".join(result)
        print(f"profile most probable kmers are {string}")


    def test_greedy_motif_search_psuedocounts(self):

        for i in range(1, 2):
            with open(f"GreedyMotifSearchPseudocounts/inputs/input_{i}.txt") as input_file:
                args = input_file.readline().strip().split(" ")
                k = int(args[0])
                genomes = input_file.readline().strip().split(" ")

            result = motif.greedy_motif_search_with_pseudocounts(genomes, k)

            with open(f"GreedyMotifSearchPseudocounts/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip().split(" ")

            print(f"profile most probable kmers are {result} compared to expected {expected}")
            self.assertEqual(expected, result)


    def test_big_greedy_motif_search_pseudocounts(self):

        with open("dataset_30306_9.txt") as input_file:
            args = input_file.readline().strip().split(" ")
            k = int(args[0])
            genomes = input_file.readline().strip().split(" ")

        result = motif.greedy_motif_search_with_pseudocounts(genomes, k)

        string = " ".join(result)
        print(f"profile most probable kmers are {string}")


    def test_randomized_motif_search(self):

        for i in range(1, 4):
            with open(f"RandomizedMotifSearch/inputs/input_{i}.txt") as input_file:
                args = input_file.readline().strip().split(" ")
                k = int(args[0])
                genomes = input_file.readline().strip().split(" ")

            best_result = []
            best_score = 1000000
            for _ in range(0, 1000):
                result = motif.randomized_motif_search(genomes, k)
                score = motif.score_motifs(result)
                if score < best_score:
                    best_result = result
                    best_score = score


            with open(f"RandomizedMotifSearch/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip().split(" ")

            print(f"best score {best_score} achieved for {best_result}")
            print(f"profile most probable kmers are {best_result} compared to expected {expected}")
            print(f"Expected score {motif.score_motifs(expected)}")
            self.assertEqual(expected, best_result)


    def test_big_randomized_motif_search(self):

        with open("dataset_30307_5.txt") as input_file:
            args = input_file.readline().strip().split(" ")
            k = int(args[0])
            genomes = input_file.readline().strip().split(" ")


        best_result = []
        best_score = 100000000
        for _ in range(0, 1000):
            result = motif.randomized_motif_search(genomes, k)
            score = motif.score_motifs(result)
            if score < best_score:
                best_result = result
                best_score = score

        print(f"best score {best_score} achieved for {best_result}")
        print("profile most probable kmers are " + " ".join(best_result))

    def test_gibbs_search(self):

        with open("GibbsSampler/inputs/input_1.txt") as input_file:
            args = input_file.readline().strip().split(" ")
            k = int(args[0])
            t = int(args[1])
            N = int(args[2])
            genomes = []
            for line in input_file:
                for genome in line.strip().split(" "):
                    genomes.append(genome)

        best_result = []
        best_score = 100000000
        for _ in range(0, 20):
            result = motif.gibbs_search(genomes, k, N)
            score = motif.score_motifs(result)
            if score < best_score:
                best_result = result
                best_score = score

        with open(f"GibbsSampler/outputs/output_1.txt") as output_file:
            expected = output_file.readline().strip().split(" ")

        print(f"best score {best_score} achieved for {best_result}")
        print(f"profile most probable kmers are {best_result} compared to expected {expected}")
        print(f"Expected score {motif.score_motifs(expected)}")
        self.assertEqual(expected, best_result)
