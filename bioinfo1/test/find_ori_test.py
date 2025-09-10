import unittest

import find_ori


class FindOriTest(unittest.TestCase):

    def test_find_ori_count(self):

        for i in range(1, 7):
            count = find_ori.find_ori_count(f"PatternCount/inputs/input_{i}.txt")
            with open(f"PatternCount/outputs/output_{i}.txt") as output_file:
                expected = int(output_file.read())
            print(f"Testing file {i}, with {count} pattern instances detected compared to {expected} expected.")
            self.assertEqual(count, expected)


    def test_big_dataset(self):

        count = find_ori.find_ori_count("./dataset_30272_6.txt")
        print(f"Testing big file with {count} pattern instances detected")


    def test_find_most_frequent(self):

        for i in range(1, 5):
            frequent_words = find_ori.find_most_frequent(f"FrequentWords/inputs/input_{i}.txt")
            with open(f"FrequentWords/outputs/output_{i}.txt") as output_file:
                expected = output_file.read().strip()
            frequent_words = " ".join(frequent_words)
            print(f"Testing file {i}, with {frequent_words} pattern instances detected compared to {expected} expected.")
            self.assertEqual(frequent_words, expected)


    def test_big_frequent_words(self):

        frequent_words = find_ori.find_most_frequent("./dataset_30272_13.txt")
        frequent_words = " ".join(frequent_words)
        print(f"Testing big file with {frequent_words} pattern instances detected")


    def test_reverse_complement(self):

        for i in range(1, 3):
            with open(f"ReverseComplement/inputs/input_{i}.txt") as input_file:
                text = input_file.read().strip()
            reverse_complement = find_ori.reverse_complement(text)
            with open(f"ReverseComplement/outputs/output_{i}.txt") as output_file:
                expected = output_file.read().strip()
            print(f"Testing file {i}, with {reverse_complement} compared to {expected} expected.")
            self.assertEqual(reverse_complement, expected)


    def test_big_reverse_complement(self):

        with open("./dataset_30273_2.txt") as input_file:
            text = input_file.read().strip()
        reverse_complement = find_ori.reverse_complement(text)
        print(f"Testing big file with {reverse_complement}")


    def test_find_ori_indices(self):

        for i in range(1, 6):
            indices_list = find_ori.find_ori_indices(f"PatternMatching/inputs/input_{i}.txt")
            indices = " ".join(map(str, indices_list))
            with open(f"PatternMatching/outputs/output_{i}.txt") as output_file:
                expected = output_file.read().strip()
            print(f"Testing file {i}, {indices_list}, with {indices} pattern instances detected compared to {expected} expected.")
            self.assertEqual(indices, expected)


    def test_big_find_ori_indices(self):

        indices_list = find_ori.find_ori_indices("./Vibrio_cholerae.txt")
        indices = " ".join(map(str, indices_list))
        print(f"Testing big file with the following pattern instances detected")
        print(indices)


    def test_find_clumps(self):

        for i in range(1, 5):
            with open(f"ClumpFinding/inputs/input_{i}.txt") as input_file:
                text = input_file.readline().strip()
                args = input_file.readline().strip().split(" ")
                z = int(args[0])
                l = int(args[1])
                t = int(args[2])

            result = find_ori.find_k_clumps(text, z, l, t)

            with open(f"ClumpFinding/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip()

            expected = set(expected.split(" "))
            self.assertEqual(result, expected)


    def test_big_find_clumps(self):

        with open(f"dataset_30274_5.txt") as input_file:
            text = input_file.readline().strip()
            args = input_file.readline().strip().split(" ")
            z = int(args[0])
            l = int(args[1])
            t = int(args[2])

        result = find_ori.find_k_clumps(text, z, l, t)

        print(result)


    def test_e_coli(self):

        with open(f"E_coli.txt") as input_file:
            text = input_file.readline().strip()
            args = input_file.readline().strip().split(" ")
            z = int(args[0])
            l = int(args[1])
            t = int(args[2])

        result = find_ori.find_k_clumps(text, z, l, t)

        print(len(result))


    def test_find_minimum_skew(self):

        for i in range(1, 6):
            with open(f"MinimumSkew/inputs/input_{i}.txt") as input_file:
                text = input_file.readline().strip()
            min_skews = find_ori.find_minimum_skew_helper(text)

            with open(f"MinimumSkew/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip()

            result = " ".join(map(str, min_skews))
            print(result)
            self.assertEqual(result, expected)


    def test_skew_test(self):

        genome = "GATACACTTCCCGAGTAGGTACTG"

        print(find_ori.find_minimum_skew_helper(genome))


    def test_approx_pattern_match_test(self):

        matches = find_ori.approximate_pattern_match("CATGCCATTCGCATTGTCCCAGTGA", "CCC", 2)
        print(len(matches))


    def test_big_find_minimum_skew(self):

        with open(f"dataset_30277_10.txt") as input_file:
            genome = input_file.readline().strip()

        min_skews = find_ori.find_minimum_skew_helper(genome)
        result = " ".join(map(str, min_skews))
        print(result)

    def test_hamming_distance(self):
        genome1 = "GGGCCGTTGGT"
        genome2 = "GGACCGTTGAC"
        genome3 = "GGACCGTCGAC"
        genome4 = "GGGCCGTTGGT"

        self.assertEqual(find_ori.compute_hamming_distance(genome1, genome2), 3)
        self.assertEqual(find_ori.compute_hamming_distance(genome1, genome3), 4)
        self.assertEqual(find_ori.compute_hamming_distance(genome1, genome4), 0)


    def test_big_hamming_distance(self):

        with open(f"dataset_30278_3.txt") as input_file:
            genome1 = input_file.readline().strip()
            genome2 = input_file.readline().strip()

        distance = find_ori.compute_hamming_distance(genome1, genome2)
        print(distance)


    def test_test_hamming_distance(self):

        genome1 = "TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC"
        genome2 = "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"

        distance = find_ori.compute_hamming_distance(genome1, genome2)
        print(distance)


    def test_approximate_pattern_match(self):

        for i in range(1, 8):
            with open(f"ApproximatePatternMatching/inputs/input_{i}.txt") as input_file:
                pattern = input_file.readline().strip()
                genome = input_file.readline().strip()
                mismatches = int(input_file.readline().strip())
            matches = find_ori.approximate_pattern_match(genome=genome, pattern=pattern, mismatches=mismatches)

            with open(f"ApproximatePatternMatching/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip()

            result = " ".join(map(str, matches))
            print(result)
            self.assertEqual(result, expected)


    def test_big_approximate_pattern_match(self):
        with open("dataset_30278_6.txt") as input_file:
            pattern = input_file.readline().strip()
            genome = input_file.readline().strip()
            mismatches = int(input_file.readline().strip())
        matches = find_ori.approximate_pattern_match(genome=genome, pattern=pattern, mismatches=mismatches)
        result = " ".join(map(str, matches))
        print(result)
        print(len(matches))


    def test_count_approximate_pattern_match(self):
        pattern = "AAAAA"
        genome = "AACAAGCTGATAAACATTTAAAGAG"
        mismatches = 2
        matches = find_ori.approximate_pattern_match(genome=genome, pattern=pattern, mismatches=mismatches)
        result = " ".join(map(str, matches))
        print(result)
        print(len(matches))


    def test_find_neighbors(self):

        for i in range(1, 5):
            with open(f"Neighborhood/inputs/input_{i}.txt") as input_file:
                pattern = input_file.readline().strip()
                mismatches = int(input_file.readline().strip())
            neighbors = find_ori.find_neighbors(pattern=pattern, d=mismatches)

            with open(f"Neighborhood/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip()

            result = set(neighbors)
            expected = set(expected.split(" "))
            print(result)
            self.assertEqual(result, expected)


    def test_neighbors_test(self):

        print(len(find_ori.find_neighbors(pattern="ACGT", d=3)))


    def test_big_neighbors(self):
        with open("dataset_30282_4.txt") as input_file:
            pattern = input_file.readline().strip()
            mismatches = int(input_file.readline().strip())
        neighbors = find_ori.find_neighbors(pattern=pattern, d=mismatches)
        print(" ".join(neighbors))


    def test_frequent_words_mismatches(self):

        for i in range(1, 5):
            with open(f"FrequentWordsMismatches/inputs/input_{i}.txt") as input_file:
                genome = input_file.readline().strip()
                params = input_file.readline().strip().split(" ")
                pattern_len = int(params[0])
                mismatches = int(params[1])
            frequent_words = find_ori.find_frequent_words_mismatches(genome=genome, pattern_len=pattern_len, mismatches=mismatches)

            with open(f"FrequentWordsMismatches/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip()

            result = set(frequent_words)
            expected = set(expected.split(" "))
            print(result)
            self.assertEqual(result, expected)


    def test_big_words_mismatches(self):
        with open("dataset_30278_9.txt") as input_file:
            genome = input_file.readline().strip()
            params = input_file.readline().strip().split(" ")
            pattern_len = int(params[0])
            mismatches = int(params[1])
        frequent_words = find_ori.find_frequent_words_mismatches(genome=genome, pattern_len=pattern_len,
                                                                 mismatches=mismatches)
        print(" ".join(frequent_words))


    def test_frequent_words_mismatches_reverse_complement(self):

        for i in range(1, 5):
            with open(f"FrequentWordsMismatchesReverseComplements/inputs/input_{i}.txt") as input_file:
                genome = input_file.readline().strip()
                params = input_file.readline().strip().split(" ")
                pattern_len = int(params[0])
                mismatches = int(params[1])
            frequent_words = find_ori.find_frequent_words_mismatches_reverse_complement(genome=genome, pattern_len=pattern_len, mismatches=mismatches)

            with open(f"FrequentWordsMismatchesReverseComplements/outputs/output_{i}.txt") as output_file:
                expected = output_file.readline().strip()

            result = set(frequent_words)
            expected = set(expected.split(" "))
            print(result)
            self.assertEqual(result, expected)


    def test_big_words_mismatches_reverse_complement(self):
        with open("dataset_30278_10.txt") as input_file:
            genome = input_file.readline().strip()
            params = input_file.readline().strip().split(" ")
            pattern_len = int(params[0])
            mismatches = int(params[1])
        frequent_words = find_ori.find_frequent_words_mismatches_reverse_complement(genome=genome,
                                                                                    pattern_len=pattern_len,
                                                                                    mismatches=mismatches)
        print(" ".join(frequent_words))


    def test_salmonella(self):

        with open("Salmonella_enterica.txt") as input_file:
            input_file.readline()
            genome = input_file.readline().strip()

        skew_mins = find_ori.find_minimum_skew_helper(genome)

        for skew_min in skew_mins:
            print(f"Examining genome at {skew_min}")

            segment = genome[skew_min:skew_min + 500]
            pattern_len = 9
            mismatches = 3
            frequent_words = find_ori.find_frequent_words_mismatches_reverse_complement(genome=segment,
                                                                                        pattern_len=pattern_len,
                                                                                        mismatches=mismatches)
            print(f"Found {len(frequent_words)} frequent words")
            print(" ".join(frequent_words))

