import math
from random import randint
from typing import Dict

import find_ori


def find_motif_enumerations(genomes: list[str], k: int, d: int) -> set[str]:

    motif_candidates = set()
    result = set()
    first = genomes[0]
    for i in range(0, len(first) - k + 1):
        motif_candidates.update(find_ori.find_neighbors(pattern=first[i:i+k], d=d))
    for motif in motif_candidates:
        candidate = True
        for genome in genomes[1:]:
            if not find_ori.approximate_pattern_match(genome, motif, d):
                candidate = False
                break
        if candidate:
            result.add(motif)
    return result


def compute_entropy(motifs: list[str]) -> float:

    length = len(motifs[0])

    entropy = 0
    for i in range(0, length):
        bases = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }
        for motif in motifs:
            bases[motif[i]] += 1
        for base in bases:
            if bases[base]:
                bases[base] /= len(motifs)
                entropy -= bases[base] * math.log(bases[base], 2)
    return entropy


def find_distances_between_pattern_and_strings(pattern: str, genomes: list[str]) -> int:

    distance = 0
    for genome in genomes:
        min_distance = len(pattern)
        for i in range(0, len(genome) - len(pattern) + 1):
            hamming = find_ori.compute_hamming_distance(pattern, genome[i:i+len(pattern)])
            if hamming < min_distance:
                min_distance = hamming
        distance += min_distance

    return distance


def find_median_string(genomes: list[str], k: int) -> list[str]:

    closest_pattern = 'A' * k
    min_distance = len(genomes[0]) * k
    neighbors = find_ori.find_neighbors(pattern=closest_pattern, d=k)
    medians = [closest_pattern]

    for pattern in neighbors:
        distance = find_distances_between_pattern_and_strings(pattern, genomes)
        if distance < min_distance:
            min_distance = distance
            medians = [pattern]
        elif distance == min_distance:
            medians.append(pattern)

    return medians


def score_kmer_probability(kmer: str, profile: Dict[str, list[float]]) -> float:

    score = 1.0
    for i in range(0, len(kmer)):
        score *= profile[kmer[i]][i]

    return score


def find_profile_most_probable_kmer(genome: str, profile: Dict[str, list[float]], k: int) -> list[str]:

    high_score = -1
    most_probable_kmers = []

    for i in range(0, len(genome) - k + 1):
        kmer = genome[i:i+k]
        score = score_kmer_probability(kmer, profile)
        if score > high_score:
            high_score = score
            most_probable_kmers = [kmer]
        elif score == high_score:
            most_probable_kmers.append(kmer)

    # print(f"Most probable kmer: {most_probable_kmers} score: {high_score}")
    return most_probable_kmers


def compute_profile(kmers: list[str]) -> Dict[str, list[float]]:

    profile = {
        'A': [ 0.0 ] * len(kmers[0]),
        'C': [ 0.0 ] * len(kmers[0]),
        'G': [ 0.0 ] * len(kmers[0]),
        'T': [ 0.0 ] * len(kmers[0]),
    }
    for i in range(0, len(kmers[0])):
        for kmer in kmers:
            base = kmer[i]
            profile[base][i] += 1
        for base in profile:
            profile[base][i] /= len(kmers)

    return profile


def compute_profile_with_pseudocounts(kmers: list[str]) -> Dict[str, list[float]]:

    profile = {
        'A': [ 1.0 ] * len(kmers[0]),
        'C': [ 1.0 ] * len(kmers[0]),
        'G': [ 1.0 ] * len(kmers[0]),
        'T': [ 1.0 ] * len(kmers[0]),
    }
    for i in range(0, len(kmers[0])):
        for kmer in kmers:
            base = kmer[i]
            profile[base][i] += 1
        for base in profile:
            profile[base][i] /= len(kmers) + 4

    return profile


def score_motifs(motifs: list[str]) -> int:
    score = 0
    for i in range(0, len(motifs[0])):
        bases = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }
        for motif in motifs:
            bases[motif[i]] += 1

        score += len(motifs) - max(bases.values())

    return score



def greedy_motif_search(genomes: list[str], k: int) -> list[str]:

    best_motifs = [ genome[0:k] for genome in genomes ]
    profile = compute_profile(best_motifs)
    best_score = score_motifs(best_motifs)

    for i in range(0, len(genomes[0]) - k + 1):
        motifs = [ genomes[0][i:i+k] ]
        print(f"Evaluating motif {motifs[0]}")
        for genome in genomes[1:]:
            profile = compute_profile(motifs)
            motif = find_profile_most_probable_kmer(genome, profile, k)[0]
            motifs.append(motif)
        score = score_motifs(motifs)
        if score > best_score:
            best_score = score
            best_motifs = motifs

    return best_motifs


def greedy_motif_search_with_pseudocounts(genomes: list[str], k: int) -> list[str]:

    best_motifs = [ genome[0:k] for genome in genomes ]
    profile = compute_profile_with_pseudocounts(best_motifs)
    best_score = score_motifs(best_motifs)

    for i in range(0, len(genomes[0]) - k + 1):
        motifs = [ genomes[0][i:i+k] ]
        print(f"Evaluating motif {motifs[0]}")
        for genome in genomes[1:]:
            profile = compute_profile_with_pseudocounts(motifs)
            motif = find_profile_most_probable_kmer(genome, profile, k)[0]
            motifs.append(motif)
        score = score_motifs(motifs)
        if score > best_score:
            best_score = score
            best_motifs = motifs

    return best_motifs


def find_most_probable_kmers(genomes: list[str], profile: Dict[str, list[float]], k: int) -> list[str]:

    return [ find_profile_most_probable_kmer(genome, profile, k)[0] for genome in genomes ]


def print_profile(profile: Dict[str, list[float]]) -> None:

    for base in profile:
        probs = profile[base]
        print(f"{base}: " + " ".join([ f"{prob:.2f}" for prob in probs ]))


def randomized_motif_search(genomes: list[str], k: int) -> list[str]:

    best_motifs = []
    for genome in genomes:
        index = randint(0, len(genome) - k)
        best_motifs.append(genome[index:index + k])
    best_profile = compute_profile_with_pseudocounts(best_motifs)
    count = 0
    # print("Initial motifs: " + " ".join(best_motifs))
    # print(f"Score of initial motifs {score_motifs(best_motifs)}")
    # print_profile(best_profile)
    while True:
        # print(f"Iteration {count}")
        count += 1
        motifs = find_most_probable_kmers(genomes, best_profile, k)
        profile = compute_profile_with_pseudocounts(motifs)
        score = score_motifs(motifs)
        # print("New motifs: " + " ".join(motifs))
        # print(f"Score of newly generated motifs {score}")
        # print_profile(profile)
        if score < score_motifs(best_motifs):
            best_motifs = motifs
            best_profile = profile
        else:
            return best_motifs


def generate_random(p: list[float]) -> int:

    num = randint(0, 1000000000)
    total = sum(p)
    normalized_p = [ 1000000000 * n / total for n in p ]

    pdf = 0
    for i, n in enumerate(normalized_p):
        pdf += n
        if pdf > num:
            return i


def find_profile_random_kmer(genome: str, profile: Dict[str, list[float]], k: int) -> str:

    kmer_counts = []

    for i in range(0, len(genome) - k + 1):
        kmer = genome[i:i+k]
        score = score_kmer_probability(kmer, profile)
        kmer_counts.append(score)

    # print("Computed kmer probabilities: " + " ".join([ f"{k}" for k in kmer_counts ]))

    random_index = generate_random(kmer_counts)
    kmer = genome[random_index:random_index+k]
    # print(f"Returning kmer {kmer} at index {random_index} of genome {genome}")

    return kmer


def gibbs_search(genomes: list[str], k: int, iterations: int) -> list[str]:

    motifs = []
    for genome in genomes:
        index = randint(0, len(genome) - k)
        motifs.append(genome[index:index + k])
    best_score = score_motifs(motifs)
    best_motifs = motifs
    print("Initial motifs: " + " ".join(motifs))
    print(f"Score of initial motifs {best_score}")
    for i in range(iterations):
        # print(f"Iteration {count}")
        random_index = randint(0, len(genomes) - 1)
        motifs.pop(random_index)
        # print(f"Iteration {i}: random index: {random_index}")
        # print("motifs: " + " ".join(motifs))
        profile = compute_profile_with_pseudocounts(motifs)

        motif = find_profile_random_kmer(genomes[random_index], profile, k)
        motifs.insert(random_index, motif)
        # print(f"New motif for index {random_index}: {motif}")
        # print("motifs: " + " ".join(motifs))
        score = score_motifs(motifs)
        if score < best_score:
            best_score = score
            best_motifs = motifs
        # print("New motifs: " + " ".join(motifs))
        # print(f"Score of newly generated motifs {score}")
        # print_profile(profile)

    print(f"Motifs after {iterations} iterations: " + " ".join(motifs))
    print(f"Score of final motifs {score_motifs(motifs)}")

    return best_motifs
