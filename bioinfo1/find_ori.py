



def find_ori_indices_helper(text: str, pattern: str) -> list[int]:

    indices = []
    for i in range(0, 1 + len(text) - len(pattern)):
        if text[i:i + len(pattern)] == pattern:
            indices.append(i)

    return indices


def find_ori_indices(input_file: str) -> list[int]:

    lines = []
    with open(input_file) as f:
        lines = f.readlines()

    if len(lines) != 2:
        raise Exception("Input file must contain exactly two lines")

    pattern = lines[0].strip()
    text = lines[1].strip()

    return find_ori_indices_helper(text, pattern)


def find_ori_count(input_file: str) -> int:

    lines = []
    with open(input_file) as f:
        lines = f.readlines()

    if len(lines) != 2:
        raise Exception("Input file must contain exactly two lines")

    text = lines[0].strip()
    pattern = lines[1].strip()

    return len(find_ori_indices_helper(text, pattern))


def find_most_frequent_helper(text: str, k: int) -> dict[str, int]:

    pattern_map = dict()

    for i in range(0, 1 + len(text) - k):
        pattern = text[i:i + k]
        value = pattern_map.get(pattern, 0)
        pattern_map[pattern] = value + 1

    return pattern_map


def find_most_frequent(input_file: str) -> list[str]:

    lines = []
    with open(input_file) as f:
        lines = f.readlines()

    if len(lines) != 2:
        raise Exception("Input file must contain exactly two lines")

    text = lines[0].strip()
    k = int(lines[1])

    return max_map(find_most_frequent_helper(text, k))


def threshold_map(pattern_map: dict[str, int], t: int) -> list[str]:
    val = []
    for k, v in pattern_map.items():
        if v >= t:
            val.append(k)

    return val


def max_map(pattern_map: dict[str, int]) -> list[str]:
    val = []
    max = 0
    for k, v in pattern_map.items():
        if v > max:
            max = v
            val = [k]
        elif v == max:
            val.append(k)

    return val


def reverse_complement(sequence: str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    result = ''
    for base in sequence[::-1]:
        result += complement[base]

    return result


def find_k_clumps(genome: str, k: int, l: int, t: int) -> set[str]:

    clump_kmers = set()
    for i in range(0, len(genome) - l + 1):
        if i % 100000 == 0:
            print(i)
        sequence = genome[i:i + l]
        frequency_map = find_most_frequent_helper(sequence, k)
        kmers = threshold_map(frequency_map, t)
        clump_kmers.update(kmers)

    return clump_kmers


def find_skew(genome: str) -> list[int]:

    lookup = { 'C': -1, 'G': 1, 'T': 0, 'A': 0}
    skew = [0]
    for base in genome:
        skew.append(skew[-1] + lookup[base])

    return skew


def find_minimum_skew_indices(skew: list[int]) -> list[int]:

    minimum = 1
    minimum_skew_indices = []
    for index, value in enumerate(skew):
        if value < minimum:
            minimum = value
            minimum_skew_indices = [index]
        elif value == minimum:
            minimum_skew_indices.append(index)

    return minimum_skew_indices


def find_minimum_skew_helper(genome: str) -> list[int]:

    skew = find_skew(genome)

    return find_minimum_skew_indices(skew)


def compute_hamming_distance(genome: str, other: str) -> int:
    hamming_distance = 0
    for g, o in zip(genome, other):
        if g != o:
            hamming_distance += 1
    return hamming_distance


def approximate_pattern_match(genome: str, pattern: str, mismatches: int) -> list[int]:

    indices = []
    for i in range(0, 1 + len(genome) - len(pattern)):
        if compute_hamming_distance(genome[i:i + len(pattern)], pattern) <= mismatches:
            indices.append(i)

    return indices

def find_neighbors(pattern: str, d: int) -> list[str]:

    neighbors = []
    bases = ('A', 'C', 'G', 'T')

    stack = [("", 0)]
    while len(stack) > 0:
        item = stack.pop()
        if len(item[0]) == len(pattern):
            neighbors.append(item[0])
        else:
            prefix = item[0]
            for base in bases:
                if base == pattern[len(prefix)]:
                    stack.append((prefix + base, item[1]))
                elif item[1] < d:
                    stack.append((prefix + base, item[1] + 1))

    return neighbors


def find_frequent_words_mismatches(genome: str, pattern_len: int, mismatches: int) -> list[str]:

    pattern_map = dict()

    for i in range(0, 1 + len(genome) - pattern_len):
        pattern = genome[i:i + pattern_len]
        neighbors = find_neighbors(pattern, mismatches)
        for neighbor in neighbors:
            value = pattern_map.get(neighbor, 0)
            pattern_map[neighbor] = value + 1

    return max_map(pattern_map)


def find_frequent_words_mismatches_reverse_complement(genome: str, pattern_len: int, mismatches: int) -> list[str]:

    pattern_map = dict()

    for i in range(0, 1 + len(genome) - pattern_len):
        pattern = genome[i:i + pattern_len]
        neighbors = find_neighbors(pattern, mismatches)
        for neighbor in neighbors:
            value = pattern_map.get(neighbor, 0)
            pattern_map[neighbor] = value + 1
            rc_neighbor = reverse_complement(neighbor)
            rc_value = pattern_map.get(rc_neighbor, 0)
            pattern_map[rc_neighbor] = rc_value + 1

    return max_map(pattern_map)


