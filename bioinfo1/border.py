# Homeland Security has identified likely crossing positions along a national border. These are given in discrete linear mile positions along the border from east to west.

# Homeland Security has also identified the optimal locations for Sentry Towers along the border, in the same format; linear positions along the border from east to west.

# Given the two lists (likely crossings and optimal Sentry Tower locations) find the minimum sensor range that, if applied to all Sentry Towers, would ensure all likely crossings are covered by Sentry Towers.

# Example 1:

# Border Crossings  = [1,2,3,4], Optimal Sentry Tower Locations  = [2]
# Output: 2

# Example 2:

# Border Crossings  = [1,3,5,7], Optimal Sentry Tower Locations  = [2,6]
# Output: 1

from typing import List
from bisect import bisect_left

from copy import deepcopy


def find_min_range(crossings: List, towers: List) -> int:
    max_distance = 0

    for crossing in crossings:
        index = bisect_left(towers, crossing)

        distance = 0
        if index > 0 and index < len(towers):
            distance = min(abs(crossing - towers[index]), abs(crossing - towers[index - 1]))
        if index == 0:
            distance = towers[index] - crossing
        elif index == len(towers):
            distance = crossing - towers[index - 1]

        if distance > max_distance:
            max_distance = distance

    return max_distance


crossings = [1, 3, 5, 7]
towers = [2, 6]

print(find_min_range(crossings, towers))
print(find_min_range([1, 2, 3, 4], [2]))
print(find_min_range([1, 2, 3, 4], [0]))
print(find_min_range([1, 2, 3, 4], [5]))
print(find_min_range(crossings, [1, 3, 5, 7]))
print(find_min_range([1, 10, 11, 21], [5, 13, 18]))


# Given the minimum sensor range you’ve determined, how would you return the list of essential tower locations? i.e. remove any towers that aren’t required to cover all likely crossings.

# Border Crossings  = [1, 5, 11], Optimal Sentry Tower Locations  = [2, 3, 6, 8]
# Output: Minimum Sensor Range: 3, Essential Towers: [2, 8] or [3, 8]

# Border Crossings = [3,5,7] Sentry Tower Locations = [1,5,7]

def remove_tower(crossings: List, towers: List, min_range: int) -> List:
    for i in range(len(towers)):
        towers_copy = deepcopy(towers)
        towers_copy.pop(i)
        potential_range = find_min_range(crossings, towers)
        if potential_range == min_range:
            return towers_copy


def find_min_towers(crossings: List, towers: List) -> List:
    min_range = find_min_range(crossings, towers)

    old_towers = towers
    new_towers = remove_tower(crossings, towers, min_range)
    while new_towers != old_towers:
        old_towers = new_towers
        new_towers = remove_tower(crossings, old_towers, min_range)

    return new_towers















