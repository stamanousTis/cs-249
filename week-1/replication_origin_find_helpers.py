# coding: utf-8

# Pattern Count (check for overlapping)
# Complexity is O(n-k+1)
def pattern_count(text, pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i:(i+len(pattern))] == pattern:
              count=count+1
    return count


# generate an array Count, where Count(i) stores Count(Text, Pattern) for Pattern = Text(i, k)
# Complexity O(n^2*k)

def naive_frequent_words(Text, k):

    Count=[]
    frequent_patterns=[]
    for i in range(0,len(Text)-k+1):
        Count.append(pattern_count(Text, Text[i:i+k]))
        #print(Text[i:i+k] + ' , ' + str(pattern_count(Text, Text[i:i+k])))

    maxCount = max(Count)
    for i in range(0,len(Text)-k+1):
        if Count[i] == maxCount:
            frequent_patterns.append(Text[i:i+k])

    return list(set(frequent_patterns))

# Use some memoisation to avoid reqeat same calculaations
# Input: Text,k
# Output: dictionary of pattern:occurences
# Complexity: O()
def frequency_table(Text, k):
    freqMap = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        if Pattern in freqMap:
            freqMap.update({Pattern:freqMap[Pattern]+1})
        else:
            freqMap.update({Pattern: 1})

    return freqMap


# MaxMap function
# Complexity O(n)
def map_max(_dict):
    max_value = next(iter(_dict.values()))
    for val in _dict.values():
        if val > max_value:
            max_value = val
    return max_value


def better_frequent_words(Text, k):
    frequent_patterns = []
    freq_map = frequency_table(Text, k)
    max_val = map_max(freq_map)
    for pattern in freq_map.keys():
        if freq_map[pattern] == max_val:
            frequent_patterns.append(pattern)
    
    return frequent_patterns



def reverse_complement(Text):
    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    Text = list(Text)
    complementary_seq = [base_complement[base] for base in Text] 
    return ''.join(complementary_seq[::-1])


def pattern_matching(pattern, genome):
    starting_positions = []
    for i in range(len(genome) - len(pattern) + 1):
        if genome[i:i+len(pattern)] == pattern:
            starting_positions.append(i) 
    return starting_positions


def find_clumps(Text, k, L, t):
    clumps = []
    for i in range(0, len(Text)-L+1):
        Window = Text[i:i+L]
        freqMap = frequency_table(Window, k)
        for pattern,occurences in freqMap.items():
            if occurences >= t:
                clumps.append(pattern)
    clump_list = list(set(clumps))
    return clump_list



def skew(Text):
    skewness = 0
    skew_lst = []
    for i in range(0, len(Text)):
        if i == 0:
            skew_lst.append(0);
        if (Text[i] == "G"):
            skewness = skewness + 1
        elif (Text[i]=='C'):
            skewness = skewness - 1
        else:
            skewness = skewness + 0

        skew_lst.append(skewness)
    return skew_lst

def skew_min_index(Text):
    skew_vec = skew(Text)
    skew_min = min(skew_vec)
    skew_min_index = []
    for i in range(len(skew_vec)):
        if skew_vec[i] == skew_min:
            skew_min_index.append(i)
        else:
            continue
    return skew_min_index


def approximate_pattern_count(Text, pattern, d):
    count = 0
    for i in range(len(Text) - len(pattern) + 1):
        if hamming_distance(Text[i:i+len(pattern)], pattern) <= d:
            count = count + 1
        else:
            continue            
    return count


def immediate_neighbors(pattern):
    nucleotide = ['A','C','G','T']
    neighborhood = []
    original_pattern = pattern
    for i in range(len(pattern)):
        pattern = original_pattern
        symbol = pattern[i]
        for j in range(len(nucleotide)):
            subst = nucleotide[j]
            if symbol != subst:
                pattern = list(pattern)
                pattern[i] = subst 
                pattern = ''.join(pattern)
                neighborhood.append(pattern)
    neighborhood.append(original_pattern)
    return list(set(neighborhood))


def neighbors(pattern, d):
    nucleotides_lst = ['A', 'C', 'G', 'T']
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return nucleotides_lst
    neighborhood = []
    suffix_neighbors = neighbors(pattern[1:],d)
    for suffix_neighbor in suffix_neighbors:
        if hamming_distance(pattern[1:],suffix_neighbor) < d:
            for nucleotide in nucleotides_lst:
                neighbor = nucleotide+suffix_neighbor
                neighborhood.append(neighbor)
        else: 
            neighbor = pattern[0]+suffix_neighbor
            neighborhood.append(neighbor)
    
    return neighborhood


def iterative_neighbors(pattern, d):
    neighborhood = list(pattern)
    for j in range(1,d):
        for neighbor in neighborhood:
            neighborhood.append(immediate_neighbors(neighbor))
            neighborhood = list(set(neighborhood))
    return neighborhood


def frequent_words_with_mismatches(Text, k, d):
    patterns = []
    freq_map = {}
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i+k]
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            if neighbor in freq_map:
                freq_map.update({neighbor: freq_map[neighbor] + 1})
            else:
                freq_map.update({neighbor : 1})
                
    max_occurs = map_max(freq_map)
    for pattern in freq_map.keys():
        if freq_map[pattern] == max_occurs:
            patterns.append(pattern)
    
    return patterns


def frequent_words_with_mismatches_rcs(Text, k , d):
    freq_map = {}
    patterns = []
    for i in range(len(Text)- k+1):
        pattern = Text[i:i+k]
        neighborhood = neighbors(pattern,d)
        for neighbor in neighborhood:
            if neighbor in freq_map:
                freq_map.update({neighbor: freq_map[neighbor] + 1})
            else:
                freq_map.update({neighbor : 1})

        neighborhood = neighbors(reverse_complement(pattern),d)
        for neighbor in neighborhood:
            if neighbor in freq_map:
                freq_map.update({neighbor: freq_map[neighbor] + 1})
            else:
                freq_map.update({neighbor : 1})

    max_occurs = map_max(freq_map)
    for pattern in freq_map.keys():
        if freq_map[pattern] == max_occurs:
            patterns.append(pattern)
            
    return patterns
