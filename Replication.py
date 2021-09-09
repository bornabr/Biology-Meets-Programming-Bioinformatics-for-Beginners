def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = PatternCount(Text, Pattern)
    return freq


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def Reverse(Pattern):
    rev = ""
    for c in Pattern:
        rev = c + rev
    return rev


def Complement(Pattern):
    com = ""
    for c in Pattern:
        if c == "A":
            com = com + "T"
        elif c == "C":
            com = com + "G"
        elif c == "G":
            com = com + "C"
        else:
            com = com + "A"
    return com


def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern)  # reverse all letters in a string
    Pattern = Complement(Pattern)  # complement each letter in a string
    return Pattern


def PatternMatching(Pattern, Genome):
    positions = []  # output variable
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions


def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(Genome[0:n//2], symbol)

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array


def SkewArray(Genome):
    Skew = []
    Skew.append(0)
    for i in range(1, len(Genome)+1):
        nucleotid = Genome[i-1]
        if nucleotid == 'G':
            Skew.append(Skew[i-1] + 1)
        elif nucleotid == 'C':
            Skew.append(Skew[i-1] - 1)
        else:
            Skew.append(Skew[i-1])
    return Skew


def MinimumSkew(Genome):
    positions = []  # output variable
    Skew = SkewArray(Genome)
    minValue = min(Skew)
    for i in range(len(Skew)):
        if Skew[i] == minValue:
            positions.append(i)
    return positions


def HammingDistance(p, q):
    count = 0
    for x, y in zip(p, q):
        if x != y:
            count += 1
    return count


def ApproximatePatternMatching(Text, Pattern, d):
    positions = []  # output variable
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions


def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count = count+1
    return count


print(HammingDistance("CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG",
      "ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT"))
