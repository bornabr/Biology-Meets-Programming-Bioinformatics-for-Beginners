def Count(Motifs):
	count = {} # initializing the count dictionary
	k = len(Motifs[0])
	for symbol in "ACGT":
		count[symbol] = []
		for j in range(k):
			 count[symbol].append(0)
	t = len(Motifs)
	for i in range(t):
		for j in range(k):
			symbol = Motifs[i][j]
			count[symbol][j] += 1
	return count


def Profile(Motifs):
	t = len(Motifs)
	count = Count(Motifs)
	profile = {key: ([value/t for value in count[key]]) for key in count}
	return profile


def Consensus(Motifs):
	k = len(Motifs[0])
	count = Count(Motifs)
	consensus = ""
	for j in range(k):
		m = 0
		frequentSymbol = ""
		for symbol in "ACGT":
			if count[symbol][j] > m:
				m = count[symbol][j]
				frequentSymbol = symbol
		consensus += frequentSymbol
	return consensus


def Score(Motifs):
	consensus = Consensus(Motifs)
	k = len(Motifs[0])
	score = 0
	for motif in Motifs:
		for j in range(k):
			if motif[j] != consensus[j]:
				score += 1
	return score


def Pr(Text, Profile):
	l = len(Text)
	p = 1
	for i in range(l):
		p *= Profile[Text[i]][i]
		if p == 0:
			return p
	return p


def ProfileMostProbableKmer(text, k, profile):
	l = len(text)
	max_p = -1
	max_profile = None
	for i in range(l-k+1):
		tmp = Pr(text[i:i+k], profile)
		if tmp > max_p:
			max_p = tmp
			max_profile = text[i:i+k]
	return max_profile


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

Profile = {
	'A': [0.4,  0.3,  0.0,  0.1,  0.0,  0.9],
	'C': [0.2,  0.3,  0.0,  0.4,  0.0,  0.1],
	'G': [0.1,  0.3,  1.0,  0.1,  0.5,  0.0],
	'T': [0.3,  0.1,  0.0,  0.4,  0.5,  0.0]
}


print(Pr("CAGTGA", Profile))
