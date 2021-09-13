import random

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


def CountWithPseudocounts(Motifs):
	count = {}  # initializing the count dictionary
	k = len(Motifs[0])
	for symbol in "ACGT":
		count[symbol] = []
		for j in range(k):
			count[symbol].append(1)
	t = len(Motifs)
	for i in range(t):
		for j in range(k):
			symbol = Motifs[i][j]
			count[symbol][j] += 1
	return count


def ProfileWithPseudocounts(Motifs):
	k = len(Motifs[0])
	count = CountWithPseudocounts(Motifs)
	profile = {}
	count_values = list(count.values())
	for key in count:
		profile[key] = []
		for i in range(k):
			profile[key].append(count[key][i]/(count_values[0][i] +
								count_values[1][i]+count_values[2][i]+count_values[3][i]))
	return profile


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
	BestMotifs = []
	for i in range(0, t):
		BestMotifs.append(Dna[i][0:k])
	n = len(Dna[0])
	for i in range(n-k+1):
		Motifs = []
		Motifs.append(Dna[0][i:i+k])
		for j in range(1, t):
			P = ProfileWithPseudocounts(Motifs[0:j])
			Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
		if Score(Motifs) < Score(BestMotifs):
			BestMotifs = Motifs
	return BestMotifs


def Motifs(Profile, Dna):
	l = len(Profile['A'])
	result = []
	for dna in Dna:
		most = dna[:l]
		most_p = Pr(most, Profile)
		for i in range(1, len(dna)-l+1):
			p = Pr(dna[i:i+l], Profile)
			if p > most_p:
				most_p = p
				most = dna[i:i+l]
		result.append(most)
	return result


def RandomMotifs(Dna, k, t):
	result = []
	for dna in Dna:
		rand = random.randint(0, len(dna)-k)
		result.append(dna[rand:rand+k])
	return result

def RandomizedMotifSearch(Dna, k, t):
	M = RandomMotifs(Dna, k, t)
	BestMotifs = M
	while True:
		Profile = ProfileWithPseudocounts(M)
		M = Motifs(Profile, Dna)
		if Score(M) < Score(BestMotifs):
			BestMotifs = M
		else:
			return BestMotifs
