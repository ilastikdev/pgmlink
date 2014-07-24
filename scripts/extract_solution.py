def extract_solution(filename):
	f = open(filename, 'rt')
	found = False
	for line in f:
		if not found:
			if "Found Solution" in line or "Found solution" in line:
				found = True
		else:
			return line.strip().split(" ")
	return None
