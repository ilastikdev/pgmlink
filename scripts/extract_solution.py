import pickle

input_file = open('mitocheck_solutions.dump','rb')

original_solution = pickle.load(input_file)
cp_solution = pickle.load(input_file)
icm_solution = pickle.load(input_file)
block_icm_solution = pickle.load(input_file)

input_file.close()


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

def extract_differences(filename):
        f = open(filename, 'rt')
        found = False
        for line in f:
                if not found:
                        if "Differences" in line:
                                found = True
                else:
                        return line.strip().split(" ")
        return None

def count_diff(a,b):
    d=0           
    for (i,j) in zip(a,b):
        if i != j:
            d+=1
    return d

def count_diff_no_disagreement(original_solution, block_icm_solution, differences):
	c = 0
	for (a,b,d) in zip(original_solution, block_icm_solution, differences):
	    if a != b:
	        if d == '0':
	            c+=1
    return c
