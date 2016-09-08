# Functions for importing vectors and matrix from the input files for ec_ssa

def import_vector(file,A):
    rubbish = file.readline()
    i = 0
    j = 0
    k = 0
    while rubbish[i] != '\n':
        if rubbish[i] == ' ':
            if A.typecode == 'l' or A.typecode == 'i':
                A[k] = int(rubbish[j:i])
            elif A.typecode == 'f' or A.typecode == 'd':
                A[k] = float(rubbish[j:i])
            j = i
            k = k + 1
        i = i + 1
        if rubbish[i] == '\n':
            if A.typecode == 'l' or A.typecode == 'i':
                A[k] = int(rubbish[j:i])
            elif A.typecode == 'f' or A.typecode == 'd':
                A[k] = float(rubbish[j:i])

def import_matrix(file,A,c,r):
    for l in range(0,r):
        rubbish = file.readline()
        i = 0
        j = 0
        k = 0
        while rubbish[i] != '\n':
            if rubbish[i] == ' ':
                if A.typecode == 'l' or A.typecode == 'i':
                    A[l*c+k] = int(rubbish[j:i])
                elif A.typecode == 'f' or A.typecode == 'd':
                    A[l*c+k] = float(rubbish[j:i])
                j = i
                k = k + 1
            i = i + 1
            if rubbish[i] == '\n':
                if A.typecode == 'l' or A.typecode == 'i':
                    A[l*c+k] = int(rubbish[j:i])
                elif A.typecode == 'f' or A.typecode == 'd':
                    A[l*c+k] = float(rubbish[j:i])
