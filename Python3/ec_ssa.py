#!/bin/python3
from aux_fun import *
from array import *
from random import *
from math import *
from ctypes import *

#### Header ####

print(" ")
print(" Greetings.")
print(" This is an SSA for electrochemical systems, developed and implemented")
print(" by Beruski, O. This script follows the description given in:")
print("  http://arxiv.org/abs/1608.07507")
print(" ")

#### Input and declaration of variables and parameters ####

# t = Simulation time and time resolution, 2x1 float
# T = Temperature, float
# S = Chosen setup, char
#     s = P for potentiostatic, S = G for galvanostatic
# E = Applied electrode potential, relevant if S = P, float
# I = Applied electric current, relevant if S = G, float
# N = Number of species, integer
# M = Number of reaction channels, integer
# X = Population number of each species, Nx1 long integer
# C = Specific transition coefficients of each channel, Mx1 float
# R = Reaction stoichiometric coefficients, MxN integer
# P = Product stoichiometric coefficients, MxN integer
# n = Electron stoichiometric coefficients, Mx1 integer
#     n[i] < 0 for reduction, n[i] > 0 for oxidation, n[i] = 0 for non-EC
# A = Transfer coefficient for each channel, Mx1 float
# E0 = Formal electrode potential for each channel, Mx1 float

# Despite giving me the chills, most variables are declared as the code
# proceeds, since N and M are dinamically allocated. Please bear with me.

print(" Input parameters are strongly recommended to be provided by a file.")
inp_choice = input(" Is the input given by a file? (y/n): ")

if inp_choice == "y":
    print(" The input will be done automatically.")
    print(" ")
    inp = input(" Please type the name of the input file: ")
    fin = open(inp,'r')
    t = array('d',range(0,2))
    import_vector(fin,t)
    fin.readline()
    T = float(fin.readline())
    fin.readline()
    S = fin.read(1)
    fin.read(1)
    if S == 'P':
        E = float(fin.readline())
    elif S == 'G':
        I = float(fin.readline())
    elif S == 'N':
        fin.readline()
    else:
        print(" ")
        print(" Problem with input, line 5. Please define the setup:")
        print("   -P = potentiostatic")
        print("   -G = galvanostatic")
        print("   -N = non-electrochemical")
    fin.readline()
    N = int(fin.readline())
    M = int(fin.readline())
    X = array('l',range(0,N))
    C = array('d',range(0,M))
    R = array('i',range(0,N*M))
    P = array('i',range(0,N*M))
    if S != 'N':
        n = array('i',range(0,M))
        A = array('d',range(0,M))
        E0 = array('d',range(0,M))
    import_vector(fin,X)
    import_vector(fin,C)
    fin.readline()
    import_matrix(fin,R,N,M)
    fin.readline()
    import_matrix(fin,P,N,M)
    if S != 'N':
        fin.readline()
        import_vector(fin,n)
        import_vector(fin,A)
        import_vector(fin,E0)
    fin.close()
else:
    print(" The input will be done manually. You have been warned.")
    print(" ")
    print(" Please provide the following parameters.")
    t = array('d',range(0,2))
    t[0] = float(input(" -Simulation time: "))
    t[1] = float(input(" -Time resolution: "))
    T = float(input(" -Temperature: "))
    S = input(" -Chosen setup (P/G/N): ")
    if S == 'P':
        E = float(input("   -Electrode potential: "))
    elif S == 'G':
        I = float(input("   -Electric current: "))
    N = int(input(" -Number of species: "))
    M = int(input(" -Number of reactions: "))
    X = array('l',range(0,N))
    C = array('d',range(0,M))
    R = array('i',range(0,N*M))
    P = array('i',range(0,N*M))
    if S != 'N':
        n = array('i',range(0,M))
        A = array('d',range(0,M))
        E0 = array('d',range(0,M))
    print(" -Initial population:")
    for i in range(0,N):
        X[i] = int(input("   -Species "+str(i+1)+": "))
    print(" -Transition coefficients:")
    for i in range(0,M):
        C[i] = float(input("   -Channel "+str(i+1)+": "))
    print(" -Reaction coefficients:")
    for i in range(0,M):
        print("   -Channel "+str(i+1)+":")
        for j in range(0,N):
            R[i*N+j] = int(input("     -Species "+str(j+1)+": "))
    print(" -Production coefficients: ")
    for i in range(0,M):
        print("   -Channel "+str(i+1)+":")
        for j in range(0,N):
            P[i*N+j] = int(input("     -Species "+str(j+1)+": "))
    if S != 'N':
        print(" -Electron coefficients:")
        for i in range(0,M):
            n[i] = int(print("   -Channel "+str(i+1)+": "))
        print(" -Transfer coefficients:")
        for i in range(0,M):
            A[i] = float(print("   -Channel "+str(i+1)+": "))
        print(" -Formal electrode potentials:")
        for i in range(0,M):
            E0[i] = float(print("   -Channel "+str(i+1)+": "))
    print(" ")
    print(" For future reference, would you like to print an input file based")
    inp_choice = input(" on the information above? (y/n): ")
    if inp_choice == 'y':
        inp = input(" Please provide the name of the file: ")
        fin = open(inp,'w')
        fin.write("{0:.2e} {1:.2e}\n".format(t[0],t[1]))
        fin.write("\n")
        fin.write("{0:.4e}\n".format(T))
        fin.write("\n")
        fin.write(S+"\n")
        if S == 'N':
            fin.write("\n")
        elif S == 'P':
            fin.write("{0:.4e}\n".format(E))
        elif S == 'G':
            fin.write("{0:.4e}\n".format(I))
        fin.write("\n")
        fin.write("{}\n".format(N))
        fin.write("{}\n".format(M))
        for i in range(0,N):
            if i == N-1:
                fin.write("{}\n".format(X[i]))
            else:
                fin.write("{} ".format(X[i]))
        for i in range(0,M):
            if i == M-1:
                fin.write("{0:.4e}\n".format(C[i]))
            else:
                fin.write("{0:.4e} ".format(C[i]))
        fin.write("\n")
        for i in range(0,M):
            for j in range(0,N):
                if j == N-1:
                    fin.write("{}\n".format(R[i*N+j]))
                else:
                    fin.write("{} ".format(R[i*N+j]))
        fin.write("\n")
        for i in range(0,M):
            for j in range(0,N):
                if j == N-1:
                    fin.write("{}\n".format(P[i*N+j]))
                else:
                    fin.write("{} ".format(P[i*N+j]))
        if S == 'P' or S == 'G':
            for i in range(0,M):
                if i == M-1:
                    fin.write("{}\n".format(n[i]))
                else:
                    fin.write("{} ".format(n[i]))
            for i in range(0,M):
                if i == M-1:
                    fin.write("{0:.2e}\n".format(A[i]))
                else:
                    fin.write("{0:.2e} ".format(A[i]))
            for i in range(0,M):
                if i == M-1:
                    fin.write("{0:.3e}\n".format(E0[i]))
                else:
                    fin.write("{0:.3e} ".format(E0[i]))
        fin.close()
    else:
        print(" Truly a shame.")

print(" ")
print(" The output is also strongly recommended to be provided in a file.")
out_choice = input(" Is the output to be given in a file? (y/n): ")

if out_choice == 'y':
    print(" You chose... wisely.")
    out = input(" Please type the name of the output file: ")
    print(" ")
else:
    print(" You chose... poorly.")
    input()

q_e = 96485.33289/6.022E23 # Fundamental charge
f = 96845.33289/(T*8.3144621) # Reduced Faraday constant
C_THRESH = 1E-6 # Threshold for convergence of the minimization procedure
D_THRESH = 1E3 # Threshold for divergence of the electrode potential

s = 1 # Sample number
ns = int(t[0]/t[1])+1 # Number of samples
tt = 0 # Total time

tau = 0 # Time for next reaction
mu = 0 # Channel fired

a0 = 0 # Total reaction propensity
a = array('d',range(0,M)) # Individual reaction propensities

ST = array('d',range(0,ns)) # Sampled times
SX = array('l',range(0,N*ns)) # Sampled population numbers

if S == 'P':
    I = 0 # Electric current
    SI = array('d',range(0,ns)) # Sampled electric current
elif S == 'G':
    E = 0 # Electrode potential
    SE = array('d',range(0,ns)) # Sampled electrode potential

d = 0 # Dummy variables
x = 0
dx = 0
dE = 0

#### The Stochastic Simulation Algorithm ####

seed('None',2)

# Writing initial conditions
for i in range(0,N): 
    SX[i] = X[i]

print(" The SSA will now begin.")

if S == 'N':
    print(" Non-electrochemical setup.")
    print(" ")
    while s < ns:
        # Resetting variables
        a0 = 0
        d = 0
        for i in range(0,M):
            a[i] = C[i]
        # Calculating reaction propensities
        for i in range(0,M):
            for j in range(0,N):
                if R[i*N+j] == 1:
                    a[i] = a[i]*X[j]
                elif R[i*N+j] == 2:
                    a[i] = a[i]*X[j]*(X[j]-1)/2
                elif R[i*N+j] == 3:
                    a[i] = a[i]*X[j]*(X[j]-1)*(X[j]-2)/6
                elif R[i*N+j] > 3:
                    print(" -Reaction order not implemented for indeces:")
                    print("  R["+str(i)+","+str(j)+"] = "+str(R[i*N+j]))
                    print(" ")
            a0 = a0 + a[i]
        # Checking for equilibrium
        if a0 == 0:
            ST[s] = tt
            for i in range(0,N):
                SX[s*N+i] = X[i]
            print(" Zero total reaction propensity.")
            print(" Terminating SSA.")
            print(" ")
            break
        # Sampling the random variables (tau,mu)
        tau = (1/a0)*log(1/uniform(0,1))
        mu = a0*uniform(0,1)
        i = 0
        while d < mu:
            d = d + a[i]
            i = i + 1
        mu = i-1
        # Applying the Markov jump
        tt = tt + tau
        for i in range(0,N):
            X[i] = X[i] - R[mu*N+i] + P[mu*N+i]
        # Checking for sampling time
        if tt >= s*t[1]:
            ST[s] = tt
            for i in range(0,N):
                SX[s*N+i] = X[i]
            s = s + 1
elif S == 'P':
    print(" Potentiostatic setup, E = {}.".format(E))
    print(" ")
    while s < ns:
        # Resetting variables
        a0 = 0
        I = 0
        d = 0
        for i in range(0,M):
            a[i] = C[i]
        # Calculating reaction propensities
        for i in range(0,M):
            for j in range(0,N):
                if R[i*N+j] == 1:
                    a[i] = a[i]*X[j]
                elif R[i*N+j] == 2:
                    a[i] = a[i]*X[j]*(X[j]-1)/2
                elif R[i*N+j] == 3:
                    a[i] = a[i]*X[j]*(X[j]-1)*(X[j]-2)/6
                elif R[i*N+j] > 3:
                    print(" -Reaction order not implemented for indeces:")
                    print("  R["+str(i)+","+str(j)+"] = "+str(R[i*N+j]))
                    print(" ")
            a[i] = a[i]*exp(f*n[i]*A[i]*(E-E0[i]))
            a0 = a0 + a[i]
            I = I + n[i]*q_e*a[i]
        # Checking for equilibrium
        if a0 == 0:
            ST[s] = tt
            SI[s] = I
            for i in range(0,N):
                SX[s*N+i] = X[i]
            print(" Zero total reaction propensity.")
            print(" Terminating SSA.")
            print(" ")
            break
        # Sampling the random variables (tau,mu)
        tau = (1/a0)*log(1/uniform(0,1))
        mu = a0*uniform(0,1)
        i = 0
        while d < mu:
            d = d + a[i]
            i = i + 1
        mu = i-1
        # Applying the Markov jump
        tt = tt + tau
        for i in range(0,N):
            X[i] = X[i] - R[mu*N+i] + P[mu*N+i]
        # Checking for sampling time
        if tt >= s*t[1]:
            ST[s] = tt
            SI[s] = I
            for i in range(0,N):
                SX[s*N+i] = X[i]
            s = s + 1
elif S == 'G':
    print(" Galvanostatic setup, I = {}.".format(I))
    print(" ")
    while s < ns:
        # Resetting variables
        a0 = 0
        E = I/fabs(I)
        dE = 1.0
        d = 0
        for i in range(0,M):
            a[i] = C[i]
        # Estimating reaction propensities
        for i in range(0,M):
            for j in range(0,N):
                if R[i*N+j] == 1:
                    a[i] = a[i]*X[j]
                elif R[i*N+j] == 2:
                    a[i] = a[i]*X[j]*(X[j]-1)/2
                elif R[i*N+j] == 3:
                    a[i] = a[i]*X[j]*(X[j]-1)*(X[j]-2)/6
                elif R[i*N+j] > 3:
                    print(" -Reaction order not implemented for indeces:")
                    print("  R[{},{}] = {}".format(i,j,R[i*N+j]))
                    print(" ")
        # Calculating electrode potential - Newton-Raphson
        while fabs(dE) >= C_THRESH:
            x = -I/q_e
            dx = 0.0
            for i in range(0,M):
                x = x + n[i]*a[i]*exp(f*n[i]*A[i]*(E-E0[i]))
                dx = dx + f*n[i]*n[i]*A[i]*a[i]*exp(f*n[i]*A[i]*(E-E0[i]))
            dE = x/dx
            E = E - dE
            if fabs(E) > D_THRESH:
                break
        # Checking for divergence
        if fabs(E) > D_THRESH:
            print(" Divergence of electrode potential.")
            print(" Terminating SSA.")
            print(" ")
            break
        # Correcting reaction propensities
        for i in range(0,M):
            a[i] = a[i]*exp(f*n[i]*A[i]*(E-E0[i]))
            a0 = a0 + a[i]
        # Checking for equilibrium
        if a0 == 0:
            ST[s] = tt
            SE[s] = E
            for i in range(0,N):
                SX[s*N+i] = X[i]
            print(" Zero total reaction propensity.")
            print(" Terminating SSA.")
            print(" ")
            break
        # Sampling the random variables (tau,mu)
        tau = (1/a0)*log(1/uniform(0,1))
        mu = a0*uniform(0,1)
        i = 0
        while d < mu:
            d = d + a[i]
            i = i + 1
        mu = i-1
        # Applying the Markov jump
        tt = tt + tau
        for i in range(0,N):
            X[i] = X[i] - R[mu*N+i] + P[mu*N+i]
        # Checking for sampling time
        if tt >= s*t[1]:
            ST[s] = tt
            SE[s] = E
            for i in range(0,N):
                SX[s*N+i] = X[i]
            s = s + 1

#### Output printing and other shenanigans ####

if out_choice == 'y':
    fout = open(out,'w')
    fout.write("# Output generated by ec_ssa.py.\n")
    if inp_choice == 'y':
        fout.write("# Input parameters given in file "+inp+".\n")
    else:
        fout.write("# Input parameters given manually.\n")
    if S == 'N':
        fout.write("# Non-electrochemical run.\n")
    elif S == 'P':
        fout.write("# Potentiostatic run.\n")
    elif S == 'G':
        fout.write("# Galvanostatic run.\n")
    fout.write("#\n")
    fout.write("# Time ")
    if S == 'P':
        fout.write(" Current ")
    elif S == 'G':
        fout.write(" Potential ")
    for i in range(0,N):
        fout.write("  S{}  ".format(i+1))
    fout.write("\n")
    for i in range(0,s):
        fout.write("{0:.4e}".format(ST[i]))
        if S == 'P':
            fout.write(" {0:.4e}".format(SI[i]))
        elif S == 'G':
            fout.write(" {0:.4e}".format(SE[i]))
        for j in range(0,N):
            fout.write(" {0:6}".format(SX[i*N+j]))
        fout.write("\n")
    fout.close()
else:
    print("Time".center(8),end="  ")
    if S == 'P':
        print(" I".center(8),end=" ")
    elif S == 'G':
        print(" E".center(8),end=" ")
    for i in range(0,N):
        print("  S{}  ".format(i+1).center(6),end=" ")
    print("\n")
    for i in range(0,s):
        print("{0:.4e}".format(ST[i]),end=" ")
        if S == 'P':
            print(" {0:.4e}".format(SI[i]),end=" ")
        elif S == 'G':
            print(" {0:.4e}".format(SE[i]),end=" ")
        for j in range(0,N):
            print(" {0:6}".format(SX[i*N+j]),end=" ")
        print("\n")
    
print(" SSA has ended accordingly.")
print(" ")
