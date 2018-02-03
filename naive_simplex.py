# The naive simplex algorithm as presented in Ch. 2
#   of our textbook (Bland's rule)
#
# usage: python3 naive_simplex.py [-v] data_file_name
# Input file format: a sequence of numbers separated 
# by space,  tab or end-of-line, starting with the number 
# of variables n, the number of constrains m, then the 
# augmented matrix of the equation constraints [A b], 
# then the coefficients of the objective function that 
# is being maximized.
   
# Note that this code is just an illustration and has not 
# been optimized for speed, and I haven't learned a lot 
# of CS so the code quality is sloppy :)

import sys

verbose=0

def print_row(row):
    out_str=' '.join([str(d) for d in row])
    return print(out_str)

def print_mat(M):
    [print_row(r) for r in M]

#turn the i-th column of M into e_k
def elim(i, k, M):
    for j in range(len(M)):
        if j!=k:
            r=[x-y*M[j][i]/M[k][i] for (x, y) in zip(M[j], M[k])]
            M[j]=r
    r=[x/M[k][i] for x in M[k]]
    M[k]=r
    return M

def simplex_from_basics(basics, o_coeffs, M):
#initialize
    n=len(M[0])-1
    m=len(M)
    cur_o_coeffs=list(o_coeffs)
    for j in range(m):
        new_o_coeffs=[x-cur_o_coeffs[basics[j]]*y for (x,y)
                      in zip(cur_o_coeffs, M[j][:-1])]
        cur_o_coeffs=new_o_coeffs
#change basic variables
    while True:
        if verbose==1:
            print_row(basics)
            print_mat(M)
            print_row(cur_o_coeffs)
            print("\n")
        flag=0
        for i in range(n):
            if cur_o_coeffs[i]>0:
                flag=1
                break
        if flag==0:
            return 0
        flag=0
        for j in range(m):
            if M[j][i]>0:
                flag=1
                break
        if flag==0:
            return 1
        opt_ratio=M[j][n]/M[j][i]
        outv=j
        for j in range(m):
            if M[j][i]>0 and (opt_ratio>M[j][n]/M[j][i] or
                              (opt_ratio==M[j][n]/M[j][i] and
                               basics[outv]>basics[j])):
                outv=j
                opt_ratio=M[j][n]/M[j][i]
        elim(i,outv,M)
        basics[outv]=i
        new_o_coeffs=[x-cur_o_coeffs[i]*y for (x,y)
                      in zip(cur_o_coeffs, M[outv][:-1])]
        cur_o_coeffs=new_o_coeffs

# read data
if(sys.argv[1]=="-v"):
    verbose=1
data_file=open(sys.argv[-1],"r")
data_lst=iter(str.split(data_file.read()))
n=int(next(data_lst))
m=int(next(data_lst))
M=[]
for i in range(m):
    c_row=[]
    for j in range(n):
        c_row+=[float(next(data_lst))]
    for j in range(m):
        if j==i:
            c_row+=[1]
        else:
            c_row+=[0]
    c_row+=[float(next(data_lst))]
    M+=[c_row]
o_coeffs2=[]
for i in range(n):
    o_coeffs2+=[float(next(data_lst))]
print_mat(M)
print_row(o_coeffs2)
if verbose==1:
    print("Phase 1")
basics=list(range(n,m+n))
o_coeffs1=[0]*n+[-1]*m
simplex_from_basics(basics, o_coeffs1, M)
for b in basics:
    if b>=n:
        print("No feasible solution.")
        exit(0)
if verbose==1:
    print("Phase 2")
for i in range(len(M)):
    M[i]=M[i][:n]+[M[i][m+n]]
r=simplex_from_basics(basics, o_coeffs2, M)
if r==1:
    print("No optimal solution.")
    exit(0)
opt=0
optv=[0]*n
for i in range(m):
    optv[basics[i]]=M[i][n]
    opt+=o_coeffs2[basics[i]]*M[i][n]
print_row(optv)
print(opt)

