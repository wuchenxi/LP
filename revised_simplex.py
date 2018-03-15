# The naive revised simplex method using the original pivoting
# rule. Note that this can not deal with degeneration properly.
#
# usage: python3 revised_simplex.py [-v] data_file_name
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

def simplex_from_basics(basics, o_coeffs, M, Binv):
#initialize    
    n=len(M[0])-1
    m=len(M)
 #change basic variables
    while True:
        if verbose==1:
            print_row(basics)
            print_mat(Binv)
            print("\n")
        #Find entering variable
        maxcoeff=0
        enter_id=0
        cbBinv=[0]*m
        for i in range(m):
            for j in range(m):
                cbBinv[i]+=o_coeffs[basics[j]]*Binv[j][i]
        for i in range(n):
            if i not in basics:
                coeff=o_coeffs[i]
                for j in range(m):
                    coeff-=cbBinv[j]*M[j][i]
                if coeff>maxcoeff:
                    maxcoeff=coeff
                    enter_id=i
        if maxcoeff==0:
            return 0
        #Now find exiting variable
        exit_id=0
        Binvb=[0]*m
        BinvA=[0]*m
        for i in range(m):
            for j in range(m):
                BinvA[i]+=Binv[i][j]*M[j][enter_id]
                Binvb[i]+=Binv[i][j]*M[j][n]
        ratio=-1
        for i in range(m):
            if BinvA[i]>0 and (ratio==-1 or Binvb[i]/BinvA[i]<ratio):
                ratio=Binvb[i]/BinvA[i]
                exit_id=i                
        if ratio==-1:
            return 1
        basics[exit_id]=enter_id
        # update matrix
        Binv[exit_id]=[x/BinvA[exit_id] for x in Binv[exit_id]]
        for i in range(exit_id):
            Binv[i]=[y-x*BinvA[i] for (x, y) in zip(Binv[exit_id],Binv[i])]
        for i in range(exit_id+1,m):
            Binv[i]=[y-x*BinvA[i] for (x, y) in zip(Binv[exit_id],Binv[i])]

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
B=[[0.0]*m for i in range(m)]
for i in range(m):
    B[i][i]=1
o_coeffs1=[0]*n+[-1]*m
simplex_from_basics(basics, o_coeffs1, M, B)
for b in basics:
    if b>=n:
        print("No feasible solution.")
        exit(0)
if verbose==1:
    print("Phase 2")
for i in range(len(M)):
    M[i]=M[i][:n]+[M[i][m+n]]
r=simplex_from_basics(basics, o_coeffs2, M, B)
if r==1:
    print("No optimal solution.")
    exit(0)
opt=0
optv=[0]*n
for i in range(m):
    for j in range(m):
        optv[basics[i]]+=B[i][j]*M[j][n]
        opt+=o_coeffs2[basics[i]]*B[i][j]*M[j][n]
print_row(optv)
print(opt)
