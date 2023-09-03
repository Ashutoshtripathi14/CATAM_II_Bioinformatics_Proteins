#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import numpy as np
import random


# In[47]:


#==========
#start of def

def s(Si,Tj):
    if(Si==Tj):
        return 0
    else:
        return 1

#end of def
#==========
#start of def

def edit_distance(S,T):
    m=len(S)
    n=len(T)
    D=[[0 for j in range(n+1)]for i in range(m+1)]
    D[0][0]= 0
    for i in range(1,m+1):
        D[i][0]=i
    for j in range(1,n+1):
        D[0][j]=j
    for i in range(1,m+1):
        for j in range(1,n+1):
            D[i][j]=min(min(D[i][j-1]+1,D[i-1][j]+1),s(S[i-1], T[j-1])+D[i-1][j-1])
    return D[m][n]
#Complexity of this is O(nm)

#end of def
#==========
#start of def

def convert_edit_sequence_to_alignment(S,T,pointers):
    m=len(S)
    n=len(T)
    i=m
    j=n
    str1=''
    str2=''
    str_match=''
    while(i>0 or j>0):
        if(pointers[i][j]=='D'):
            str1=S[i-1]+str1
            str2=T[j-1]+str2
            if(S[i-1]==T[j-1]):
                str_match='M'+str_match
            else:
                str_match='R'+str_match
            i-=1
            j-=1
        elif(pointers[i][j]=='L'):
            str1=' '+str1
            str2=T[j-1]+str2
            str_match='I'+str_match
            j-=1
        else:
            str1=S[i-1]+str1
            str2=' '+str2
            str_match='D'+str_match
            i-=1
    print(str_match)
    print(str1)
    print(str2)
    return str1, str2
        

#end of def 
#==========
#start of def

def convert_alignment_into_steps(str1, str2):
    pos=1
    step=1
    for i in range(len(str1)):
        if(step>50):
            break
        if(str1[i]==' '):
            print("\"Step "+ str(step)+ " : Insert "+ str(str2[i])+"\"")
            step+=1
            pos+=1
        elif(str2[i]==' '):
            print("\"Step "+ str(step)+ " : Delete "+ str(str1[i])+"\"")
            step+=1
        elif(str1[i]!=str2[i]):
            print("\"Step "+ str(step)+ " : Replace "+ str(str1[i])+" with "+str(str2[i])+"\"")
            step+=1
            pos+=1
        else:
            print("\"Step "+ str(step)+" : Match "+ str(str1[i])+"\"")
            step+=1
            pos+=1
        print()

#end of def 
#==========
#start of def

def edit_distance_with_pointer(S,T):
    m=len(S)
    n=len(T)
    D=[[0 for j in range(n+1)]for i in range(m+1)]
    pointers=[['' for j in range(n+1)] for i in range(m+1)]
    D[0][0]= 0
    for i in range(m+1):
        D[i][0]=i
        pointers[i][0]='U'
    for j in range(1,n+1):
        D[0][j]=j
        pointers[0][j]='L'
    for i in range(1,m+1):
        for j in range(1,n+1):
            D[i][j]=min(min(D[i][j-1]+1,D[i-1][j]+1),s(S[i-1], T[j-1])+D[i-1][j-1]) 
            if(i==1 and j==1):
                pointers[i][j]='D'
            elif(D[i][j]==D[i-1][j]+1):
                pointers[i][j]='U'
            elif(D[i][j]==D[i][j-1]+1):
                pointers[i][j]='L'
            else:
                pointers[i][j]='D'
    print(D[m][n])
    print("=====")
    str1, str2 = convert_edit_sequence_to_alignment(S,T,pointers)
    print("=====")
    convert_alignment_into_steps(str1, str2)
    


#end of def 
#==========
#start of def

def ParseProteinsTxt():
    '''
    This is a simple function to parse the Fasta files and return the amino acid sequence it reads.
    '''
    file_location='II-9-3-2022-proteins.txt'
    dict_of_sequences={}
    sequence_curr=""
    curr_sequence_name=""
    with open(file_location, "r") as fh:
        for line in fh:
            if(line.startswith("Protein")):
                curr_line=''.join(list(line))
                curr_sequence_name=curr_line[-2]
                dict_of_sequences[curr_sequence_name]=""
            elif(line.startswith("#")):
                if(sequence_curr!=""):
                    dict_of_sequences[curr_sequence_name]=sequence_curr
                sequence_curr=""
                curr_sequence_name=""
            else:
                seq_on_line = ''.join(list(line))
                seq_on_line=seq_on_line[:-1]
                sequence_curr+=seq_on_line
    return dict_of_sequences

#end of def 
#==========
#start of def

def ParseBLOSUM():
    file_location='II-9-3-2022-blosum.txt'
    blosum_matrix={}
    with open(file_location, "r") as fh:
        count=0
        aa_list=""
        for line in fh:
                if(count==0):
                    seq_on_line = ''.join(list(line.split(' ')))
                    seq_on_line=seq_on_line[:-1]
                    aa_list=seq_on_line
                else:
                    seq_on_line = ''.join(list(line))
                    seq_on_line=seq_on_line[:-1]
                    if(len(seq_on_line)>0):
                        nums_temp_=line[1:-1].split(' ')
                        nums_temp=[]
                        for curr_rem in nums_temp_:
                            if(len(curr_rem)!=0):
                                nums_temp.append(curr_rem)
                        nums=[int(n) for n in nums_temp]
                        curr_AA_row=seq_on_line[0]
                        for i in range(len(aa_list)): 
                            blosum_matrix[(curr_AA_row,aa_list[i])]=nums[i]
                count+=1
    return blosum_matrix

#end of def
#==========
#start of def

def new_s(Si,Tj):
    return blosum_matrix[(Si,Tj)]

#end of def
#==========
#start of def

def max_score_with_scoring_matrix(S,T):
    m=len(S)
    n=len(T)
    D=[['' for j in range(n+1)]for i in range(m+1)]
    pointers=[['' for j in range(n+1)] for i in range(m+1)]
    D[0][0]= 0
    for i in range(m+1):
        D[i][0]=i*-8
        pointers[i][0]='U'
    for j in range(1,n+1):
        D[0][j]=j*-8
        pointers[0][j]='L'
    for i in range(1,m+1):
        for j in range(1,n+1):
            D[i][j]=max(max(D[i][j-1]-8,D[i-1][j]-8),new_s(S[i-1], T[j-1])+D[i-1][j-1]) 
            if(i==1 and j==1):
                pointers[i][j]='D'
            elif(D[i][j]==D[i-1][j]-8):
                pointers[i][j]='U'
            elif(D[i][j]==D[i][j-1]-8):
                pointers[i][j]='L'
            else:
                pointers[i][j]='D'
    print(D[m][n])
    print("=====")
    str1, str2 = convert_edit_sequence_to_alignment(S,T,pointers)
    print("=====")
    convert_alignment_into_steps(str1, str2)

#Complexity of this is O(nm)   


#end of def
#==========
#start of def

def scoring_for_gaps(S,T,u, flag_for_printing):
    '''
    V(i,j)=max(E(i,j), F(i,j), G(i,j))
    E(i,j)=max_{0<=k<j}(V(i,k) + w(j-k))
    F(i,j)=max_{0<=k<i}(V(k,j) + w(i-k))
    G(i,j)=V(i-1,j-1)+s(Si,Tj)

    normally: O(mn+ n+m)
    w(l)={0; l=0 \\ u; l>=1}
    u<0
    E(i,j)=max_{0<=k<=j-1}(V(i,k)+w(j-k))=max(E(i,j-1), V(i,j-1)+w(1)) (only insertions)
    F(i,j)=u+max_{0<=k<=i-1}(V(k,j)+w(i-k))=max(F(i-1,j), V(i-1,j)+w(1)) (only deletions)
    G(i,j)=V(i-1,j-1)+s(Si,Tj) (replacement on last character)
    V(i,j)=max(E(i,j), F(i,j), G(i,j))
    V(0,0)=0
    E(0,0)=0
    F(0,0)=0
    V(i,0)=u(i>0)
    V(0,j)=u(j>0)
    E(i,0)=u (i>0)
    F(0,j)=u(j>0)
    E(0,j)=0 (j>0)
    F(i,0)=0 (i>0)
    '''
    m=len(S)
    n=len(T)
    V=[[0 for j in range(n+1)] for i in range(m+1)]
    E=[[0 for j in range(n+1)]for i in range(m+1)]
    F=[[0 for j in range(n+1)]for i in range(m+1)]
    G=[[0 for j in range(n+1)]for i in range(m+1)]
    pointers=[['' for j in range(n+1)] for i in range(m+1)]
    for i in range(m+1):
        V[i][0]=u
        E[i][0]=2*u
        pointers[i][0]='U'
    for j in range(n+1):
        V[0][j]=u
        F[0][j]=2*u
        pointers[0][j]='L'
    V[0][0]=0
    E[0][0]=0
    F[0][0]=0
    for i in range(1,m+1):
        for j in range(1,n+1):
            E[i][j]=max(E[i][j-1], V[i][j-1]+u)
            F[i][j]=max(F[i-1][j], V[i-1][j]+u)
            G[i][j]=V[i-1][j-1]+new_s(S[i-1],T[j-1])
            V[i][j]=max(E[i][j], F[i][j], G[i][j])
            if(V[i][j]==G[i][j]):
                pointers[i][j]='D'
            elif(V[i][j]==E[i][j]):
                pointers[i][j]='L'
            else:
                pointers[i][j]='U'
    if(flag_for_printing):
        print(V[m][n])
        print("=====")
        str1, str2 = convert_edit_sequence_to_alignment(S,T,pointers)
        print("=====")
        convert_alignment_into_steps(str1, str2)
    return V[m][n]
#end of def
#==========
#start of def

def new_new_s(Si,Tj):
    if(Si==Tj):
        return 1
    else:
        return -1

#end of def
#==========
#start of def

def scoring_for_gaps_new(S,T,u, flag_for_printing):
    m=len(S)
    n=len(T)
    V=[[0 for j in range(n+1)] for i in range(m+1)]
    E=[[0 for j in range(n+1)]for i in range(m+1)]
    F=[[0 for j in range(n+1)]for i in range(m+1)]
    G=[[0 for j in range(n+1)]for i in range(m+1)]
    for i in range(m+1):
        V[i][0]=u
        E[i][0]=2*u
    for j in range(n+1):
        V[0][j]=u
        F[0][j]=2*u
    V[0][0]=0
    E[0][0]=0
    F[0][0]=0    
    for i in range(1,m+1):
        for j in range(1,n+1):
            E[i][j]=max(E[i][j-1], V[i][j-1]+u)
            F[i][j]=max(F[i-1][j], V[i-1][j]+u)
            G[i][j]=V[i-1][j-1]+new_new_s(S[i-1],T[j-1])
            V[i][j]=max(E[i][j], F[i][j], G[i][j])
    return V[m][n]
#end of def
#==========
#start of def

def montecarloestimation(N, n, u, p_):
    Ui=['' for x in range(N)]
    Vi=['' for x in range(N)]
    V_gap=[0 for x in range(N)]
    rng = np.random.default_rng()
    for i in range(N):
        pred=rng.random()
        U=rng.choice(['a','b'],n, p=[p_,1-p_])
        Ui[i]=''.join(U)
        V=rng.choice(['a','b'],n, p=[p_,1-p_])
        Vi[i]=''.join(V)
        V_gap[i]=scoring_for_gaps_new(Ui[i], Vi[i], u, False)
    mean=sum(V_gap)/N
    return (mean/n)

#end of def
#==========
#start of def

def scoring_for_alignment(S,T,u_del, u_ins):
    '''
    V_s(i,j)=V_suf(i,j)=max(0,V_s(i-1,j)+s(Si,-), Vs(i,j-1)+s(-,Tj), Vs(i-1,j-1)+s(Si,Tj))
    s(Si,-)=u_del
    s(-,Tj)=u_ins
    s being used is new_s from Section 2 of the project
    v_sub(S,T)=max_i,j {V_s(i,j)}
    Assume u_ins and u_del are constants (ie s(Si,-) and s(-,Tj) are constants for all Si,Tj)
    '''
    m=len(S)
    n=len(T)
    V_s=[[0 for j in range(n+1)] for i in range(m+1)]
    vsub=0
    for i in range(m+1):
        V_s[i][0]=0
    for j in range(n+1):
        V_s[0][j]=0
    V_s[0][0]=0
    for i in range(1,m+1):
        for j in range(1,n+1):
            V_s[i][j]=max(0, V_s[i][j-1]+u_ins, V_s[i-1][j]+u_del, V_s[i-1][j-1]+new_s(S[i-1],T[j-1]))
            vsub=max(vsub,V_s[i][j])
    return vsub

#end of def
#==========


# In[ ]:


-2+(-2)+12+18=26 


# In[42]:


scoring_for_gaps('SSSSCCCS','SSSSCCC',-2, True)
#DMMMRMMR


# In[40]:


# rng = np.random.default_rng()
# p_=1/2
# U=rng.choice(['a','b'],10,p=[p_,1-p_])
# print(U)


# In[ ]:


for n in range(100,5000):
    print(montecarloestimation(10, n, -3, 1/2))


# In[ ]:


x=[]


# In[24]:


print(montecarloestimation(3, 5000, -3, 1/2))


# In[32]:


print(montecarloestimation(4, 2000, -3, 1/2))


# In[31]:


print(montecarloestimation(500, 10, -3, 1/2))


# In[50]:


print(montecarloestimation(500, 10, -3, 1/2))


# In[51]:


print(montecarloestimation(500, 100, -3, 1/2))


# In[52]:


print(montecarloestimation(500, 1000, -3, 1/2))


# In[35]:


0.37112-0.092


# In[38]:


0.27912*(5/4)+0.092


# In[2]:


0.42913799999999996-0.37112


# In[3]:


0.05801799999999996/0.27912000000000003


# In[4]:


0.09+0.28*(1/(4/5))


# In[26]:


print(montecarloestimation(500, 10000, -3, 1/2))


# In[5]:


blosum_matrix=ParseBLOSUM()
proteins=ParseProteinsTxt()


# In[5]:


edit_distance("shesells","seashells")


# In[43]:


edit_distance_with_pointer("shesells","seashells")


# In[15]:


wow="MRRRMRDDDDMMMRMMRMMRMRRRRMRRMMRMMMMRMMMMRDMIMRMMRMRRRRMRRRIMDMRRRDMMRMMMRRMMRRMMRMMRMRMRMMMMMRMMMMMMMMRMRIMIMDDMMMRRRRMRRRMRRRMRMRDDMRMRRRRRRRRRRMRRMRMMMRMMRM"
print(len(wow))
print(wow[:50])


# In[17]:


wow="MGLSDGEWQLVLKVWGKVEGDLPGHGQEVLIRLFKTHPETLEK FDKFKGLKTEDEMK ASADLKKHGGTVLTALGNILKKKGQHEAELKPLAQSHATKHKISIK F LEYISEAIIHVLQSKHSADFGADAQAAMGKALELFRNDMAAKYKEFGFQG"
print(wow[:50])


# In[44]:


edit_distance_with_pointer(proteins['A'],proteins['B'])


# In[ ]:





# In[18]:


wow="MDDRMDDRRRMMMRMMRMMRMRRRRMRRMMRMMMMRMMMMRRRMRMMRMRRRRMDRRRRMRRRRMMRMMMRRMMRRMMRMMRMRMRMMMMMRMMMMMMMMRMRRRRRMMMRRRRMRRRMRRRMRMRDDMRMRRRRRRRRRRMRRMRMMMRMMRM"
print(wow[:50])


# In[45]:


max_score_with_scoring_matrix(proteins['A'],proteins['B'])


# In[8]:


max_score_with_scoring_matrix(proteins['A'],proteins['B'])


# In[9]:


scoring_for_gaps('AMMA','MAMA',-12, True)


# In[19]:


wow="MDDDDDDDDDDDDDDDDDDDDDMRRRMRRRRRMRRDRMRRRRMRRMMRDDDDDRDRDMRMMDDDDMRMMRMRRMRMMRMRRRMRMRRRMMMMMMMMMMMMMMMMMMMMRMMRRMMMMMRRMRMRRRMRRMRRMMRMMRMMRMMMRMMMMMMMMMRMMRMMRMMRMMRMMMRMMMRRMMRRMMMMMMMRMRMMRRMRRMMRMMMMMMMMMMMRMMMMMMRMMMMMRMMMRMMMRMMRMMMMMMRMMRMRMMMRMMMMMMRMRMRMRMRRRRRMMMMRMRMRRRMMMMMMRMRRMMRMMMMRRMRMMRMMRMMRRMMMMMMMMMRRMRRMMMRMRRMRMMRMMMMMRMRMRMMRMRMRMMMMMMRMMMMMMMMMMMMMMMRMMMMRMMMMMRRMMRMMMRRMIIIIIIIIIRMMIIIIMMIIRIIMIIIIIIIIMIMRRRRRMRMMRRRDDDDMDDDDR"
wow[:50]


# In[46]:


scoring_for_gaps(proteins['C'],proteins['D'],-12, True)


# In[7]:


scoring_for_gaps(proteins['C'],proteins['D'],-12, True)


# In[10]:


scoring_for_gaps(proteins['C'],proteins['D'],-12, True)


# In[25]:


scoring_for_alignment(proteins['C'],proteins['D'],-2, -2)


# In[ ]:




