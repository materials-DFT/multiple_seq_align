#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:21:09 2021

@author: kianaamaral

install lxml parser

"""
from bs4 import BeautifulSoup
import pandas as pd 
import os
import copy

# target protein pdb id (e.g. 1ai2)
targetpdb = ''

# prot: protein name (e.g. Adk)
prot = ''

# path: where alignment xml files are saved 
path = ''

# xmlList: list of pdb ID's of target's homologs
'''In order to run this script it is important that the first alignment in the xmlList has the longest protein sequence '''

xmlList = []



def get1Dsequence(xmlfile,path):
    '''Parses xml file and saves AA sequence in a list'''
    chain1seq = ''
    chain2seq = ''  
    infile = open(path + targetpdb + '_' +  xmlfile + '.xml')
    contents = infile.readlines()
    content = "".join(contents) 
    soup = BeautifulSoup(content, "lxml")
    soupString = soup.get_text()
    #print(soupString)
    for i in range(len(soupString)):
        if soupString[i:i+4] == "Note":
            alignment = soupString[:i]
            for i in range(len(alignment)):
                if alignment[i:i+5] == 'Chain':
                    alignment = alignment[i:]
                    break
    
    for i in range(len(alignment)):
        if alignment[i:i+7] == 'Chain 1':
            temp1 = alignment[i+11:]
            for i in range(len(temp1)):
                if temp1[i].isalpha() or temp1[i] == '-':
                    start = i
                    break
                
            chain1 = temp1[start:]
            for i in range(len(chain1)):
                if chain1[i] == '\n': 
                    end = i 
                    break
                
            chain1 = chain1[:end]
            chain1seq = chain1seq + chain1
        elif alignment[i:i+7] == 'Chain 2':
            temp2 = alignment[i+11:]
            for i in range(len(temp2)):
                if temp2[i].isalpha()or temp2[i] == '-' and temp2[i+1].isnumeric() == False:
                    start = i
                    break
                
            chain2 = temp2[start:]
            for i in range(len(chain2)):
                if chain2[i] == '\n': 
                    end = i 
                    break
            chain2 = chain2[:end]
            chain2seq = chain2seq + chain2  
            
    chain1seq = Convert(chain1seq)
    chain2seq = Convert(chain2seq)
    
    return chain1seq, chain2seq 
 
           
def matchstart(targ,rawlist,a):
    '''Aligns beginning of protein'''
    rawstart = 0
    for ni,i in enumerate(rawlist[a]):
        if i!='~':
            rawstart = ni
            break
    #print (rawstart,rawlist[a][0+rawstart:5+rawstart])
    for k in range(len(targ)):
        if rawlist[a][0+rawstart:5+rawstart]==targ[k:k+5]:
            break
    #print (rawlist[a][k:5+k],targ[k:k+5])
    while rawlist[a][k:k+5]!=targ[k:k+5]:
        rawlist[a].insert(0,'~')
        rawlist[a+1].insert(0,'~')
  
        
                        
def twoNewlistInsert(newlist,a,b):
    '''Takes two target that indicate insertions from the homolog and 
    aligns them using *, which directly modifies the aa list.
    Adding inserts into the homologs'''
    startfrom = -1
    newlistb = newlist[b]
    homologb = newlist[b+1]
    newlista = newlist[a]
    homologa = newlist[a+1]
    for nk, k in enumerate(newlistb):
        #if nk < 100: print(k,nk, startfrom)
        if k!=str(a):
            if newlistb[nk-1]==str(a):
                startfrom =startfrom - 1
                #print (nk,k,startfrom) 
            for nj, j in enumerate(newlista):
                if nj > startfrom:
                    if k==j:
                        startfrom = nj
                        break
                    elif k=='~' or j=='~':
                        startfrom = nj
                        break
                    elif k<"A" and j<"A":
                        startfrom = nj
                        break
                    elif k=="-": #insert into a
                        newlista.insert(nj, '-')
                        homologa.insert(nj,"-")
                        startfrom = nj
                        break
                    elif j=="-": #inserts into b
                        #print (j, nk, nj, k)
                        newlistb.insert(nk,str(a))
                        homologb.insert(nk,"-")
                        startfrom = nj
                        
                    
def comparelists(lists,rangest,rangeend):
    for i in range(rangest,rangeend):
        toprint = []
        for j in range(len(lists)):
            toprint.append(lists[j][i])
        print(i,toprint)
        
def listToString(s): 
    str1 = "" 
     
    return (str1.join(s))

def Convert(string):
    list1=[]
    list1[:0]=string
    return list1

        
###Run code here###
  
rawlist = [] 
for i in range(len(xmlList)):
    target, homolog = get1Dsequence(xmlList[i],path)
    rawlist.append(target)
    rawlist.append(homolog)
    

final = rawlist[:2] 
insertion = [] 

for i in range(len(rawlist)):
     if i%2 == 0 and i != 0: 
        for k in range(len(final[0])):
            if rawlist[i][0:5] == final[0][k:k+5]: 
                start = k 
                break  
                
        for k in range(start):
            rawlist[i].insert(0, '~')
            rawlist[i+1].insert(0, '-')
            
newlist = copy.deepcopy(rawlist)
lenraw = len(newlist)

for c in range(2,lenraw,2):
    rcurr = newlist[c]
    matchstart(newlist[0],newlist,c)
    print("running newlist",c)
    for d in range(c-2,-2,-2):
        prevr = newlist[d]
        twoNewlistInsert(newlist,d,c)
        
final = []
finalall = []

newlist[0].insert(0, targetpdb)
newlist[1].insert(0,xmlList[0])

for i in range(2, len(newlist),2):
   newlist[i].insert(0, xmlList[int(i/2)]) 
   newlist[i+1].insert(0,xmlList[int(i/2)])
   
final.append(newlist[0])
final.append(newlist[1])
finalall.append(newlist[0])
finalall.append(newlist[1])

for j in range(2,len(newlist),2):
    final.append(newlist[j+1])
    finalall.append(newlist[j])
    finalall.append(newlist[j+1])
    
df = pd.DataFrame(final)
filename = targetpdb + ' full multiple sequence alignment.xlsx'    
df.to_excel(filename, index=False)
       
for k in range(len(final)): 
    string = listToString(final[k][1:]) 
    final[k] = string   
    
textfile = open( prot +  ".fasta", "w")
textfile.write('>' + targetpdb + "\n")
textfile.write(final[0] + "\n")


for i in range(1,len(final)):                              
    textfile.write('>' + xmlList[i-1] + "\n")
    textfile.write(final[i] + "\n")
    
textfile.close()          
        
    
print("done")


        
