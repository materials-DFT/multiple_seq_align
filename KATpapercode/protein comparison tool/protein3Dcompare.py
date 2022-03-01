"""
Created on Sun Sep 27 14:40:55 2020

@author: kianaamaral
"""
import requests
from bs4 import BeautifulSoup
import pandas as pd
from pymol import cmd
import os

print('Enter pdbDirectory (This is where the input pdb files are saved): ')
pdbDirectory = input()
print('Enter outputDirectory (This is where the output .cif files will be saved): ')
outputDirectory = input() 


class protein3Dcompare:
    
    def __init__(self):
        pass
            
    def testClass(self,mystring):
        print(mystring)
        
        
####### UNIPROT DATA #######      
        
    def webscrape_gene(self,uniprot_ids):
        LIST_OF_GENES = [] 
        for items in uniprot_ids:
            # WEB SCRAPES GENE 
            URL = 'https://www.uniprot.org/uniprot/' + items
            page = requests.get(URL)
            soup = BeautifulSoup(page.content, 'html.parser')
            results = soup.find_all('a')
            for result in results:
                if 'eco' in result.text: 
                    s = result.text
                    gene = s.strip('eco: ')
                    LIST_OF_GENES.append(gene)
                    break
                
        return LIST_OF_GENES
    
    def get_geneDict(self,uniprot_ids): 
        gene_dict = {}
        for items in uniprot_ids:
            # WEB SCRAPES GENE 
            URL = 'https://www.uniprot.org/uniprot/' + items
            page = requests.get(URL)
            soup = BeautifulSoup(page.content, 'html.parser')
            results = soup.find_all('a')
            for result in results:
                if 'eco' in result.text: 
                    s = result.text
                    gene = s.strip('eco: ')
                    gene_dict[items] = gene
                    break
        
        return gene_dict
         
    def pdbidToUniprot(self, pdbid): 
        URL = 'https://www.rcsb.org/structure/' + pdbid
        page = requests.get(URL)
        soup = BeautifulSoup(page.content, 'html.parser')
        uniprot = 'none'
        results = soup.find_all('a')
        for result in results:
            u = result.text
            if len(u) == 6 and u[1].isdigit():
                uniprot = u
                break
        
            
        return uniprot 
    
    def getOrganism(self, uniprot):
        if uniprot == 'none': 
            organism = 'none'
        else: 
            URL = 'https://www.uniprot.org/uniprot/' + uniprot 
            page = requests.get(URL)
            soup = BeautifulSoup(page.content, 'html.parser')
            result = soup.find_all('em')
            organism = str(result)
            organism = organism[5:]
            organism = organism[:-6]
            
        return organism
    
    def getName(self, pdbid):
        URL = 'https://www.rcsb.org/structure/' + pdbid
        page = requests.get(URL)
        soup = BeautifulSoup(page.content, 'html.parser')
        result = soup.find_all("span", id= "structureTitle")  
        new_string = str(result)
        index = new_string.find('>')
        new_string = new_string[index+1:]
        index = new_string.find('<')
        name = new_string[:index]

        return name
    
    def getLigand(self, pdbid):
        URL = 'https://www.rcsb.org/structure/' + pdbid
        page = requests.get(URL)
        soup = BeautifulSoup(page.content, 'html.parser')
        results = soup.find_all("a") 
        ligands = ''
        for result in results: 
            if 'Query' in result.text:
                string = result.text
                index = string.find('on')
                new_string = string[index+2:]
                ligands = ligands + new_string + ','
        ligands = ligands[:-1]

        return ligands 


                
    def getTaxonomy(self, uniprot):
        if uniprot == 'none': 
            taxonomyList = 'none'
        else:
            taxonomyList = []
            hiddenTaxonlist = []
            URL = 'https://www.uniprot.org/uniprot/' + uniprot
            page = requests.get(URL)
            soup = BeautifulSoup(page.content, 'html.parser')
            results = soup.find(id = 'names_and_taxonomy')
            Newresults= results.find_all('span')
            for result in Newresults:
                if 'hiddenTaxon' in str(result):
                    resultFound = False
                    final = result.text
                    while resultFound == False:
                        if final[-1].isalpha():
                            final = final
                            resultFound = True
                        else:
                            final = final[:-1] 
                    hiddenTaxonlist.append(final)
                    
            Newerresults= results.find_all('a') 
            for result in Newerresults:
                if 'taxonomy' in str(result):
                    taxonomyList.append(result.text)
            taxonomyList = taxonomyList[2:]
            for item in taxonomyList:
                if item in hiddenTaxonlist:
                    taxonomyList.remove(item)
                    
        return taxonomyList
    
    def getPubmed(self, pdbid):
        try:
            URL = 'https://www.rcsb.org/structure/' + pdbid
            page = requests.get(URL)
            soup = BeautifulSoup(page.content, 'html.parser')
            results = soup.find(id = 'pubmedLinks')
            Newresults= results.find_all('a')
            for result in Newresults:
                if 'querySearchLink' in str(result):
                    pubmedID = result.text       
                    PubmedLink = 'https://pubmed.ncbi.nlm.nih.gov/' + pubmedID +'/'
        except AttributeError:
            PubmedLink = ''
        
        return PubmedLink
    
    def getSequence(self, uniprot):
        if uniprot == 'none': 
            seq = 'none'
        else:
            URL = 'https://www.uniprot.org/uniprot/' + uniprot
            page = requests.get(URL)
            soup = BeautifulSoup(page.content, 'html.parser')
            results = soup.find(id = 'entrySequence')
            Newresults= results.find_all('pre')
            for result in Newresults:
                seq = result.text 
                
        return seq
    
    def getLength(self, pdbid):
        URL = 'https://www.rcsb.org/structure/' + pdbid
        page = requests.get(URL)
        soup = BeautifulSoup(page.content, 'html.parser')
        results = soup.find(id="MacromoleculeTable")
        Newresults= results.find_all('td')
        for result in Newresults:
            if result.text.isdigit():
                length = result.text
            
        return length 
    
    def webscrapeData(self, pdbid):
        featureslist = ['PDB', 'Uniprot', 'Organism', 'Name', 'Ligands', 'Global Stoichiometry', 'Sequence','Sequence Length', 'PubMed link', 'Taxonomy'] 
        pdb = pdbid
        uni = self.pdbidToUniprot(pdbid)
        org = self.getOrganism(uni)
        n = self.getName(pdbid)
        l = self.getLigand(pdbid)
        g = self.get_global_stoich(pdbid)
        s = self.getSequence(uni)
        sl = self.getLength(pdbid)
        p = self.getPubmed(pdbid)
        t = self.getTaxonomy(uni)
        
        df1 = pd.DataFrame([[pdb,uni,org,n,l,g,s,sl,p,t]], columns = featureslist )
        
        return df1
            
        
    def getHomologData(self, genename, targetprotein, pdbidlist):
        df = pd.DataFrame()
        targetdf = self.webscrapeData(targetprotein)
        df = df.append(targetdf)
        
        for pdbid in pdbidlist: 
            df1 = self.webscrapeData(pdbid)
    
            df = df.append(df1)
            
        filename = genename +'data.xlsx'    
        df.to_excel(filename, index=False)    
        return df 
    
####### FATCAT METHODS #######
        
    def importGlob(self, pdb1, pdb2):
        
        # Running FATCAT alignment, saving .cif file in output directory  
        os.system("bash runFATCAT.sh -pdb1 " + pdb1 + " -pdb2 " + pdb2 + " -pdbFilePath " + pdbDirectory + " -autoFetch -flexible -outputPDB -outFile " + outputDirectory + "pdbtemporary.pdb")
                                                                                           
                                
        proteinz = outputDirectory + 'pdbtemporary.pdb'
        
        cmd.load(proteinz, "proteinz")
        cmd.split_states("proteinz")
        cmd.delete("proteinz")
        cmd.set_name("proteinz_0001", "target")
        cmd.set_name("proteinz_0002", "homolog")
        cmd.color("red", "target")
        cmd.save(str(outputDirectory) + str(pdb1)+ "_" + str(pdb2) + ".cif")
        
        temppdb = outputDirectory + 'pdbtemporary.pdb'
        os.remove(temppdb)
        print("deleted temporary pdb")
        print("done")

    def shortFatcat(self, pdb1, pdb2):

        os.system("bash runFATCAT.sh -pdb1 " + pdb1+ " -pdb2 " + pdb2 + " -pdbFilePath " + pdbDirectory + " -autoFetch -flexible -printFatCat -show3d")
    
        
    def specificSingleChainAlignment(self, pdb1, pdb2, chain1, chain2):
        
        # Running FATCAT alignment, saving .cif file in output directory  
        os.system("bash runFATCAT.sh -pdb1 " + pdb1 + "." + chain1 + " -pdb2 " + pdb2 + "." + chain2 + " -pdbFilePath " + pdbDirectory + " -autoFetch -flexible -outputPDB -outFile " + outputDirectory + "pdbtemporary.pdb")
      
                                                                                                                
        proteinz = outputDirectory + 'pdbtemporary.pdb'
        
        cmd.load(proteinz, "proteinz")
        cmd.split_states("proteinz")
        cmd.delete("proteinz")
        cmd.set_name("proteinz_0001", "target")
        cmd.set_name("proteinz_0002", "homolog")
        cmd.color("red", "target")
        cmd.save(outputDirectory + pdb1 + "." + chain1 + "_" + pdb2 + "." + chain2 + ".cif")

        temppdb = outputDirectory + 'pdbtemporary.pdb'
        os.remove(temppdb)
        print("deleted temporary pdb")
        print("done")
        
    
    def ChainByChain(self, pdb1, pdb2):
        cmd.fetch(pdb1)
        for ch1 in cmd.get_chains(pdb1):
            cmd.delete(pdb1)
            cmd.fetch(pdb2)
            for ch2 in cmd.get_chains(pdb2):
                cmd.delete(pdb2)
                self.specificSingleChainAlignment(pdb1,pdb2,ch1, ch2)
                
    def fatcat_savexml(self, pdb1, pdb2):  
        
        os.system("bash runFATCAT.sh -pdb1 " + pdb1 + " -pdb2 " + pdb2 + " -pdbFilePath " + pdbDirectory + " -autoFetch -flexible -printFatCat -printXML " + outputDirectory + pdb1 + "_" + pdb2 + ".xml")
        print("done")
    ,
    

