# Installing & Running FATCAT 1.0 (required to run ImportGlob method)
Running FATCAT 1.0 (include citation) will allow users to run the flexible protein alignment tool needed to run some of the methods 
Recommended to verify that FATCAT is running properly prior to running methods such as ‘ImportGlob’

Materials
Download the folder ‘ProteinComparisonTool’ from Github 

Before Beginning 
Have Anaconda & JAVA installed
Anaconda: https://docs.anaconda.com/anaconda/install/windows/ 
JAVA: https://www.java.com/en/download/
(For Windows OS) have WSL installed
(WSL installation manual) https://docs.microsoft.com/en-us/windows/wsl/install-manual#step-3---enable-virtual-machine-feature 

Running FATCAT
in the command prompt enter:
bash runFATCAT.sh -pdb1 pdbID(enter pdb id here) -pdb2 pdbID(enter 2nd pdbID here) -pdbFilePath /(enter path to pdb files here)/ -autoFetch -printFatCat -flexible

Ex:
bash runFATCAT.sh -pdb1 1ake -pdb2 2ake -pdbFilePath /Users/name/Documents/pdb_file_folder/ -autoFetch -printFatCat -flexible


### protein3dcompare class

# Installing
Clone the Protein Comparison Tool repository

# Dependencies
Install the latest version of Java: https://www.java.com/en/download/
Install Pymol: https://pymol.org/2/

Make sure the following libraries are installed: 
requests 
BeautifulSoup
pandas

# Main Methods
*Enter pdbDirectory and outputDirectory as prompted

getHomologData(genename, targetprotein, pdbidlist)

Returns a dataframe and excel spreadsheet ("genenamedata.xlsx") including pdb id, uniprot id, organism, pdb entry name, global stoichiometry, length, associated pubmed ID(s), and Taxonomy for each pdb in pdbidlist

# FATCAT Methods
importGlob(pdb1, pdb2)
Runs FATCAT and saves alignment as a cif file

ChainByChain(pdb1,pdb2)
Runs FATCAT between every chain in both PDB inputs and returns cif file for every alignment

fatcat_savexml(pdb1,pdb2)
Runs FATCAT and saves alignment as an xml file

