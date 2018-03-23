#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 11:18:03 2018

@author: sven

Mainly copy & paste from the excellent documentation from pyteomics.
"""

from pyteomics import fasta, parser, mass, achrom, electrochem, auxiliary
import pandas as pd
import itertools
import numpy as np
import time

PROTON = 1.00727646658
MODS = {"cm":[57.021464, "H3C2NO"],
        "ox":[15.994915 , "O"]}

def digest(sequences, name_col, min_length=6, missed_cleavage=2,
           protease="trypsin"):
    """
    
    Digests a give FASTA sequence collection into peptides. The protease, 
    and minimum length can be configured via the options.
    
    Parameters:
    --------------------------
    sequences: pyteomics FileReader,
              Filereader object from an FASTA file.
              
    name_col: str,
                identifier column that is added to the final dataframe.
    min_length: int,
                minimal peptide length for sequences to be included in the 
                final dataframe.
    protease: str,
                Text identifier for a the desired protease
                
    Returns:
    ---------------
    df: df,
         A dataframe with the columns, sequence, proteine and DB
    """
    #store peptide (sequence) and proteins (by an integer id)
    peptides = []
    proteins = []
    
    for description, sequence in sequences:
        #get desc and cleavage products
        desc = fasta.parse(description)
        new_peptides = parser.cleave(sequence, parser.expasy_rules[protease],
                                     missed_cleavages=missed_cleavage,
                                     min_length=min_length)
               
        #store new data
        peptides.extend(new_peptides)
        proteins.extend([desc["id"]]*len(new_peptides))
    
    
    peptide_df = pd.DataFrame([peptides, proteins]).transpose()
    peptide_df.columns = ["Peptidesequence", "Protein"]
    peptide_df["DB"] = name_col
    return(peptide_df)


def get_mass(peptide, mass_dic={}, fixed={"C":57.021464}):
    """
    Compute mass of a peptide either from the sequence or from a dictionary
    look-ip
    """

    if peptide in mass_dic:
        pep_mass = mass_dic[peptide]

    else:
        #add modification masses
        add = 0
        for fixed_mod in fixed:
            add = peptide.count(fixed_mod) * fixed[fixed_mod]
        
        #compute pepmass
        pep_mass = mass.fast_mass(peptide)+ add
        mass_dic[peptide] = pep_mass
    return(pep_mass)
    

def compute_charged_mass(pepmass, z):
    """
    Computes the mz, given a peptide mass and a desired charge state-
    
    Parameters:
    ------------------
    pepmass: float,
             mass of the peptide
    z: int,
        charge state

    returns:
    ----------
    mz: float,
        mz
    """
    return((pepmass+(PROTON*z)) / z)
    
    
def generate_peptidepairs(df, cl_mass, idname):
    """
    Computes all pairwise combinations of peptides.
    """
    a = time.time()
    index = df.index.tolist()
    n = len(index)
    max_entries = int(n*(n+1.)/2)
    
    #init data
    ID_1 = np.zeros(max_entries)
    ID_2 = np.zeros(max_entries)
    masses = np.zeros(max_entries)
    mz1 =  np.zeros(max_entries)
    mz2 =  np.zeros(max_entries)
    mz3 =  np.zeros(max_entries)
    mz4 =  np.zeros(max_entries)
    
    #store redundant masses in dictionary
    mass_dic = {}

    #this will create all A-B combinations but not A-A
    idx = 0
    ten_percent = 0.1*max_entries
    fobj = open("{}_masses.csv".format(idname), "w")
    fobj.write("ID1,ID2,Mass,mz(z1),mz(z2),mz(z3),mz(z4)\r\n")
    print ("{}_masses.csv".format(idname))
    print ("{} total sequences to compute".format(max_entries))
    for i in itertools.combinations(index, 2):
        if idx % ten_percent == 0:
            print ("{}% done.".format(np.round(idx*100/max_entries), 2))
        ID_1[idx] = int(i[0])
        ID_2[idx] = int(i[1])
        pep1 = df.loc[i[0]]["Peptidesequence"]
        pep2 = df.loc[i[1]]["Peptidesequence"]
        
        #check if mass was previously computed
        mass1 = get_mass(pep1, mass_dic)
        mass2 = get_mass(pep2, mass_dic)
        cl_mass = mass1 + mass2 + cl_mass
        
        #store mz
        fobj.write("{},{},{},{},{},{},{}\r\n".format(
                int(i[0]), int(i[1]), cl_mass, 
                compute_charged_mass(cl_mass, 1),
                compute_charged_mass(cl_mass, 2),
                compute_charged_mass(cl_mass, 3),
                compute_charged_mass(cl_mass, 4)))
#        masses[idx] = cl_mass
#        mz1[idx] = compute_charged_mass(cl_mass, 1)
#        mz2[idx] = compute_charged_mass(cl_mass, 2)
#        mz3[idx] = compute_charged_mass(cl_mass, 3)
#        mz4[idx] = compute_charged_mass(cl_mass, 4)
        idx += 1

    res_df = pd.DataFrame([ID_1, ID_2, masses, mz1, mz2, mz3, mz4]).transpose()
    res_df.columns = ["ID1", "ID2", "mass", "mz(z1)", "mz(z2)", "mz(z3)", "mz(z4)"]
    res_df['ID1']=res_df['ID1'].astype(int)
    res_df['ID2']=res_df['ID2'].astype(int)
    b = time.time()
    print (np.round((b-a)/60.,2))
    return(res_df)
    
#set input params
idname = "test"
fasta_file = "test.fasta"
min_length = 6
protease = "trypsin"
crosslinker_mass = 138.06807961

#read database
target = fasta.read(fasta_file)
decoys = fasta.read(fasta_file)

#generate peptides
target_df = digest(target, "Target", min_length=min_length, protease=protease).head(100)
decoy_df = digest(decoys, "Decoy", min_length=min_length, protease=protease).head(100)
df = pd.concat([target_df, decoy_df])
df = df.reset_index()
df.to_csv("{}_PeptidesDigest.csv".format(idname))
pairs = generate_peptidepairs(df, crosslinker_mass, idname)

#df_masses = pd.read_csv("test_masses.csv")
