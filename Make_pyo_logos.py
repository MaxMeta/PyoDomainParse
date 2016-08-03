# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 17:10:35 2016

@author: owenje
"""
import xml.etree.ElementTree as ET
import os
from Bio import SeqIO
import re

def PullOutProts(directory,OutName, TermsToMatch = ["pyoverdine", "pvd","Pvd","PVD","Pyoverdine", "PYOVERDINE"],extensions = {".faa"}):
    """
    searches a folder full of multi fasta files, finds all sequences matching
    one or more of the terms in TermsToMatch and writes these to a multifasta
    output file.
    """
    names = os.listdir(directory)           
    in_list = [directory+name for name in names if name[-4:] in extensions]
    Terms_RE = re.compile("|".join(TermsToMatch))
    matched = []
    for genome in in_list:
        Parsed = SeqIO.parse(open(genome,"rU"),"fasta")
        for record in Parsed:
            if Terms_RE.search(record.description):
                matched.append(record)     
    with open(OutName,"w") as F:
        SeqIO.write(matched,F,"fasta")

#PullOutProts("directory/containing/protein/fastas","output_name.faa")
#ran above on proteins from all pseudomonas genomes then passed these to antismash

def ParseXML(in_list):
    """
    Takes a list of filenames for antismash XML outputs. Parses these and returns
    a list of tuples containing two items: 
           1) PCP sequences from that XML, and 
           2) a dictionary of infor for all domains in which the key is the domain number
    The output list can be passed to BinByDomain to bin PCP sequenes by their
    downstream interaction domain
    """
    OutList = []
    for XML in in_list:
        print XML
        PCP_Seqs = []
        All_domains = {}
        try:
            parsed = ET.parse(XML)
        except Exception:
            print XML +" is broken"
            continue
        domains = parsed.findall('.//domain')
        for domain in domains:
            number = int(domain.attrib['id'])
            DomainType = domain.find('label').text
            location = domain.find('location')
            geneID = int(location.find('gene').find('geneid').text)
            All_domains[number] = {'source gene':geneID, 'domain type':DomainType}
            if DomainType == 'PCP':
                seq = location.find('protein').find('sequence').text
                PCP_Seqs.append({'sequence':seq,'source gene':geneID, 'domain number':number})
        OutList.append((PCP_Seqs,All_domains))
    return OutList


def BinByDomain(OutList):
    """
    Takes the output from Parse XML and outputs a ditionary where each PCP sequence
    is binned by it's downstream interaction partner (if Cis) of binned as Trans
    if no interaction domain is downstream on the same protein. When the down-
    stream domain has no label in the XML file, the PCP is binned as other
    Currrently set up for C, A, E, TE domain in Cis
    """
    C_down = []
    A_down = []
    E_down = []
    TE_down = []
    Trans = []
    Other = []
    for tup in OutList:
        PCP_Seqs = tup[0]
        All_domains = tup[1]
        for PCP in PCP_Seqs:
            try:#necessary to avoid KeyError if last seq is a PCP
                if All_domains[PCP['domain number'] +1]['source gene'] == PCP['source gene']:
                    if All_domains[PCP['domain number'] +1]['domain type'] == 'C':
                        C_down.append(PCP['sequence'])
                    elif All_domains[PCP['domain number'] +1]['domain type'] == 'A':
                        A_down.append(PCP['sequence'])            
                    elif All_domains[PCP['domain number'] +1]['domain type'] == 'E':
                        E_down.append(PCP['sequence'])
                    elif All_domains[PCP['domain number'] +1]['domain type'] == 'TE':
                        TE_down.append(PCP['sequence'])
                    else: Other.append((PCP['sequence']))
                else: Trans.append(PCP['sequence'])
            except KeyError: 
                Trans.append(PCP['sequence'])
    return {'CD':C_down,'AD':A_down,'ED':E_down,'TE':TE_down,'TR':Trans,'OT':Other}



def SeqDict2Fasta(SeqDict,Outname,include = {'CD','AD','TE','ED'}):
    """
    Takes the dictionary ouptut from BinByDomain and creates a multi-Fasta file
    in which each sequence has a header of the format: >DomainType_sequence_number
    For example: >TE_Sequence_1
    """   
    F = open(Outname,'w')
    i = 1
    for Domain in SeqDict:
        if Domain in include:
            prefix = Domain
            for Seq in SeqDict[Domain]:
                F.write( ">%s_Sequence_%s\n"%(prefix, i))
                F.write(Seq+'\n')
                i+=1
    F.close()
        

def WriteSubAlignments(AlignFile,TargetSet = {'CD','TE','ED'}):
    """
    Takes an alignment file and writes seperate alignments for each target domain
    in the set TargetSet. These can then be passed to weblogo to generate sequence
    logo images
    """
    AF = open(AlignFile,'r')
    tmp_dict = {target:[]for target in TargetSet}
    data = 'place holder'
    to_write = []
    target = 'CD'
    while data:
        data = AF.readline()
        #print data
        if len(data) >1 and data[0] == '>':
            if data != 'place holder' and target in TargetSet:
                tmp_dict[target].append(''.join(to_write))
            to_write = []
            target = data[1:3]
        to_write.append(data)
    #return tmp_dict
    for target in tmp_dict:
        fname = "%s_sub_align.txt"%target
        with open(fname,'w') as F:
            for line in tmp_dict[target]:
                F.write(line)
    AF.close()
        

def RunAll(directory,Fasta_Name,Cluster_Name,Align_Name):
    """
    Wrapper function that:
    -Runs ParseXML() to pull out PCPs, 
    -Runs BinByDomain to categorize these by theirdownstream interaction partner 
    -runs SeqDict2Fasta() to write the fasta outputs
    -Clusters the fasta outputs using usearch to dereplicate, default is  id = 0.95
    -Aligns the dereplicated files using muscle. 
    -Writes sub alignments for each target domain using WrtieSubAlignments(). 
    """
    names = os.listdir(directory)
    in_list = [directory+name for name in names if name[-4:] ==".xml"]    
    OutList = ParseXML(in_list)
    SeqDict = BinByDomain(OutList)
    SeqDict2Fasta(SeqDict,Fasta_Name)
    os.system("usearch -cluster_fast %s -id 0.95 -centroids %s"%(Fasta_Name,Cluster_Name))
    os.system("muscle -in %s -out %s"%(Cluster_Name,Align_Name))
    WriteSubAlignments(Align_Name)

directory = "/path/tp/directory/containing/xml_files"
Fasta_Name = "PCP_fasta_file.fasta"
Cluster_Name = "Dereplicated_PCPs.fasta"
Align_Name = "Aligned_PCPs.txt"
#os.chdir(directory)
#RunAll(directory,Fasta_Name ,Cluster_Name,Align_Name)
#Resutling subalignments are then used to generate sequence logos.
