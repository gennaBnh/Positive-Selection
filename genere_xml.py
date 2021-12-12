# Part of this  code is taken from https://groups.google.com/forum/#!topic/etetoolkit/cYTHHsL21KY
# written by Jaime Huerta-Cepas

"""
Ecrit un  arbre au format XML .
Usage:
  genere_xml.py <treeFile> <alignmentFile> <resultsFile>

Positional arguments:
  treeFile              family name / phylogenic trees in newick format
  alignmentFile         multiple alignment in fasta format
  resultsFile           statistics

"""

import argparse
import sys
import re
import random
import os
from Bio import Phylo
from io import StringIO
import time
from ete3 import Phyloxml, phyloxml

from lxml import etree
from xml.etree import ElementTree
from xml.dom import minidom

from lxml.etree import XMLParser, parse

# from docopt import docopt
import sqlite3
import zlib
import base64


# args = docopt(__doc__)
# tree = args["<treeFile>"]
# alignment = args["<alignmentFile>"]
# results = args["<resultsFile>"]
# arguments = docopt(__doc__, version='1.0.0rc2')

parser = argparse.ArgumentParser(description='Generate an XML file from \
    a tree, an alignment and statistics.')
parser.add_argument('-t', '--tree', dest='treeFile', action='store',\
    required=True,\
    help='file containing a tree in Newick format')
parser.add_argument('-a', '--align', dest='alignmentFile', action='store',\
    required=True,\
    help='file containing an alignment in FASTA format')
parser.add_argument('-s', '--stats', dest='resultsFile', action='store',\
    required=True,\
    help='file containing the results')
parser.add_argument('-c', '--col', dest='statcol', action='store', type=int,\
    default=1,\
    help='index of the results column to use')
parser.add_argument('-n', '--nostat', dest='nostat', action='store', type=float,\
    default=-1.0,\
    help='value to use in case there is no statistic associated\
    with a site in the sequence')
parser.add_argument('-o', '--output', dest='output', action='store',\
    help='name of the output XML file (if not specified, the XML will have \
    the same name as the tree file)')
args = parser.parse_args()


# def dna_to_prot(dna_seq:str):
def nuc_acid_to_prot(dna_seq:str):
    matches = {
        'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AGA':'R', 'AGC':'S', 'AGG':'R', 'AGT':'S',
        'ATA':'I', 'ATC':'I', 'ATG':'M', 'ATT':'I',

        'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAT':'H',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',

        'GAA':'E', 'GAC':'D', 'GAG':'E', 'GAT':'D',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',

        'TAA':'*', 'TAC':'Y', 'TAG':'*', 'TAT':'Y',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TGA':'*', 'TGC':'C', 'TGG':'W', 'TGT':'C',
        'TTA':'L', 'TTC':'F', 'TTG':'L', 'TTT':'F',

        'UAA':'*', 'UAC':'Y', 'UAG':'*', 'UAU':'Y',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UGA':'*', 'UGC':'C', 'UGG':'W', 'UGU':'C',
        'UUA':'L', 'UUC':'F', 'UUG':'L', 'UUU':'F',

        '---':'-'
    }
    aa = ''
    # if len(dna_seq)%3 == 0:
    codons = [dna_seq[i:i+3].upper() for i in range(0, len(dna_seq), 3)]
    for codon in codons:
        if len(codon) == 3:
            # codon = codon.replace('U', 'T')
            try:
                aa += matches[codon]
            except KeyError:
                print(f"Unknown codon: {codon}")
                return False
        else:
            print(f"Ignored: {codon} (not a codon)")
    return aa
    # else:
    #     print("Inorrect sequence length")
    #     return False


def loadAlignment(alignmentFile):
    """Loads a FASTA alignment."""
    alignmentDict = {}
    seq_idMax = 0
    with open(alignmentFile, 'r') as af:
        for line in af:
            if line[0] == '>': # FASTA header
                seq_id = line.strip('\n')[1:]
                if len(seq_id) > seq_idMax:
                    seq_idMax = len(seq_id)
                alignmentDict[seq_id] = ''
            else: # sequence
                aligned_seq = line.strip('\n')
                alignmentDict[seq_id] = aligned_seq
    return alignmentDict, seq_idMax


def loadResults(resultsFile, statcol=1, nostat_value=-1.0):
    """Loads statistics. <statcol> defines the 0-based index
    of the column to use (1 by default)."""
    resultsDict = {}
    with open(resultsFile, 'r') as f:
        i = 0
        for line in f:
            if i > 0: # avoid column headers
                line = line.strip('\n').split('\t')
                try:
                    site = int(line[0])
                    res = float(line[statcol])
                except ValueError:
                    print(f"Conversion failed in line\n\t{line}")
                else:
                    while not i == site: # in case there is no statistic for
                    # a site, define an aberrant value
                        print(f'Site {i} missing')
                        resultsDict[i] = nostat_value
                        i += 1
                    resultsDict[i] = res
            i += 1

    resultsText = ''
    for key, item in resultsDict.items():
        resultsText += f'{key}:{item}'
        if key != max(resultsDict.keys()):
            resultsText += ' '
    return resultsText


def loadDico(fileDico):
    """
    Fonction qui prend un fichier d'alignement en entree et retourne un dictionnaire des noms d'espece des sequences
    """

    file = open(fileDico,"r")
    speciesDico = {}
    for line in file:
        if line[0]==">":
            speciesDico[line.split('>')[1].split('\n')[0]] = line.split('>')[1].split('\n')[0]
            # tline = re.split('_',line)
            # #key = f"{tline[0]} {tline[1]}"
            # key = line.split('\n')[0]
            # fin = tline[4].split("\n")[0]
            # speciesDico[key.split(">")[1]] = f"{tline[2]} {tline[3]} {fin}"
    file.close()
    return speciesDico


def createPhyloXML(fam,newick):
    #Parse and return exactly one tree from the given file or handle
    if not ':' in newick:
        nv_arbre = ""
        for i in range(len(newick)):
            if (newick[i]==',' and newick[i-1]!=')') or (newick[i]==')' and newick[i-1]!=')'):
                nv_arbre+=":0.7"
                nv_arbre+=newick[i]
            elif (newick[i]==',' and newick[i-1]==')') or (newick[i]==')' and newick[i-1]==')'):
                nv_arbre+="30:0.4"
                nv_arbre+=newick[i]
            else:
                nv_arbre+=newick[i]
        newick = nv_arbre
    # print(newick)
    handle = StringIO(newick)
    trees = Phylo.read(handle, 'newick')
    #Write a sequence of Tree objects to the given file or handle
    rd = str(random.randint(0,1000))
    Phylo.write([trees], 'tmpfile-'+rd+'.xml', 'phyloxml')
    file = open('tmpfile-'+rd+'.xml', 'r')
    #Copie tous les objets dans une variable et supprime le fichier créé
    text = file.read()
    file.close()
    os.remove('tmpfile-'+rd+'.xml')
    #
    p = XMLParser(huge_tree=True)
    text = text.replace("phy:", "")

    text = re.sub("b'([^']*)'", "\\1", text)
    text = re.sub('branch_length_attr="[^"]+"', "", text)
    header = "<phyloxml>"

    text = re.sub('<phyloxml[^>]+>', header, text)
    text = text.replace('Phyloxml', 'phyloxml')
    tree = etree.fromstring(text,parser=p)
    # ajout du nom d'arbre
    treename = etree.Element("name")
    treename.text = fam
    ins = tree.find('phylogeny')
    ins.append(treename)

    clade = tree.xpath("/phyloxml/phylogeny/clade")
    subtree = tree.xpath("/phyloxml")
    nbfeuille = 0
    famspecies = {}

    for element in clade[0].iter('clade'): # pour chaque id de séquence
        enom=element.find('name')
        if (enom is not None) :
            nbfeuille = nbfeuille + 1
            cds = enom.text
            sp = dico.get(cds)
            if (not  sp):
                print ("undefined species for "+ cds)
                sp = "undefined"
            famspecies[sp] = 1

            ## Affectation de la séquence pour l'ID en cours
            seq_alg = alignmentDict.get(cds)
            if not seq_alg:
                print ("undefined alignment for "+ cds)
                seq_alg = ""

            evrec = etree.Element("eventsRec")
            leaf = etree.Element("leaf")
            if 'crossdico' in globals():
                crossref = etree.Element("crossref")
                if cds in crossdico[0]:
                    response = crossdico[0][cds]
                    tabbuf = response.split("|")
                    print (tabbuf)
                    for buf in tabbuf:
                        print ("process "+buf)
                        tabcross = buf.split(":")
                        if len(tabcross) > 1 :
                            crossref.set(tabcross[0], tabcross[1])
            leaf.set('speciesLocation', sp)
            if 'seqdefdico' in globals():
                if cds in seqdefdico:
                    leaf.set('defintiion', seqdefdico[cds])



            ## Ajout des séquences aux feuilles
            leaf.set('dnaAlign', seq_alg) # ajout de la séquence en nucléotides
            # leaf.set('aaAlign', translate.dna_to_prot(seq_alg)) # ajout de la séquence en acides aminés
            # leaf.set('aaAlign', dna_to_prot(seq_alg)) # ajout de la séquence en acides aminés
            leaf.set('aaAlign', nuc_acid_to_prot(seq_alg)) # ajout de la séquence en acides aminés

            if 'crossdico' in globals():
                leaf.append(crossref)
            evrec.append(leaf)
            element.append(evrec)
    print ("Number of leaves : ")
    print (nbfeuille)
    nbspecies = len(famspecies)
    print ("Number of species : ")
    print (nbspecies)

    ## Ajout des résultats
    resultsElement = etree.Element('statistics')
    resultsElement.set('results', resultsText)

    LengthMaxSeqID = etree.Element('maxSeqIdLength')
    LengthMaxSeqID.text = str(maxSeqIdLength)

    treesize =  etree.Element("size")
    treesize.set('leaves',str(nbfeuille))
    treesize.set('species',str(nbspecies))
    e=subtree[0].find('phylogeny')
    e.append(treesize)
    e.append(LengthMaxSeqID)
    e.append(resultsElement) # ajout de la balise contenant les résultats
    text =  minidom.parseString(ElementTree.tostring(subtree[0])).toprettyxml()
    # remove blank lines
    cleantext = "\n".join([ll.rstrip() for ll in text.splitlines() if ll.strip()])
    return cleantext


sys.setrecursionlimit(15000)
# print(sys.getrecursionlimit())

print ("Loading alignment... ")
alignmentDict =  loadAlignment(args.alignmentFile)[0]
maxSeqIdLength = loadAlignment(args.alignmentFile)[1]
print ("OK")

print ("Loading results... ")
resultsText =  loadResults(args.resultsFile, args.statcol, args.nostat)
print ("OK")

#Creates empty phyloxml document
#project = Phyloxml()   a decommenter si on veut un fichier xml unique

# Loads Species name dico
dico = loadDico(args.alignmentFile)

# Loads newick tree
treefile = open(args.treeFile,"r")
print(args.output)
if args.output:
    output_name = args.output
else:
    if '.' in args.treeFile:
        output_name = args.treeFile[::-1].split('.', 1)[1][::-1]
    else:
        output_name = args.treeFile
xmloutputfile = open(output_name+".xml","w")

for line in treefile:
    tline = re.split(' ',line)
    if len(tline) > 1:
        newick=tline[1]
        fam=tline[0]
    else:
        newick = tline[0]
        fam = ''
    phyloxmltree = createPhyloXML(fam,newick)
    xmloutputfile.write(phyloxmltree)
    print ("Famille "+fam+" OK")

treefile.close()
xmloutputfile.close()
