#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 16:20:18 2015
@author: luca zoccarato

Useful to implement fasta36 tools in your qiime-like pipeline
- grep the best hit for each querried seq using evalue/identity/score criteria
- in case of equality, it parses the taxonomy of those hits and retains the last common ancestor

USAGE: parse_fasta36_output.py X Y W Z
X -> GGSEARCH36 text output file (with -m BB option)
Y -> filtering criteria evalue/identity/score
W -> e-value threshold
Z -> Y/N whether the taxonomy is present in the header of the reference database or not
e.g. "parse_fasta36_output.py rep_seq_gls36_tax.txt identity 1 N"
"""

import re
from sys import argv

fasta36 = open(argv[1], "U").readlines()
filtering_criteria = argv[2]
e_threshold = float(argv[3])
tax_yn = argv[4]

output_filename = '.'.join(argv[1].split('.')[:-1])+'_parsed.'+argv[1].split('.')[-1]
out = open(output_filename, 'w')

# get list of assigned otus
seq_id = []
for line in fasta36:
    if re.match('Query=', line): # identify row with queried sequence
        seq_id.extend([re.sub('Query= ','',line.strip())])

# retrive taxonomies, evalues/scores and taxonomy labels
tax_id = [[] for x in range(len(seq_id))]
tax_hits = [[] for x in range(len(seq_id))]
evalue = [[] for x in range(len(seq_id))]  
score = [[] for x in range(len(seq_id))]
identity = [[] for x in range(len(seq_id))] 
count = -1
for line in range(len(fasta36)):
    if re.match('Query=', fasta36[line]): # query counter
        count += 1
    if re.match('>', fasta36[line]): # best taxonomic hits for query
        tax = fasta36[line].split()
        if len(tax) == 1:  # solve issue with long tax_id when the tax_hit is in the next row
            tax.extend([fasta36[line+1].strip()]) 
        if tax_yn == "Y":
            tax_hits[count].extend([tax[1]])
        elif tax_yn == "N":
            tax_hits[count].extend(["no taxonomy in ref DB"])
        tax_id[count].extend([re.sub('>','',tax[0])])
    elif re.search('Expect = ', fasta36[line]):
        evalue[count].extend([float(fasta36[line].split('Expect = ')[-1].strip())]) # for e-value
        score[count].extend([float(fasta36[line].split('Score = ')[-1].split(' ')[0].strip())]) # for score
    elif re.search('Identities', fasta36[line]): # identities for taxonomic hits
        identity[count].extend([float(fasta36[line].split('(')[1].split('%')[0])])

# select filtering criteria
if filtering_criteria == "evalue":
    pick_criteria = evalue
elif filtering_criteria == "identity":
    pick_criteria = identity
elif filtering_criteria == "score":   
    pick_criteria = score
else:
    print "Please enter as filtering criteria: evalue/identity/score"

# parse output
for otus in range(len(seq_id)):
    if filtering_criteria == "evalue" and pick_criteria[otus][0] < e_threshold or filtering_criteria in ("identity","score") and pick_criteria[otus][0] > e_threshold:
        # when there is only ONE hit 
        if len(pick_criteria[otus]) == 1:
            out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (seq_id[otus], tax_hits[otus][0], evalue[otus][0], identity[otus][0], score[otus][0], tax_id[otus][0]))
    
        elif len(pick_criteria[otus]) > 1:
            # when there is a best hit 
            if pick_criteria[otus][1] > min(pick_criteria[otus]): # there is a best hit
                out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (seq_id[otus], tax_hits[otus][0], evalue[otus][0], identity[otus][0], score[otus][0], tax_id[otus][0]))
          
            # when there are more best hits - equality
            else:
                equal_ids = [i for i, x in enumerate(pick_criteria[otus]) if x == min(pick_criteria[otus])] #ids of equivalent hits
                # fix unkown issue that gives +1 number of evalues (e.g. for 10 hits give 11 evalues for that otu)
                equal_ids = [z for z in equal_ids if z <= 10]
                                
                # equivalent taxonomies
                tax_hits_new = [] 
                for t in equal_ids:
                    tax_hits_new.extend([tax_hits[otus][t]])
                # parse equival tax to keep only minimum common ancestor
                parsed_tax = [] 
                for letter in range(len(min(tax_hits_new, key=len))): 
                    n_letter = []                    
                    for x in tax_hits_new:
                        n_letter.extend([x[letter]])
                    if all(x==n_letter[0] for x in n_letter):
                        parsed_tax += n_letter[0]
                    else:
                        break
                parsed_tax = ''.join(parsed_tax)
                if len(parsed_tax) == 0:
                    parsed_tax = 'totally_different_ancestor'
                
                # average evalue
                av_evalue = []
                for t in equal_ids:
                    av_evalue.extend([evalue[otus][t]])
                av_evalue = sum(av_evalue)/len(av_evalue)
                
                # average identity
                av_identity = []
                for t in equal_ids:
                    av_identity.extend([identity[otus][t]])
                av_identity = sum(av_identity)/len(av_identity)
                
                # average score
                av_score = []
                for t in equal_ids:
                    av_score.extend([score[otus][t]])
                av_score = sum(av_score)/len(av_score)
                
                # equivalent tax identifiers
                eq_tax_id = [] 
                for t in equal_ids:
                    eq_tax_id.extend([tax_id[otus][t]])
                eq_tax_id = '\t'.join(eq_tax_id)
                
                out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (seq_id[otus], parsed_tax, av_evalue, av_identity, av_score, eq_tax_id))
