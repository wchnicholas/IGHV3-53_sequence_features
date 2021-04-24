#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from scipy import stats
from collections import Counter

def CsvWithHeader2Hash(fitfile):
  H = {}
  infile = open(fitfile,'r')
  countline = 0
  header = []
  for line in infile.xreadlines():
    countline += 1
    line = line.rstrip().rsplit(",")
    if countline == 1: header = line; continue
    mut = line[0].replace(' ','')
    H[mut] = {}
    for i in range(1,len(line)): H[mut][header[i]] = line[i]
  infile.close()
  return H

def extract_variant(Enrich_dict):
  CDRH3_list = {'non-enriched':[],'enriched':[]}
  for CDRH3 in Enrich_dict.keys():
    bind1 = float(Enrich_dict[CDRH3]['IGGK1-9_bind_rep1_RF'])
    bind2 = float(Enrich_dict[CDRH3]['IGGK1-9_bind_rep2_RF'])
    bind  = np.mean([bind1, bind2])
    CDRH3 = CDRH3[2::]
    if bind > 0:
      CDRH3_list['enriched'].append(CDRH3)
    else:
      CDRH3_list['non-enriched'].append(CDRH3)
      if len(CDRH3) == 9:
        print CDRH3
  return (CDRH3_list)

def out_CDRH3_length(CDRH3_list, length_file):
  outfile = open(length_file, 'w')
  outfile.write("\t".join(['class','CDRH3_length','count'])+"\n")
  CDRH3_length_enriched_dict    = Counter([len(CDRH3) for CDRH3 in CDRH3_list['enriched']])
  CDRH3_length_nonenriched_dict = Counter([len(CDRH3) for CDRH3 in CDRH3_list['non-enriched']])
  for CDRH3_length in CDRH3_length_enriched_dict.keys():
    count = CDRH3_length_enriched_dict[CDRH3_length]
    outfile.write("\t".join(map(str,['enriched',CDRH3_length, count]))+"\n")
  for CDRH3_length in CDRH3_length_nonenriched_dict.keys():
    count = CDRH3_length_nonenriched_dict[CDRH3_length]
    outfile.write("\t".join(map(str,['nonenriched',CDRH3_length, count]))+"\n")
  outfile.close()

def main():
  Enrich_dict = CsvWithHeader2Hash('result/CDRH3_Lib_RFindex.csv')
  length_file = 'result/CDRH3_length_dist.tsv'
  CDRH3_list  = extract_variant(Enrich_dict)
  out_CDRH3_length(CDRH3_list, length_file)
  CDRH3_len_enriched = [len(CDRH3) for CDRH3 in CDRH3_list['enriched']]
  CDRH3_len_nonenriched = [len(CDRH3) for CDRH3 in CDRH3_list['non-enriched']]
  print stats.ttest_ind(CDRH3_len_enriched, CDRH3_len_nonenriched)

if __name__ == "__main__":
  main()
