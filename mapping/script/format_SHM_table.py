#!/usr/bin/python
import os
import sys
from collections import Counter

def format_SHM(SHM_raw):
  if 'none' in SHM_raw:
    return []
  else:
    SHM_list = SHM_raw.rsplit(',')
    SHM_new_list = []
    for SHM in SHM_list:
      SHM = SHM.replace(' ','')
      if SHM[-2].isalpha():
        SHM = SHM[0]+str(int(SHM[1:-2]))+SHM[-2].lower()+SHM[-1]
      else:
        SHM = SHM[0]+str(int(SHM[1:-1]))+SHM[-1]
      SHM_new_list.append(SHM) 
    return SHM_new_list

def parse_SHM_file(filename):
  infile = open(filename, 'r')
  ab_info_dict = {}
  for line in infile.readlines():
    if 'Antibody' in line: continue
    entry = line.rstrip().rsplit("\t")
    Ab_name = entry[0]
    HC_germline = entry[1]
    SHM         = entry[2]
    CDRH3 = entry[3]
    LC_germline = entry[4].rsplit('*')[0]
    if '-' in CDRH3:
      continue
    if 'sequence' in SHM: 
      continue
    if 'length' in CDRH3:
      CDRH3_len = int(CDRH3.rsplit('=')[1])
    else:
      CDRH3_len = len(CDRH3)
    SHM = format_SHM(SHM)
    ab_info_dict[Ab_name] = {'LC':LC_germline, 
                             'CDRH3_length':CDRH3_len,
                             'SHM':SHM}
  return ab_info_dict

def SHM_positioning(SHM):
  if SHM[-2].isalpha():
    return int(SHM[1:-2])
  else:
    return int(SHM[1:-1])

def output_heatmap_format(ab_info_dict, SHM_list, outfile):
  print ("writing: %s" % outfile)
  outfile = open(outfile,'w')
  header  = ['Ab', 'LC',  'CDRH3_length']+SHM_list
  outfile.write("\t".join(header)+"\n")
  for Ab in ab_info_dict.keys():
    LC           = ab_info_dict[Ab]['LC']
    CDRH3_length = ab_info_dict[Ab]['CDRH3_length']
    ab_SHMs      = ab_info_dict[Ab]['SHM']
    out = [Ab, LC, CDRH3_length]
    for SHM in SHM_list:
      if SHM in ab_SHMs:
        out.append(1) 
      else:
        out.append(0)
    outfile.write("\t".join(map(str,out))+"\n")
  outfile.close()

def main():
  filename = "data/VH3-53_RBD_antibodies.tsv"
  outfile  = "result/AntibodyDataSM_heatmap.tsv"
  ab_info_dict = parse_SHM_file(filename)  
  SHM_list = []
  [SHM_list.extend(ab_info_dict[ab]['SHM']) for ab in ab_info_dict.keys()]
  SHM_list = sorted(list(set(SHM_list)), key=lambda x:SHM_positioning(x))
  output_heatmap_format(ab_info_dict, SHM_list, outfile)
  length_dict = Counter([ab_info_dict[ab]['CDRH3_length'] for ab in ab_info_dict.keys()])
  print (sum(length_dict.values()))
  for length in sorted(length_dict.keys(), reverse=True):
    print (length, length_dict[length])

if __name__ == "__main__":
  main()
