#!/usr/bin/python
import os
import sys
import glob

def ReadingData(filename,trial,datatype):
  DataDict = {}
  infile = open(filename,'r')
  printing = 0
  for line in infile.xreadlines():
    if printing == 1:
      line = line.rstrip().rsplit("\t")
      if len(line) < 3: continue
      if trial == 'Time1': time = line[0]; sign = line[1]; fit = line[2]
      if trial == 'Time2': time = line[3]; sign = line[4]; fit = line[5]
      if trial == 'Time3': time = line[6]; sign = line[7]; fit = line[8]
      if trial == 'Time4': time = line[9]; sign = line[10]; fit = line[11]
      if trial == 'Time5': time = line[12]; sign = line[13]; fit = line[14]
      DataDict[time] = fit if datatype == 'fit' else sign
    if 'Time1' == line[0:5]: printing = 1
  infile.close()
  return DataDict

def CompileData(AllData, outfile):
  outfile = open(outfile,'w')
  SampleIDs = AllData.keys()
  outfile.write("\t".join(['Time','Signal','SampleID'])+"\n")
  for ID in SampleIDs:
    for time in AllData[ID].keys():
      outfile.write("\t".join([time,AllData[ID][time],ID])+"\n")
  outfile.close()

def wrapper(In_Folder,Out_Folder,Exp,trial,datatype):
  filenames = glob.glob(In_Folder+'/'+Exp+'_*.txt')
  outfile   = Out_Folder+'/'+Exp+'_'+datatype+'_All.compile'
  AllData   = {}
  for filename in filenames:
    print "\tReading file: %s" % filename
    ID = os.path.basename(filename).replace('.txt','')
    AllData[ID] = ReadingData(filename,trial,datatype)
  CompileData(AllData, outfile)
  print "Compiled Data into %s" % outfile

def main():
  In_Folder    = 'data'
  Out_Folder   = 'result'
  #targets   = ['nCoVRBD','SARSRBD','nCoVSpike']
  #binders   = ['CR3022IgG','CR3022Fab']
  targets   = ['RBD']
  binders   = ['CC12.3-F58Y', 'CC12.3-WT', 'COV107-23-WT', 'COV107-23-swap', 'COVD21-C8-WT', 'COVD21-C8-swap',
               'COV107-23-F58Y','COVD21-C8-F58Y','COVA2-20-WT','COVA2-20-Y58F']
  Exps      = [target+'_'+binder for target in targets for binder in binders]
  datatypes  = ['fit','sign']
  Trials     = ['Time1']*len(Exps)
  for Trial, Exp in zip(Trials,Exps):
    for datatype in datatypes:
      wrapper(In_Folder,Out_Folder,Exp,Trial,datatype)

if __name__ == "__main__":
  main()
