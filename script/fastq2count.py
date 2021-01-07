#!/usr/bin/python
import csv
from Bio import SeqIO
from collections import Counter

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def Ref2list(fasta):
    ref_records = SeqIO.parse(fasta,"fasta")
    Ref_ls =[]
    for record in ref_records:
        seq = str(record.seq)
        Ref_ls.append(seq)
    return Ref_ls


def Read2Mut2Count(fastq, Fprimer, Rprimer):
  print "reading %s" % fastq
  FPlth = len(Fprimer)
  RPlth = len(Rprimer)
  Rrecords = SeqIO.parse(fastq,"fastq")
  fastq_ls = []
  for record in Rrecords:
      ID = str(record.id)
      seq = str(record.seq)[FPlth:-RPlth]
      fastq_ls.append(seq)
  fastq_dict = Counter(fastq_ls)
  return fastq_dict

def main():
    F_primer = 'ACCTACAGATGAATTCTCTTAGGGCAGAAGATACCGCCGTCTACTACTGC'
    R_primer = 'GGGCCTTTTGTAGAAGCTGAACTCACAGTGACGGTAGTCCCTTGTCCCCA'
    ref = 'data/trimmed_CDRH3_pool.fa'
    ref_key = Ref2list(ref)

    for i in range(1,11):
        fastq = 'data/B38_CDRH3_lib/assembled_fastq/S'+ str(i) +'.assembled.fastq'
        dict = Read2Mut2Count(fastq,F_primer,R_primer)
        with open('result/S'+ str(i) + '_variant_count2.csv', 'w') as f:
            for key in ref_key:
                f.write("%s, %s, %s, %s\n" % (key, translation(key), len(key), dict[key]))
if __name__ == "__main__":
  main()

