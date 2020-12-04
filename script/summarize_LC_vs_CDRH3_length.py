#!/usr/bin/python
from collections import defaultdict

def LC_vs_CDRH3_count(in_filename, HC_germlines):
	infile = open(in_filename,'r')
	summary_dict = defaultdict(int)
	for line in infile.readlines():
		if "Antibody_ID" in line:
			continue
		else:
			entry = line.rstrip().rsplit("\t")
			HC_germline = entry[1]
			CDRH3 = entry[2]
			LC_germline = entry[3]
			Antigen_target = entry[5]
			if '-' in CDRH3 or 'length' in CDRH3:
				continue
			elif HC_germline in HC_germlines:
				CDRH3_len = len(CDRH3)-2
				ID = LC_germline+'__'+str(CDRH3_len)
				summary_dict[ID] += 1
	infile.close()
	return summary_dict

def write_output(summary_dict, out_filename):
	print ("writing file %s" % out_filename)
	outfile = open(out_filename,'w')
	outfile.write('LC_germline'+"\t"+"CDRH3_length"+"\t"+"count"+"\n")
	for ID in summary_dict.keys():
		LC_germline = ID.rsplit('__')[0]
		CDRH3_len = ID.rsplit('__')[1]
		count = str(summary_dict[ID])
		outfile.write(LC_germline+"\t"+CDRH3_len+"\t"+count+"\n")
	outfile.close()

def main():
	in_filename = "data/SARS-CoV-2_Abs_v12.tsv"
	out_filename = "result/LC_germline_vs_CDRH3_length.tsv"
	HC_germlines = ['IGHV3-53', 'IGHV3-66']
	summary_dict = LC_vs_CDRH3_count(in_filename, HC_germlines)
	write_output(summary_dict, out_filename)

if __name__ == "__main__":
	main()
