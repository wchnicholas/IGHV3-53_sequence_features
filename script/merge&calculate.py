import numpy as np
import pandas as pd

def cal_RFindex(count_DF,column,input):
    count_DF[column+'_RFindex'] = np.log10(((count_DF[column]+1)/(count_DF[column].sum(0)+143))/((count_DF[input]+1)/(count_DF[input].sum(0)+143)))
    return count_DF

def cal_mean(count_DF,mean,column1,column2):
    count_DF[mean] = (count_DF[column1]+count_DF[column2])/2
    return count_DF
def main():
    S1_DF = pd.read_table('result/S1_variant_count2.csv', header=None, sep= ',')
    S1_DF.columns = ['DNA','AA','DNA_LENGTH','S1count']
    count_DF = pd.DataFrame(S1_DF.loc[:,'AA'])
    for i in range(1,6):
        filename = 'result/S'+ str(i) + '_variant_count2.csv'
        readfile = pd.read_table(filename, header=None, sep= ',')
        readfile.columns = ['DNA','AA','DNA_LENGTH','S'+str(i)+'count']
        sample_DF = readfile.loc[:,'S'+str(i)+'count']
        count_DF = count_DF.join(sample_DF)
        print "processing file: %s" % (filename,)
    count_DF.columns = ['AA','input', 'exp_rep1','exp_rep2','bind_rep1','bind_rep2']#make CDRH3 variant count table

    cal_RFindex(count_DF,'exp_rep1','input') #calculating and adding a new column of the exp_rep1 RFindex
    cal_RFindex(count_DF,'exp_rep2','input')
    cal_RFindex(count_DF,'bind_rep1','input')
    cal_RFindex(count_DF,'bind_rep2','input')
    cal_mean(count_DF,'exp_average_RFindex','exp_rep1_RFindex','exp_rep2_RFindex')
    cal_mean(count_DF,'bind_average_RFindex','bind_rep1_RFindex','bind_rep2_RFindex')
    count_DF.to_csv('result/CDR3_LIB_index.csv') #save dataframe to csv

if __name__ == "__main__":
  main()
