import os, sys, json, csv
from logger import logger
log = logger().log


def isolate_subset_from_vectorfile(vectorFile, idList):

    # initialize smaller vector file
    outfile = csv.writer(open('../../../example_data/eur_small_vectorfile', 'w'), delimiter = '\t')
    outfile.writerow(['FID1'] + ['IID1'] + ['FID2'] + ['IID2'] + ['MOST_LIKELY_REL'] + ['LIKELIHOOD_VECTOR'] + ['IBD0'] + ['IBD1'] + ['IBD2'] + ['PI_HAT'] + ['-1'] + ['POSSIBLE_RELS'] + ['LIKELIHOOD_CUTOFF'])     # headers
    #outfile.close()

    # save list of IDS from famfile to memory
    idlist = []
    with open(idList, 'r') as famfile:
        for line in famfile:
            idlist.append(line.split(' ')[0])
    log(f'ID list populated, length {len(idlist)}', 'success')
    
    # iterate through BIG vector file and only add lines whose IDs correspond to ones that are in the idlist
    with open(vectorFile, 'r') as vectorfile, open('../../../example_data/eur_small_vectorfile', 'a') as outfile2:
        for line in vectorfile:
            if line.split('\t')[0] != "FID1":
                if line.split('\t')[0] in idlist or line.split('\t')[2] in idlist:
                    # print (line.split('\t'))
                    outfile2.write(line)
                    # sys.exit()
    log('Done!', 'success')            


def gen_new_vectorfile(smallVectorFile, ersafile, idlist):

    # Given a set of IIDs that were used in the PRIMUS run,
    # extract data on these and create a new matrix based on the likelihood vector one 
    
    # basically we want to keep the probability vector file as intact as possible 
    # and just replace the likelhood estimates with ERSA ones

    with open(smallVectorFile, 'r') as v:
        for line in v:
            if line.split('\t')[0] != "FID1":
                print (line.split('\t'))
                sys.exit()

    # at this point, i need to figure out how to replace the STRING likelihood values already in the vectorfile with INT values from ersa file 
    # once this is done, we can read this new file back into primus using the load_likelihood_vectors_from_file() method



if __name__ == "__main__":

    vectorFile = '../../../example_data/prob_vector_file_og'
    smallVectorFile = '../../../example_data/eur_small_vectorfile'

    eur_ersafile = '../../../example_data/ersa_eur.ersa'
    afr_erafile = '../../../example_data/ersa_afr.ersa'

    idList = '../../../example_data/eur_test.fam'

    # ------------------------------------------------------------
    
    #isolate_subset_from_vectorfile(vectorFile, idList)

    gen_new_vectorfile(smallVectorFile, eur_ersafile, idList) 