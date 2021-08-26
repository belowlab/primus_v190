import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f1', '--fid1', help = 'fid1', dest = 'fid1')
parser.add_argument('-f2', '--fid2', help = 'fid2', dest = 'fid2')
args = parser.parse_args()

# ---------------------------------------------------------------

def convert(fid1, fid2):

    # 1. read in ersa file 
    # 2. get ersa estimate (integer) from the line that has the two individuals
    # 3. convert ersa integer into equivalent PI_HAT value
    # 4. PRINT out the PI_HAT value so that PRIMUS can use it in place of the old one for vector formation
    
    return 

def check_file(fid1, fid2):

    with open('/data100t1/home/grahame/projects/compadre/primus/primus_old/example_data/eur_small_vectorfile', 'r') as f:
        for line in f:
            if line.split('\t')[0] != "FID1":
                fid_1 = line.split('\t')[0]
                fid_2 = line.split('\t')[2]
                if fid_1 == fid1 and fid_2 == fid2:
                    print ('success')
                    return
        
    #print ('failure')
    return 



if __name__ == "__main__":

    if 'fid1' in args and "fid2" in args:

        #convert(args.fid1, args.fid2)
        check_file(args.fid1, args.fid2)

    else:

        print ('f up')



