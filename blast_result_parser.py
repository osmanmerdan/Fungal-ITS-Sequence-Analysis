#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Author: Osman Merdan, MD
# Uludag University Medical School Department of Medical Microbiology 
# Last Updated : 6 July 2023
# Bioinformatic analysis pipeline for the project:
# Retrospective determination of species distribution and antifungal 
# susceptibility of Mucorales order fungi isolated from clinical specimens.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

import pandas as pd
import sys

# Read argument from command line
file = sys.argv[1]


# 1. Load blast output and name columns accordingly ------------------------------------------------------
# This step also add species name column to the output
dat = pd.read_table(file, delimiter='\t', header=None)
dat.columns = ['qacc', 'saccver', 'stitle', 'qlen', 'slen', 'length', 'qstart', 'qend',
                'sstart', 'send', 'qcovs', 'pident', 'nident', 'mismatch', 'gaps', 'bitscore', 'evalue']
dat.loc[:,'ssciname'] = dat.loc[:,'stitle'].str.split().str[:2].str.join(' ')
dat=dat.reindex(columns=['qacc', 'saccver','ssciname', 'stitle', 'qlen', 'slen', 'length', 'qstart', 'qend',
                'sstart', 'send', 'qcovs', 'pident', 'nident', 'mismatch', 'gaps', 'bitscore', 'evalue'])
# 2. Rule 1 ---------------------------------------------------------------------------------------------
# Any hit with pident>=96 and qcovs>=95 and evalue of 0.0 is identified as species level identification
# Normally we expect to see only one HSP passing these filters but some exceptions occur
# Take a copy of the original dataframe
R1 = dat.copy()

# Filter for Rule 1 criteria
R1 = R1.loc[(R1['pident']>=96) & (R1['qcovs']>=95) & (R1['evalue']==0), : ]
R1['idlevel']=''
R1['idgenus']=''
R1['idspecies']=''
print(f"******************************************************")
# If there is multiple HSP; qacc and species values must match in different rows
# Drop duplicate occurences (rows with same qacc and species names)
subset_cols1 = ['qacc', 'ssciname']
Multiple_HSP = R1.loc[R1.duplicated(subset=subset_cols1), ['qacc']]['qacc']

if len(Multiple_HSP)>0:
    print(f">> WARNING -- Samples with multiple HSP for same species: ")
    for i in Multiple_HSP:
        print(i)

R1.drop_duplicates(subset=subset_cols1, keep='first', inplace=True)

# Some mucorales species names in database have identical reference sequences
# Thats why hits have same 'length', 'slen', 'start', 'end', 'nident' values for those entries.
# Like Rhizopus azygosporus/ microsporus
# For those samples take the first hit drop duplicates
subset_cols2 = ['qacc', 'length', 'slen', 'qstart', 'qend', 'nident', 'qcovs']
Same_Alignment_Multiple_Species = R1.loc[R1.duplicated(subset=subset_cols2), ['qacc']]['qacc']

if len(Same_Alignment_Multiple_Species)>0:
    print(f">> WARNING -- Samples having same alignment with more than one species: ")
    print(f"In the parsed_blast_output file only one of those species names will be reported.")
    for i in Same_Alignment_Multiple_Species.drop_duplicates():
        print(i)

R1.drop_duplicates(subset=subset_cols2, keep='first', inplace=True)

# If there are multiple hits for different species which are passing filters
# Those results carry same qacc but other columns vary significantly
# Those samples must be tagged as unidentified, then drop duplicates which has same qacc values
subset_cols3 = ['qacc']
Multiple_Hits_Passing_Filters = R1.loc[R1.duplicated(subset=subset_cols3), ['qacc']]['qacc']

if len(Multiple_Hits_Passing_Filters)>0:
    print(f">> WARNING -- Samples having multiple species hits passing the filters: ")
    print(f"In the parsed_blast_output file only the best hit (according to pident and qcovs) will be reported")
    for i in Multiple_Hits_Passing_Filters.drop_duplicates():
        print(i)

#R1.loc[R1.duplicated(subset=subset_cols3), ['idspecies', 'idgenus', 'idlevel']] = 'unidentified'

# Name not-unidentified samples
R1.loc[R1['idgenus']!='unidentified','idgenus'] = R1.loc[R1['idgenus']!='unidentified','ssciname'].str.split().str[0]
R1.loc[R1['idgenus']!='unidentified','idspecies'] = R1.loc[R1['idgenus']!='unidentified','ssciname']
R1.loc[R1['idgenus']!='unidentified','idlevel'] ='species' 


R1=R1.drop_duplicates(subset= subset_cols3, keep='first')

# 3. Rule 2 ---------------------------------------------------------------------------------------------
# Remaining hits are considered unidentified.
R2 = dat.loc[
    (~dat['qacc'].isin(R1['qacc'])),:
    ].drop_duplicates(subset='qacc', keep='first')
if not R2.empty:
    R2.loc[:,'idlevel'] = 'unidentified'
    R2.loc[:,'idgenus'] = 'unidentified'
    R2.loc[:,'idspecies'] = 'unidentified'

# 4. Print information --------------------------------------------------------------------------------
print(f"******************************************************")
print(f"Unidentified number of hits: {len(R2)}")
print(f"Species level identification: {len(R1)}")
print(f"Total best hits evalueated: {len(R1)+len(R2)}")
print(f"******************************************************")
# 6. Concatenate and Write Results --------------------------------------------------------------------------------
results = pd.concat([R1,R2])
results.to_csv(path_or_buf='./blast-results/parsed_blast_result.tsv',
               sep='\t', header=True, index=False)