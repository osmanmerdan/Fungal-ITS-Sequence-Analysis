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
# 2. Rule ---------------------------------------------------------------------------------------------
# Any hit with pident>=96 and qcovs>=95 and evalue of 0.00 is identified as species level identification
# There may be multiple hits because of multiple hsp's.
# Thats why duplicate occurences of qacc is discarded and best hits are saved as R1.
R1 = dat.loc[
    (dat['pident']>=96) & (dat['qcovs']>=95) & (dat['evalue']==0),:
    ].drop_duplicates(subset='qacc', keep='first')
R1.loc[:,'idlevel'] = 'species'
R1.loc[:,'idgenus'] = R1.loc[:,'ssciname'].str.split().str[0]
R1.loc[:,'idspecies'] = R1.loc[:,'ssciname']


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