import os
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

#extract reference sequence from the reference fasta file

ref_seq = SeqIO.read('hiv1cref.fasta', 'fasta').seq

#parse the CDS features from GFF file - defining and calling function
''' Adjusts start to 0-based (python slicing), retains end ( Python slices are end-exclusive). Extract all pol-encoding sequences (CDS) for the relevant reference genome from the GFF file '''

def parse_gff(gff_file):            #function definition
    cds = []                        #initialise an empty list       
    with open(gff_file) as f:       #opening the gff file  
        for line in f:              #iterating through file lines
            if line.startswith('#') or not line.strip(): continue       #skip comment and blank lines
            fields = line.strip().split('\t')                           #splitting each line into fields
            if len(fields) < 9 or fields[2].lower() not in ["cds", "gene"]: continue            #skip if line does not have all the necessary 9 fields and filter by CDS features only
            attrs = dict(pair.split('=',1) for pair in fields[8].split(';') if '=' in pair)     #split the attributes field into pair values
            gene = attrs.get('gene', attrs.get('Name', '')).lower()                             #extract gene name 
            if gene == 'protease' or gene == 'reversetranscriptase' or gene == 'integrase':     #condition for gene name - [here it is for pol region only (depends on usage intent)]
             cds.append({
                'gene': gene,
                'start': int(fields[3])-1,
                'end': int(fields[4]),
                'strand': fields[6]
                 })                                                                             #append the necessary features to the CDS list
    return cds                                                                                  #returns a list of CDS features as dictionaries

cds_features = parse_gff('hiv1cref.gff3')       #calling the function

print(cds_features)                             #printing the CDS list
print()

# annotation function
''' Process the rows from tsv files to generate the reference and alternate amino acids and the respective position. Returns the annotation column values to the calling statement to help generate variant position and list'''

def annotate_variant(row):                                           #receives individual rows from the tsv files
    pos = int(row['POS']) -1                                         #ensure base position as integer value and converts to 0-based for python indexing from 1-based gff/vcf format
    ref = row['REF'].upper()                                         #extract reference base from row and uppercase them
    alt = row['ALT'].upper()                                         #extracts alternate base from row and uppercase them
    for cds in cds_features:                                         #loop for iterating over CDS features
        if cds['start'] <= pos < cds['end']:                         #ensuring variant lies within the selected CDS regions
            strand = cds['strand']
            if strand == '+':                                        #ensuring the variant is on the positive strand 
                offset = pos-cds['start']                            #genomic offset - coordinate of the base downstream of the CDS start position
                codon_idx = offset // 3                              #codon index - calculates which codon the base lies in
                codon_start = cds['start'] + codon_idx * 3           #absolute starting genomic coordinate of the codon
                ref_codon = ref_seq[codon_start:codon_start+3]       #extract original codon from reference genome using codon_start position; +3 to ensure last base is added (python is end-exclusive)
                codon_list = list(str(ref_codon))                    #converts immutable string to mutable list format
                codon_baseidx = pos - codon_start                    #index of the codon's base which is mutated (0,1 or 2)
                codon_list[codon_baseidx] = alt                      #replace the mutated base in the reference codon list with the alternate base
                alt_codon = ''.join(codon_list)                      #rejoin the mutated codon 
            
            if len(ref_codon) != 3:                                                 #codon integrity check (incomplete codon, etc.)
                return pd.Series({'POS_AA_m': '-', 'MUT_AA_m': 'invalid-codon'})
            
            try:                    
                ref_aa = str(Seq(str(ref_codon)).translate())                       #translates reference codon
                alt_aa = str(Seq(alt_codon).translate())                            #translates alternate codon
                aa_pos = codon_idx + 1                                              #generate amino acid position and convert it to 1-based from python's 0-based
                mut_aa = f"{ref_aa}{aa_pos}{alt_aa}"                                #generate the mutation/variant column with reference aa, position of the aa, and alternate aa
                return pd.Series({'POS_AA_m': aa_pos, 'MUT_AA_m': mut_aa})          #return column values
            except Exception:                                                       #executed if translation fails due to some error
                return pd.Series({'POS_AA_m': '-', 'MUT_AA_m': 'translation-error'})
    return pd.Series({'POS_AA_m': '-', 'MUT_AA_m': 'non-encoding'})                 #marks non-encoding if variant is not within the selected CDS regions

# process tsv files and generate updated dataframe

for tsv_file in glob.glob("*ivar.tsv"):                         #loop for all input files (has to follow a specific naming pattern)
    print(f"Processing: {tsv_file}")                            #Progress feedback
    df = pd.read_csv(tsv_file, sep='\t')                        #read the input tsv file into a pandas dataframe
    aa_cols = df.apply(annotate_variant, axis=1)                #generate annotation columns by applying the called function by iterating over each row of the dataframe 
    outdf = pd.concat([df, aa_cols], axis=1)                    #concatenate the original dataframe and annotated columns as a new dataframe
    outfile = tsv_file.replace(".tsv", "_annotated.tsv")        #constructing the output filename for the annotated dataframe
    outdf.to_csv(outfile, sep='\t', index=False)                #writing annotated dataframe to new tsv file with specified name
    print(f"Processing complete: {tsv_file} -> {outfile}")      #Progress feedback
    print()
