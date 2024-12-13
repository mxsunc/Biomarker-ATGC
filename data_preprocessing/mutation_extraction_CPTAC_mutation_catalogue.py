#%%
import pandas as pd
import numpy as np
import pickle
import pyranges as pr

file_path ="..."

usecols = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'case_id']

cptac_maf = pd.read_csv(file_path+'cohortMAF.ucec',usecols=usecols,sep='\t',
                        low_memory=False)
cptac_maf = cptac_maf.dropna(subset=['Start_Position', 'End_Position'])

indels=["INS","DEL"] # or "SNP"
cptac_maf = cptac_maf[cptac_maf["Variant_Type"].isin(indels)]

# df of counts via groupby
non_syn = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Nonstop_Mutation']
cptac_counts = cptac_maf[['Variant_Classification', 'Tumor_Sample_Barcode']].groupby('Tumor_Sample_Barcode').apply(lambda x: pd.Series([len(x), (x['Variant_Classification'].isin(non_syn)).sum()], index=['all_counts', 'non_syn_counts']))
cptac_counts['non_syn_tmb'] = cptac_counts['non_syn_counts'] / 31.85
cptac_counts.reset_index(inplace=True)

# chromosome information
chromosomes = {}
for i in list(range(1, 23))+['X', 'Y']:
    with open(file_path+'/' + 'chromosomes/' + ('chr' + str(i) + '.txt')) as f:
        chromosomes[str(i)] = f.read()

# reference genome
gff = pd.read_csv(file_path+'/Homo_sapiens.GRCh38.112.gff3',
                  sep='\t',
                  names=['chr', 'unknown', 'gene_part', 'start', 'end', 'unknown2', 'strand', 'unknown3', 'gene_info'],
                  usecols=['chr','gene_part', 'start', 'end', 'gene_info'],
                  low_memory=False)

gff_cds_pr = pr.PyRanges(gff.loc[(gff['gene_part'] == 'CDS') & gff['chr'].isin(chromosomes), ['chr', 'start', 'end', 'gene_info']].astype({'start': int, 'end': int}).rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'})).merge()
gff_exon_pr = pr.PyRanges(gff.loc[(gff['gene_part'] == 'exon') & gff['chr'].isin(chromosomes), ['chr', 'start', 'end', 'gene_info']].astype({'start': int, 'end': int}).rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'})).merge()
del gff

# make index column for merging
cptac_maf['index'] = cptac_maf.index.values

maf_pr = pr.PyRanges(cptac_maf.loc[:, ['Chromosome', 'Start_Position', 'End_Position', 'index']].rename(columns={'Start_Position': 'Start', 'End_Position': 'End'}))

chromosome_sizes = pd.read_csv(file_path+'/chromInfo.hg38.tsv', sep='\t', header=None)
chromosome_sizes.columns = chromosome_sizes.iloc[0]
chromosome_sizes = chromosome_sizes.drop(chromosome_sizes.index[0])
chromosome_sizes= chromosome_sizes[chromosome_sizes.index < 25] 
chrom_names = chromosome_sizes.hg38h0.tolist()
new_names = [x[3:] for x in chrom_names]
chromosome_sizes["hg38h0"] = new_names
chromosome_sizes['size'] = chromosome_sizes['size'].astype(int)
chromosome_sizes = dict(zip(chromosome_sizes['hg38h0'], chromosome_sizes['size']))
desired_order = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
sorted_chromosomes = sorted(chromosome_sizes.keys(), key=lambda x: desired_order.index(x))
cumulative_sizes = {}
total_genome_size = 0
for chrom in sorted_chromosomes:
    na = 'chr'+chrom
    cumulative_sizes[na] = total_genome_size
    total_genome_size += chromosome_sizes[chrom]

def calculate_genome_wide_position(row):
    chrom = row['Chromosome']
    pos = row['Start_Position']
    return cumulative_sizes[chrom] + pos

cptac_maf['Genome_Wide_Position'] = cptac_maf.apply(calculate_genome_wide_position, axis=1)
cptac_maf['Normalized_Position'] = cptac_maf['Genome_Wide_Position'] / total_genome_size

# get mutation context, set to sequence length of 200 to get longer indels
def variant_features(maf, ref_length=200, alt_length=200, five_p_length=200, three_p_length=200):
    refs = []
    alts = []
    five_ps = []
    three_ps = []
    if ref_length % 2 != 0:
        ref_length += 1
        print('Your ref length was not even, incrementing by 1.')
    if alt_length % 2 != 0:
        alt_length += 1
        print('Your alt length was not even, incrementing by 1.')

    for index, row in enumerate(maf.itertuples()):
        Ref = row.Reference_Allele
        Alt = row.Tumor_Seq_Allele2
        Chr = str(row.Chromosome)[3:]
        Start = row.Start_Position
        End = row.End_Position
        if pd.isna(Alt):
            print(str(index)+' Alt is nan')
            Ref = np.nan
            Alt = np.nan
            context_5p = np.nan
            context_3p = np.nan
        else:
            if len(Ref) > ref_length:
                Ref = Ref[:int(ref_length / 2)] + Ref[-int(ref_length / 2):]
            else:
                while len(Ref) < ref_length:
                    Ref += '-'
            if len(Alt) > alt_length:
                Alt = Alt[:int(alt_length / 2)] + Alt[-int(alt_length / 2):]
            else:
                while len(Alt) < alt_length:
                    Alt += '-'
            if row.Reference_Allele == '-':
                ##the TCGA coordinates for a null ref are a little weird
                assert Start-five_p_length >= 0
                context_5p = chromosomes[Chr][Start-five_p_length:Start]
                context_3p = chromosomes[Chr][Start:Start+three_p_length]
            else:
                assert Start-(five_p_length+1) >= 0
                context_5p = chromosomes[Chr][Start-(five_p_length+1):Start-1]
                context_3p = chromosomes[Chr][End:End+three_p_length]
        refs.append(Ref)
        alts.append(Alt)
        five_ps.append(context_5p)
        three_ps.append(context_3p)
    return refs, alts, five_ps, three_ps

cptac_maf['Ref'], cptac_maf['Alt'], cptac_maf['five_p'], cptac_maf['three_p'] = variant_features(cptac_maf)
cptac_maf.drop(columns=['index'], inplace=True)
pickle.dump(cptac_maf, open(file_path + 'cptac_ucec_maf_table_id_mutation_catalogue.pkl', 'wb'),protocol=4)
