import pandas as pd
import numpy as np
import pickle
import pyranges as pr
import requests
import json

file_path ="..."

# load tcga clinical file 
tcga_sample_table = pd.read_excel(file_path + 'TCGA-CDR-SupplementalTableS1.xlsx').iloc[:, 1:]
tcga_sample_table['histological_type'].fillna('', inplace=True)

# load pathology annotations from asb
tumor_types = pd.read_csv(file_path + 'tumor_types_NCI-T.csv',sep='\t')

tumor_types.fillna('', inplace=True)
tcga_sample_table = pd.merge(tcga_sample_table, tumor_types[['type', 'histological_type', 'NCI-T Label', 'NCI-T Code']], how='left', on=['type', 'histological_type'])

# maf file with mutations
mc3_file_name = file_path+'/mc3.v0.2.8.CONTROLLED.maf'
usecols = ['Hugo_Symbol', 'Hugo_Symbol', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'STRAND', 'Variant_Classification', 'Variant_Type', 'Consequence', 'Reference_Allele', 'Tumor_Seq_Allele2', 't_ref_count', 't_alt_count', 'Tumor_Sample_Barcode', 'CONTEXT', 'FILTER', 'CDS_position']
tcga_maf = pd.read_csv(mc3_file_name, sep='\t', usecols=usecols, low_memory=False)
indels=["INS","DEL"] # or "SNP"
tcga_maf = tcga_maf[tcga_maf["Variant_Type"].isin(indels)] 
tcga_maf = tcga_maf.loc[tcga_maf['Chromosome'] != 'MT']

# df of counts via groupby
non_syn = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Nonstop_Mutation']
tcga_counts = tcga_maf[['Variant_Classification', 'Tumor_Sample_Barcode']].groupby('Tumor_Sample_Barcode').apply(lambda x: pd.Series([len(x), (x['Variant_Classification'].isin(non_syn)).sum()], index=['all_counts', 'non_syn_counts']))
tcga_counts['non_syn_tmb'] = tcga_counts['non_syn_counts'] / 31.85
tcga_counts.reset_index(inplace=True)
tcga_counts['bcr_patient_barcode'] = tcga_counts['Tumor_Sample_Barcode'].str.extract(r'^([^-]+-[^-]+-[^-]+)-')
tcga_counts['bcr_sample_barcode'] = tcga_counts['Tumor_Sample_Barcode'].str.extract(r'^([^-]+-[^-]+-[^-]+-[^-]+)-')

tcga_sample_table = pd.merge(tcga_sample_table, tcga_counts, how='right', on='bcr_patient_barcode')

# not every kit covered the entire exome
cases = list(tcga_sample_table['bcr_patient_barcode'])

cases_endpt = 'https://api.gdc.cancer.gov/cases'
responses = []
step = 100
count = 0
X = True
while X:
    print(count)
    if (count+1)*step >= len(cases):
        value = cases[count*step:]
        X = False
    value = cases[count*step: (count+1)*step]
    count += 1
    filt = {"op": "in",
            "content": {
                "field": "cases.submitter_id",
                "value": value
            }
            }
    params = {'filters': json.dumps(filt), "expand": "files.analysis.metadata.read_groups", 'size': '100'}
    response = requests.get(cases_endpt, params=params).json()
    responses.append(response)

flattened_responses = []
for response in responses:
    for i in response["data"]["hits"]:
        flattened_responses.append(i)

# map a sample to its barcode
case_to_barcode = {i: j for i, j in zip(tcga_sample_table['bcr_patient_barcode'], tcga_sample_table['Tumor_Sample_Barcode'])}
cancer_dict = {i: j for i, j in zip(tcga_sample_table['bcr_patient_barcode'], tcga_sample_table['type'])}

hits = {}
for case in case_to_barcode:
    for j in flattened_responses:
        for k in j['submitter_sample_ids']:
            if case in k:
                hits[case] = j
                break

data = {}
for case in case_to_barcode:
    for i in hits[case]['files']:
        for k in i['analysis']['metadata']['read_groups']:
            if k['library_strategy'] == 'WXS':
                Y = False
                sample = k['experiment_name']
                if sample[-2:] == '-1':
                    sample = sample[:-2]
                if sample == case_to_barcode[case]:
                    key = case_to_barcode[case]
                    Y = True
                else:
                    sample = sample[5:]
                    temp_sample = case_to_barcode[case][5:]
                    if sample == temp_sample:
                        key = case_to_barcode[case]
                        Y = True
                    if Y== False:
                        temp_sample = case_to_barcode[case][5:20]
                        if sample == temp_sample:
                            key = case_to_barcode[case]
                            Y = True
                if Y == True:
                    key = key[:15]
                    data[key] = data.get(key, {'centers': [], 'kits': [], 'beds': []})
                    data[key].update([('centers', data[key]['centers'] + [k['sequencing_center']]),\
                                         ('kits', data[key]['kits'] + [k['target_capture_kit_name']]),\
                                         ('beds', data[key]['beds'] + [k['target_capture_kit_target_region']])])


bad_kits=['Gapfiller_7m','NimbleGen Sequence Capture 2.1M Human Exome Array']
bad_beds=['https://bitbucket.org/cghub/cghub-capture-kit-info/raw/c38c4b9cb500b724de46546fd52f8d532fd9eba9/BI/vendor/Agilent/tcga_6k_genes.targetIntervals.bed',
'https://bitbucket.org/cghub/cghub-capture-kit-info/raw/c38c4b9cb500b724de46546fd52f8d532fd9eba9/BI/vendor/Agilent/cancer_2000gene_shift170.targetIntervals.bed']

bad_samples = []
null_samples = []

for sample in tcga_sample_table['Tumor_Sample_Barcode']:
    sample = sample[:15]
    if sample not in data:
        null_samples.append(sample)
    else:
        X = False
        for kit, bed in zip(data[sample]['kits'], data[sample]['beds']):
            if not kit:
                null_samples.append(sample)
            else:
                for sub_kit, sub_bed in zip(kit.split('|'), bed.split('|')):
                    if sub_kit not in bad_kits:
                        if sub_bed not in bad_beds:
                            X = True
                            break
        if X == False:
            bad_samples.append(sample)

# add columns to the sample table
tcga_sample_table['Exome_Covered'] = ~tcga_sample_table['Tumor_Sample_Barcode'].str[:15].isin(bad_samples + null_samples)
tcga_sample_table['Exome_Unknown'] = tcga_sample_table['Tumor_Sample_Barcode'].str[:15].isin(null_samples)

chromosomes = {}
for i in list(range(1, 23))+['X', 'Y']:
    with open(file_path + 'chromosomes' + '/' + ('chr' + str(i) + '.txt')) as f:
        chromosomes[str(i)] = f.read()

# reference genome
gff = pd.read_csv(file_path + 'Homo_sapiens.GRCh37.87.gff3',
                  sep='\t',
                  names=['chr', 'unknown', 'gene_part', 'start', 'end', 'unknown2', 'strand', 'unknown3', 'gene_info'],
                  usecols=['chr','gene_part', 'start', 'end', 'gene_info'],
                  low_memory=False)


gff_cds_pr = pr.PyRanges(gff.loc[(gff['gene_part'] == 'CDS') & gff['chr'].isin(chromosomes), ['chr', 'start', 'end', 'gene_info']].astype({'start': int, 'end': int}).rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'})).merge()
gff_exon_pr = pr.PyRanges(gff.loc[(gff['gene_part'] == 'exon') & gff['chr'].isin(chromosomes), ['chr', 'start', 'end', 'gene_info']].astype({'start': int, 'end': int}).rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'})).merge()
del gff

tcga_maf['index'] = tcga_maf.index.values

maf_pr = pr.PyRanges(tcga_maf.loc[:, ['Chromosome', 'Start_Position', 'End_Position', 'index']].rename(columns={'Start_Position': 'Start', 'End_Position': 'End'}))
genie = pd.read_csv(file_path + 'genomic_information.txt', sep='\t', low_memory=False)
panels = genie.SEQ_ASSAY_ID.unique()
panel_df = pd.DataFrame(data=panels, columns=['Panel'])
repeats = pd.read_csv(file_path + 'simpleRepeat.txt', sep='\t', low_memory=False, header=None, usecols=[1, 2, 3])
repeats[1] = repeats[1].str.replace('chr', '')
repeats.rename(columns={1: 'Chromosome', 2: 'Start', 3: 'End'}, inplace=True)
repeats_pr = pr.PyRanges(repeats.loc[repeats['Chromosome'].isin(chromosomes)]).merge()

total_sizes = []
cds_sizes = []
exon_sizes = []
panel_prs = []

for panel in panels:
    print(panel)
    panel_pr = pr.PyRanges(genie.loc[(genie['SEQ_ASSAY_ID'] == panel) & genie['Chromosome'].isin(chromosomes), 'Chromosome':'End_Position'].rename(columns={'Start_Position': 'Start', 'End_Position': 'End'})).merge()
    total_sizes.append(sum([i + 1 for i in panel_pr.lengths()]))
    cds_sizes.append(sum([i + 1 for i in panel_pr.intersect(gff_cds_pr).lengths()]))
    exon_sizes.append(sum([i + 1 for i in panel_pr.intersect(gff_exon_pr).lengths()]))
    panel_prs.append(panel_pr)

grs = {k: v for k, v in zip(['repeat', 'CDS', 'exon'] + list(panels), [repeats_pr, gff_cds_pr, gff_exon_pr] + panel_prs)}
result = pr.count_overlaps(grs, pr.concat({'maf': maf_pr}.values()))
result = result.df

tcga_maf = pd.merge(tcga_maf, result.iloc[:, 3:], how='left', on='index')

panel_df['total'] = total_sizes
panel_df['cds'] = cds_sizes
panel_df['exon'] = exon_sizes

agilent_df = pd.read_csv(file_path + 'whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed', sep='\t', low_memory=False, header=None)
kit_pr = pr.PyRanges(agilent_df.rename(columns={0: 'Chromosome', 1: 'Start', 2: 'End'})).merge()
kit_total = sum([i + 1 for i in kit_pr.lengths()])
kit_cds = sum([i + 1 for i in kit_pr.intersect(gff_cds_pr).merge().lengths()])
kit_exon = sum([i + 1 for i in kit_pr.intersect(gff_exon_pr).merge().lengths()])

panel_df = panel_df.append({'Panel': 'Agilent_kit', 'total': kit_total, 'cds': kit_cds, 'exon': kit_exon}, ignore_index=True)

chromosome_sizes = pd.read_csv(file_path + 'chromInfo.hg19.tsv', sep='\t', header=None)
chromosome_sizes.columns = chromosome_sizes.iloc[0]
chromosome_sizes = chromosome_sizes.drop(chromosome_sizes.index[0])
chromosome_sizes= chromosome_sizes[chromosome_sizes.index < 25] 
chrom_names = chromosome_sizes.hg19g0.tolist()
new_names = [x[3:] for x in chrom_names]
chromosome_sizes["hg19g0"] = new_names
chromosome_sizes['size'] = chromosome_sizes['size'].astype(int)
chromosome_sizes = dict(zip(chromosome_sizes['hg19g0'], chromosome_sizes['size']))
desired_order = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
sorted_chromosomes = sorted(chromosome_sizes.keys(), key=lambda x: desired_order.index(x))
cumulative_sizes = {}
total_genome_size = 0
for chrom in sorted_chromosomes:
    cumulative_sizes[chrom] = total_genome_size
    total_genome_size += chromosome_sizes[chrom]

def calculate_genome_wide_position(row):
    chrom = row['Chromosome']
    pos = row['Start_Position']
    return cumulative_sizes[chrom] + pos

tcga_maf['Genome_Wide_Position'] = tcga_maf.apply(calculate_genome_wide_position, axis=1)
tcga_maf['Normalized_Position'] = tcga_maf['Genome_Wide_Position'] / total_genome_size

# get mutation context, set to sequence length of 200 to get longer indels
def variant_features(maf, ref_length=200, alt_length=200, five_p_length=200, three_p_length=200):
    refs = []
    alts = []
    five_ps = []
    three_ps = []
    if ref_length % 2 != 0:
        ref_length += 1
    if alt_length % 2 != 0:
        alt_length += 1
    for index, row in enumerate(maf.itertuples()):
        Ref = row.Reference_Allele
        Alt = row.Tumor_Seq_Allele2
        Chr = str(row.Chromosome)
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

tcga_maf['Ref'], tcga_maf['Alt'], tcga_maf['five_p'], tcga_maf['three_p'] = variant_features(tcga_maf)
tcga_maf.drop(columns=['index'], inplace=True)
pickle.dump(tcga_maf, open(file_path + 'tcga_maf_table_id_mutation_catalogue.pkl', 'wb'),protocol=4)