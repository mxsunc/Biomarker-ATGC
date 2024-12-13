#%%
import logomaker
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_sequence_logos_for_cluster(sequences_5p, sequences_ref, sequences_alt, sequences_3p):
    """Get the sequence logos of all four sequence blocks."""
    fig, axes = plt.subplots(1, 4, figsize=(20, 5), tight_layout=True)
    sequence_sets = [sequences_5p, sequences_ref, sequences_alt, sequences_3p]
    titles = ["5' Sequence", "Reference Sequence", "Alteration Sequence", "3' Sequence"]
    for ax, sequences, title in zip(axes, sequence_sets, titles):
        frequencies = sequences_to_df(sequences)
        logomaker.Logo(frequencies, ax=ax)
        ax.set_title(title)
        ax.set_ylabel('Frequency')
        ax.set_xlabel('Position')
        ax.set_ylim(0,2)
    plt.show()

def sequences_to_df(sequences):
    """Convert a list of sequences into a DataFrame suitable for LogoMaker."""
    # Initialize a DataFrame to store nucleotide counts
    counts = pd.DataFrame(0, index=np.arange(len(sequences[0])), columns=['A', 'C', 'G', 'T'])
    
    # Count nucleotides at each position
    for seq in sequences:
        for i, nucleotide in enumerate(seq):
            if nucleotide in counts.columns:
                counts.loc[i, nucleotide] += 1

    counts += 1

    frequencies = counts.divide(counts.sum(axis=1), axis=0)
    ic = np.log2(4) + np.sum(frequencies * np.log2(frequencies), axis=1)
    ic_df = frequencies.multiply(ic, axis='index')
    
    return ic_df

def generate_mutation_classes():
    """Generates SBS mutation classes."""
    mutation_classes = []
    basic_mutations = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    sequence_contexts = ['A', 'C', 'G', 'T']
    for basic_mutation in basic_mutations:
        for five in sequence_contexts:
            for three in sequence_contexts:
                mutation_classes.append(five+basic_mutation+three)
    return mutation_classes

# classify SBS mutations coarse grained
mutation_classes = generate_mutation_classes()
mutation_types = {
    'C>A': 'C>A', 'C>G': 'C>G', 'C>T': 'C>T',
    'T>A': 'T>A', 'T>C': 'T>C', 'T>G': 'T>G'}        

def classify_mutations(five_prime, ref_base, alt_base, three_prime):
    """Classification of SBS into 96 features."""
    valid_bases = 'ACGT'
    if '-' in (five_prime, ref_base, alt_base, three_prime):
        return 'Other' 
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if ref_base in 'GA' and all(base in valid_bases for base in [five_prime, ref_base, alt_base, three_prime]):
        # Use reverse complement for the trinucleotide context
        rev_five_prime = complement[three_prime]
        rev_three_prime = complement[five_prime]
        rev_ref_base = complement[ref_base]
        rev_alt_base = complement[alt_base]
        mutation_class = rev_three_prime + rev_ref_base + '>' + rev_alt_base + rev_five_prime
    elif all(base in valid_bases for base in [five_prime, ref_base, alt_base, three_prime]):
        mutation_class = five_prime + ref_base + '>' + alt_base + three_prime 
    else:
        return 'Other'  
    return mutation_class 

# sequence stored in numerical format
def decode_sequence(encoded_seq):
    # Define a mapping from numbers to nucleotides
    num_to_nucleotide = {0: '-', 1: 'A', 2: 'T', 3: 'C', 4: 'G'}

    # Decode each sequence
    decoded_seqs = []
    for seq in encoded_seq:
        decoded_seq = ''.join([num_to_nucleotide.get(num, 'N') for num in seq[:, 0]])
        decoded_seqs.append(decoded_seq)
    return decoded_seqs 

def reverse_complement(mutation):
    # Dictionary to map each nucleotide to its complement
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # Split the mutation into the context (original) and the change (mutated)
    mut = mutation[1:4]
    ref = complement[mut[0]]
    alt = complement[mut[2]]
    f_base = complement[mutation[-1]]  # Middle base that is mutated
    t_base = complement[mutation[0]] # New base post mutation

    return f_base + ref + '>' + alt + t_base


def find_microhomology_bool(deletion_seq, prefix, suffix):
    """Method checks if there is a micohomology between mutation and sequence context."""
    deletion_seq = deletion_seq.replace("-", "")  # Remove gaps from deletion sequence
    max_length = min(len(deletion_seq), len(prefix), len(suffix))
    longest_beginning = ""
    longest_end = ""
    
    # Check for longest beginning match with suffix
    for i in range(1, max_length + 1):
        if deletion_seq[:i] == suffix[:i]:
            longest_beginning = deletion_seq[:i]
        else:
            break  

    # Check for longest end match with prefix
    for i in range(1, max_length + 1):
        if deletion_seq[-i:] == prefix[-i:]:
            longest_end = deletion_seq[-i:]
        else:
            break  
    print(longest_beginning)
    print(longest_end)

    # Ensure the match is partial and not the entire deletion sequence
    if (longest_beginning and len(longest_beginning) < len(deletion_seq)) or (longest_end and len(longest_end) < len(deletion_seq)):
        if longest_beginning:
            print(f"Longest microhomology at the beginning: {longest_beginning} (Length: {len(longest_beginning)})")
        if longest_end:
            print(f"Longest microhomology at the end: {longest_end} (Length: {len(longest_end)})")
        return True
    else:
        return False
    
def find_microhomology_int(deletion_seq, prefix, suffix):
    """Method returns how long the micohomology between mutation and sequence context is."""
    deletion_seq = deletion_seq.replace("-", "")  # Remove gaps from deletion sequence
    max_length = min(len(deletion_seq), len(prefix), len(suffix))
    longest_beginning = ""
    longest_end = ""
    
    # Check for longest beginning match with suffix
    for i in range(1, max_length + 1):
        if deletion_seq[:i] == suffix[:i]:
            longest_beginning = deletion_seq[:i]
        else:
            break  # Break as soon as there is no match

    # Check for longest end match with prefix
    for i in range(1, max_length + 1):
        if deletion_seq[-i:] == prefix[-i:]:
            longest_end = deletion_seq[-i:]
        else:
            break  # Break as soon as there is no match
    print(max(len(longest_beginning),len(longest_end)))
    return max(len(longest_beginning),len(longest_end))


def classify_indels(five,ref,alt,three, num_mut):
    """Classifying indels into coarse grained indel features (16)."""
    indel_mutation_dict = {"1DelC":0,"1DelT":0,"1InsC":0,"1InsT":0,"2DelRep":0,"3DelRep":0,"4DelRep":0,"5+DelRep":0,"2InsRep":0,"3InsRep":0,"4InsRep":0,"5+InsRep":0,"2MH":0,"3MH":0,"4MH":0,"5+MH":0}
    for i in range(len(five.tolist())):
        row_five = five.tolist()[i]
        row_three = three.tolist()[i]
        row_ref = ref.tolist()[i] 
        row_alt = alt.tolist()[i] 
        if set(row_alt) == {'-'}: #Deletion case
            if row_ref[1] == '-': #Deletion of length 1
                if row_ref[0] == "G" or "C":
                    indel_mutation_dict["1DelC"] +=1
                elif row_ref[0] == "T" or "A":
                    indel_mutation_dict["1DelT"] +=1
            elif row_ref[2] == '-':               #Deletion of length 2
                if find_microhomology_bool(row_ref,row_five,row_three):
                    indel_mutation_dict["2MH"] +=1
                else: 
                    indel_mutation_dict["2DelRep"] +=1
            elif row_ref[3] == '-':               #Deletion of length 3
                if find_microhomology_bool(row_ref,row_five,row_three):
                    indel_mutation_dict["3MH"] +=1
                else:
                    indel_mutation_dict["3DelRep"] +=1
            elif row_ref[4] == '-':               #Deletion of length 4
                if find_microhomology_bool(row_ref,row_five,row_three):
                    indel_mutation_dict["4MH"] +=1
                else:
                    indel_mutation_dict["4DelRep"] +=1
            elif any(row_ref[j] == '-' for j in range(5, len(row_ref))) or '-' not in row_ref:    #Deletion of length 5+
                if find_microhomology_bool(row_ref,row_five,row_three):
                    indel_mutation_dict["5+MH"] +=1
                else:
                    indel_mutation_dict["5+DelRep"] +=1

        elif set(row_ref) == {'-'}: # Insertion case
            if row_alt[1] == '-': #Insertion of length 1
                if row_alt[0] == "G" or "C":
                    indel_mutation_dict["1InsC"] +=1
                elif row_alt[0] == "T" or "A":
                    indel_mutation_dict["1InsT"] +=1
            elif row_alt[2] == '-': #Insertion of length 2
                indel_mutation_dict["2InsRep"] +=1
            elif row_alt[3] == '-': #Insertion of length 3
                indel_mutation_dict["3InsRep"] +=1
            elif row_alt[4] == '-': #Insertion of length 4
                indel_mutation_dict["4InsRep"] +=1
            elif any(row_alt[j] == '-' for j in range(5, len(row_alt))) or '-' not in row_alt: #Insertion of length 5
                indel_mutation_dict["5+InsRep"] +=1
    for key in indel_mutation_dict.keys():
        indel_mutation_dict[key] =  indel_mutation_dict[key] / (num_mut)
    return indel_mutation_dict

def classify_indels_detailed(five,ref,alt,three, num_mut):
    """Classifying indels into detailed indel features (83)."""
    indel_mutation_dict = {"1DelC1":0,"1DelC2":0,"1DelC3":0,"1DelC4":0,"1DelC5":0,"1DelC6":0,"1DelT1":0,"1DelT2":0,"1DelT3":0,"1DelT4":0,"1DelT5":0,"1DelT6":0,"1InsC0":0,"1InsC1":0,"1InsC2":0,"1InsC3":0,"1InsC4":0,"1InsC5":0,"1InsT0":0,"1InsT1":0,"1InsT2":0,"1InsT3":0,"1InsT4":0,"1InsT5":0,"2DelRep1":0,"2DelRep2":0,"2DelRep3":0,"2DelRep4":0,"2DelRep5":0,"2DelRep6":0,"3DelRep1":0,"3DelRep2":0,"3DelRep3":0,"3DelRep4":0,"3DelRep5":0,"3DelRep6":0,"4DelRep1":0,"4DelRep2":0,"4DelRep3":0,"4DelRep4":0,"4DelRep5":0,"4DelRep6":0,"5DelRep1":0,"5DelRep2":0,"5DelRep3":0,"5DelRep4":0,"5DelRep5":0,"5DelRep6":0,"2InsRep0":0,"2InsRep1":0,"2InsRep2":0,"2InsRep3":0,"2InsRep4":0,"2InsRep5":0,"3InsRep0":0,"3InsRep1":0,"3InsRep2":0,"3InsRep3":0,"3InsRep4":0,"3InsRep5+":0,"4InsRep0":0,"4InsRep1":0,"4InsRep2":0,"4InsRep3":0,"4InsRep4":0,"4InsRep5":0,"5InsRep0":0,"5InsRep1":0,"5InsRep2":0,"5InsRep3":0,"5InsRep4":0,"5InsRep5":0,"2MH1":0,"3MH1":0,"3MH2":0,"4MH1":0,"4MH2":0,"4MH3":0,"5+MH1":0,"5+MH2":0,"5+MH3":0,"5+MH4":0,"5+MH5":0}
    for i in range(len(five.tolist())):
        row_five = five.tolist()[i]
        row_three = three.tolist()[i]
        row_ref = ref.tolist()[i] 
        row_alt = alt.tolist()[i] 
        if set(row_alt) == {'-'}: #Deletion case
            if row_ref[1] == '-': #Deletion of length 1
                if row_ref[0] == "G" or "C":
                    classi = "1DelC"
                    lenC = len([char for char in row_three if char == 'C' and row_three.startswith('C' * (row_three.index(char) + 1))])
                    lenG = len([char for char in row_three if char == 'G' and row_three.startswith('G' * (row_three.index(char) + 1))])
                    lenadd = max(lenC,lenG) + 1
                    if lenadd > 6 :
                        lenadd = 6
                    indel_mutation_dict[classi+str(lenadd)] +=1
                elif row_ref[0] == "T" or "A":
                    classi = "1DelT"
                    lenT = len([char for char in row_three if char == 'T' and row_three.startswith('T' * (row_three.index(char) + 1))])
                    lenA = len([char for char in row_three if char == 'A' and row_three.startswith('A' * (row_three.index(char) + 1))])
                    lenadd = max(lenT,lenA) + 1
                    if lenadd > 6 :
                        lenadd = 6
                    indel_mutation_dict[classi+str(lenadd)] +=1
            elif row_ref[2] == '-':               #Deletion of length 2
                if find_microhomology_bool(row_ref,row_five,row_three):
                    indel_mutation_dict["2MH1"] +=1
                else: 
                    classi = "2DelRep"
                    main_string = row_three
                    sub_string = row_ref[:2]
                    count = sum(main_string[i:i+len(sub_string)] == sub_string for i in range(0, len(main_string) - len(sub_string) + 1, len(sub_string)))
                    if count > 5 :
                        count = 5
                    indel_mutation_dict[classi+str(count+1)] +=1
            elif row_ref[3] == '-':               #Deletion of length 3
                if find_microhomology_bool(row_ref,row_five,row_three):
                    indel_mutation_dict["3MH"+str(find_microhomology_int(row_ref,row_five,row_three))] +=1
                else:
                    classi = "3DelRep"
                    main_string = row_three
                    sub_string = row_ref[:3]
                    count = sum(main_string[i:i+len(sub_string)] == sub_string for i in range(0, len(main_string) - len(sub_string) + 1, len(sub_string)))
                    if count > 5:
                        count = 5
                    indel_mutation_dict[classi+str(count+1)] +=1
            elif row_ref[4] == '-':               #Deletion of length 4
                if find_microhomology_bool(row_ref,row_five,row_three):
                    indel_mutation_dict["4MH"+str(find_microhomology_int(row_ref,row_five,row_three))] +=1
                else:
                    classi = "4DelRep"
                    main_string = row_three
                    sub_string = row_ref[:4]
                    count = sum(main_string[i:i+len(sub_string)] == sub_string for i in range(0, len(main_string) - len(sub_string) + 1, len(sub_string)))
                    if count > 5 :
                        count = 5
                    indel_mutation_dict[classi+str(count+1)] +=1
            elif any(row_ref[j] == '-' for j in range(5, len(row_ref))) or '-' not in row_ref:    #Deletion of length 5+
                if find_microhomology_bool(row_ref,row_five,row_three):
                    if find_microhomology_int(row_ref,row_five,row_three)>5:
                        indel_mutation_dict["5+MH5"] +=1
                    else:
                        indel_mutation_dict["5+MH"+str(find_microhomology_int(row_ref,row_five,row_three))] +=1
                else:
                    classi = "5DelRep"
                    main_string = row_three
                    sub_string = row_ref[:5]
                    count = sum(main_string[i:i+len(sub_string)] == sub_string for i in range(0, len(main_string) - len(sub_string) + 1, len(sub_string)))
                    if count > 5 :
                        count = 5
                    indel_mutation_dict[classi+str(count+1)] +=1

        elif set(row_ref) == {'-'}: # Insertion case
            if row_alt[1] == '-': #Insertion of length 1
                if row_alt[0] == "G" or "C":
                    classi = "1InsC"
                    lenC = len([char for char in row_three if char == 'C' and row_three.startswith('C' * (row_three.index(char) + 1))])
                    lenG = len([char for char in row_three if char == 'G' and row_three.startswith('G' * (row_three.index(char) + 1))])
                    lenadd = max(lenC,lenG)
                    if lenadd > 5:
                        lenadd = 5
                    indel_mutation_dict[classi+str(lenadd)] +=1
                elif row_alt[0] == "T" or "A":
                    classi = "1InsT"
                    lenT = len([char for char in row_three if char == 'T' and row_three.startswith('T' * (row_three.index(char) + 1))])
                    lenA = len([char for char in row_three if char == 'A' and row_three.startswith('A' * (row_three.index(char) + 1))])
                    lenadd = max(lenT,lenA)
                    if lenadd > 5:
                        lenadd = 5
                    indel_mutation_dict[classi+str(lenadd)] +=1
            elif row_alt[2] == '-': #Insertion of length 2
                    classi = "2InsRep"
                    main_string = row_three
                    sub_string = row_ref[:2]
                    count = sum(main_string[i:i+len(sub_string)] == sub_string for i in range(0, len(main_string) - len(sub_string) + 1, len(sub_string)))
                    if count > 5 :
                        count = 5
                    indel_mutation_dict[classi+str(count)] +=1
            elif row_alt[3] == '-': #Insertion of length 3
                    classi = "3InsRep"
                    main_string = row_three
                    sub_string = row_ref[:3]
                    count = sum(main_string[i:i+len(sub_string)] == sub_string for i in range(0, len(main_string) - len(sub_string) + 1, len(sub_string)))
                    if count > 5 :
                        count = 5
                    indel_mutation_dict[classi+str(count)] +=1
            elif row_alt[4] == '-': #Insertion of length 4
                    classi = "4InsRep"
                    main_string = row_three
                    sub_string = row_ref[:4]
                    count = sum(main_string[i:i+len(sub_string)] == sub_string for i in range(0, len(main_string) - len(sub_string) + 1, len(sub_string)))
                    if count > 5 :
                        count = 5
                    indel_mutation_dict[classi+str(count)] +=1
            elif any(row_alt[j] == '-' for j in range(5, len(row_alt))) or '-' not in row_alt: #Insertion of length 5
                    classi = "5InsRep"
                    main_string = row_three
                    sub_string = row_ref[:5]
                    count = sum(main_string[i:i+len(sub_string)] == sub_string for i in range(0, len(main_string) - len(sub_string) + 1, len(sub_string)))
                    if count > 5 :
                        count = 5
                    indel_mutation_dict[classi+str(count)] +=1
    for key in indel_mutation_dict.keys():
        indel_mutation_dict[key] =  indel_mutation_dict[key] / (num_mut)
    return indel_mutation_dict