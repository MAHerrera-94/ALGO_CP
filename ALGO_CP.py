import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO, SeqRecord
from Bio.Seq import Seq
from collections import Counter
import random

# Set random seed
random_seed = random.randint(0, 10000)
np.random.seed(random_seed)
random.seed(random_seed)

# Load the Multiple Sequence Alignment (MSA) file in FASTA format
alignment = AlignIO.read("AcpP_input_MSA.fasta", "fasta")
sequences = [str(record.seq) for record in alignment]
num_positions = alignment.get_alignment_length()

# Dictionaries containing physicochemical properties of amino acids
vdW_volume = {
    'A': 88.6, 'C': 108.5, 'D': 111.1, 'E': 138.4, 'F': 189.9,
    'G': 60.1, 'H': 153.2, 'I': 166.7, 'K': 168.6, 'L': 166.7,
    'M': 162.9, 'N': 114.1, 'P': 112.7, 'Q': 143.8, 'R': 173.4,
    'S': 89.0, 'T': 116.1, 'V': 140.0, 'W': 227.8, 'Y': 193.6
}

hydropathy_values = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
    'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2
}

isoelectric_point = {
    'A': 6.00, 'R': 10.76, 'N': 5.41, 'D': 2.77, 'C': 5.07, 'Q': 5.65,
    'E': 3.22, 'G': 5.97, 'H': 7.59, 'I': 6.02, 'L': 5.98, 'K': 9.74,
    'M': 5.74, 'F': 5.48, 'P': 6.30, 'S': 5.68, 'T': 5.60, 'W': 5.89,
    'Y': 5.66, 'V': 5.96
}

# Function to calculate properties for each amino acid in a sequence
def calculate_properties(seq):
    properties = []
    for aa in seq:
        if aa == '-':
            properties.append([np.nan] * 3)
        else:
            hydropathy = hydropathy_values[aa]
            pi = isoelectric_point[aa]
            volume = vdW_volume.get(aa, np.nan)
            properties.append([hydropathy, pi, volume])
    return properties

# Function to average properties over aligned positions
def average_properties(properties):
    avg_properties = np.nanmean(properties, axis=0)
    return avg_properties

# Initialize storage for position-wise properties
position_properties = []

# Calculate and average properties for each position in the alignment
for pos in range(num_positions):
    column = [seq[pos] for seq in sequences]
    properties = calculate_properties(column)
    avg_properties = average_properties(properties)
    position_properties.append(avg_properties)

# Convert the list of position properties to a DataFrame for easier manipulation
properties_df = pd.DataFrame(position_properties, columns=["Hydropathy", "pI", "Volume"])

# Filter the alignment to find blocks of 50 or more residues without missing values
blocks = []
current_block = []
for pos in range(num_positions):
    if not np.isnan(properties_df.iloc[pos]).any():
        current_block.append(pos)
    else:
        if len(current_block) >= 50:
            blocks.append(current_block)
        current_block = []
if len(current_block) >= 50:
    blocks.append(current_block)

# Combine the filtered blocks into a final dataset
filtered_df = pd.concat([properties_df.iloc[block] for block in blocks])
filtered_positions = filtered_df.index



#### Route Ac #####
# Redefine position frequencies based on filtered positions
position_frequencies = []
for pos in range(num_positions):
    if pos in filtered_positions:
        column = [seq[pos] for seq in sequences if seq[pos] != '-']
        if len(column) >= 100:
            freq = Counter(column)
            position_frequencies.append(freq)
        else:
            position_frequencies.append(None)
    else:
        position_frequencies.append(None)

# Identify highly conserved positions in the filtered positions (default = 0.9 frquency)
conserved_positions = [pos for pos, freq in enumerate(position_frequencies) if freq and max(freq.values()) / sum(freq.values()) > 0.9]



#### Route Ap #####
# Define user-set weights for physicochemical properties
user_defined_weights = {
    "Hydropathy": 0.22,
    "pI": 0.44,
    "Volume": 0.33,
}


#### Weighting Coefficient r for Au ####
# Define weight_r (value between 0-1) for the creation of the unified probability distribution Au. Values closer to 1 favour Ap-guided sequence design.
weight_r = 0.50

#Do NOT adjust:
weight_Ac = (1.00 - weight_r)


#### Sequence Generation ####
# Generate new sequences based on the user-defined weights. Set the number of sequences to generate (default = 5000)
num_novel_sequences = 5000
novel_sequences = []

for _ in range(num_novel_sequences):
    new_sequence = ''
    for i in range(num_positions):
        if i not in filtered_positions:
            continue  # Skip positions that were not used in the training data
        if i in conserved_positions:
            new_sequence += sequences[0][i]  # Preserve strictly conserved positions relative to MSA
        elif position_frequencies[i] is not None:
            residues = list(position_frequencies[i].keys())
            frequencies = list(position_frequencies[i].values())
            avg_properties = position_properties[i]
            physicochemical_scores = []
            for aa in residues:
                property_scores = []
                for j, feature in enumerate(filtered_df.columns[:-1]):
                    if feature == "Hydropathy":
                        property_value = hydropathy_values[aa]
                    elif feature == "pI":
                        property_value = isoelectric_point[aa]
                    elif feature == "Volume":
                        property_value = vdW_volume[aa]

                    
                    difference = abs(avg_properties[j] - property_value)
                    property_score = user_defined_weights[feature] / (difference + 1e-6)  # Add a small constant to avoid division by zero
                    property_scores.append(property_score)
                total_property_score = sum(property_scores)
                physicochemical_scores.append(total_property_score)
            
            # Normalize physicochemical scores to generate Ap distribution:
            total_physicochemical_score = sum(physicochemical_scores)
            route_Ap = [score / total_physicochemical_score for score in physicochemical_scores]

            # Adjust Ap by weighting coefficient r:
            weighted_Ap = [score * weight_r for score in route_Ap]
            
            # Generate Ac distribution:
            total_frequency = sum(frequencies)
            route_Ac = [freq / total_frequency for freq in frequencies]

            # Adjust Ac by weighting coefficient (1-r):
            weighted_Ac = [score * weight_Ac for score in route_Ac]

            # Combine Ac and Ap:
            combined_scores = [weighted_Ap[j] + weighted_Ac[j] for j in range(len(residues))]

            # Normalize combined scores to generate unified probability profile Au:
            Au = sum(combined_scores)
            probabilities = [score / Au for score in combined_scores]

            # Randomly select an amino acid based on Au:
            new_residue = np.random.choice(residues, p=probabilities)
            new_sequence += new_residue
        else:
            new_sequence += '-'  # Placeholder for positions not meeting criteria
    novel_sequences.append(new_sequence)

# Remove gaps from generated sequences
novel_sequences = [seq.replace('-', '') for seq in novel_sequences]



##### Sequence Properties and Data Export #####

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import ProtParam

# Function to calculate the percentage identity
def calculate_percentage_identity(reference, sequence):
    matches = 0
    count = 0
    for a, b in zip(reference, sequence):
        if a == 'X' or b == 'X':
            continue
        if a == b:
            matches += 1
        count += 1
    return (matches / count) * 100 if count > 0 else 0


# Function to calculate GRAVY
def calculate_gravy(sequence):
    prot_param = ProtParam.ProteinAnalysis(sequence)
    return prot_param.gravy()

# User-specified reference sequence (EcAcpP used as default)
reference_sequence = "MSTIEERVKKIIGEQLGVKQEEVTNNASFVEDLGADSLDTVELVMALEEEFDTEIPDEEAEKITTVQAAIDYINGHQA"

# Calculate GRAVY and percentage ID for each novel sequence
gravies = []
percentage_ids = []

for seq in novel_sequences:
    gravies.append(calculate_gravy(seq))
    percentage_ids.append(calculate_percentage_identity(reference_sequence, seq))
    

# Create a DataFrame with sequences and their computed properties
sequences_df = pd.DataFrame({
    'Sequence': novel_sequences,
    'GRAVY': gravies,
    'Percentage ID': percentage_ids
})

# Write the sequences to a separate FASTA file
fasta_records_unranked = [SeqRecord(Seq(seq), id=f'seq_{i + 1}', description='ALGO-CP Sequence') for i, seq in enumerate(novel_sequences)]
SeqIO.write(fasta_records_unranked, "ALGO_CP_sequences.fasta", "fasta")

# Write the ranked sequences and their properties to a CSV file
sequences_df.to_csv("ALGO_CP_sequences.csv", index=False)

print("FASTA written.")
print("CSV written.")
print("Job Complete!")