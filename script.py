mport numpy as np
import math
import random
from Bio import Entrez, SeqIO


def genbank_sample(entrez_email, search_term, sample_count=50):
    """Queries GenBank to create a realistic gene pool."""
    Entrez.email = entrez_email
    gene_sample = []

    print("Fetching sequences from GenBank...")
    search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=sample_count * 2)
    search_results = Entrez.read(search_handle)
    id_list = search_results["IdList"]
    search_handle.close()

    for access_id in id_list:
        try:
            with Entrez.efetch(db="nucleotide", id=access_id, rettype="gb", retmode="text") as handle:
                record = SeqIO.read(handle, "genbank")
                for feature in record.features:
                    if feature.type == "CDS":
                        seq = str(feature.location.extract(record).seq).upper()
                        if len(seq) % 3 == 0 and seq.startswith("ATG"):
                            gc = (seq.count('G') + seq.count('C')) / len(seq)
                            if 0.35 <= gc <= 0.65:
                                gene_sample.append(seq)
                                break 
        except Exception as e:
            continue

    if not gene_sample:
        return []
        
    chromosome_sim = random.choices(gene_sample, k=10)
    return chromosome_sim


def mutate_sequence(chromosome, mutation_rate=0.01):
    """Applies ONLY SNP mutations to keep sequence lengths strictly identical for Tajima's D."""
    homolog_seq = list("".join(chromosome)) # Converting to a list makes index replacement much faster
    bases = ['A', 'T', 'C', 'G']
    
    total_mutations = int(len(homolog_seq) * mutation_rate)
    
    for _ in range(total_mutations):
        idx = random.randint(0, len(homolog_seq) - 1)
        current_base = homolog_seq[idx]
        new_base = random.choice([b for b in bases if b != current_base])
        homolog_seq[idx] = new_base
        
    return "".join(homolog_seq)

def generate_population_fasta(entrez_email, search_term, pop_size=5, output_file="simulated_population.fasta"):
    """Generates a base chromosome and mutates it for N individuals, saving to FASTA."""
    base_chromosome = genbank_sample(entrez_email, search_term)
    
    if not base_chromosome:
        print("Failed to generate base chromosome.")
        return

    print(f"Generating population of {pop_size} individuals...")
    
    try:
        with open(output_file, "w") as f:
            for i in range(pop_size):
                # We use a lower mutation rate so sequences can still be aligned easily
                mutated_seq = mutate_sequence(base_chromosome, mutation_rate=0.005) 
                
                f.write(f">Individual_{i+1}\n")
                # Wrap sequence to 80 characters per line standard
                for j in range(0, len(mutated_seq), 80):
                    f.write(mutated_seq[j:j+80] + "\n")
                    
        print(f"Success! Saved population to {output_file}")
    except Exception as e:
        print(f"Error writing file: {e}")


email = "email" 
query = "Escherichia coli[Orgn] AND complete genome[Title]"
generate_population_fasta(email, query, pop_size=5, output_file="simulated_population.fasta")

def parse_fasta(input_file):
    sequences = {}
    current_id = None
    current_seq = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue 
            
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]  
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = "".join(current_seq)

    return sequences

def calculate_k2p_distance(seq1, seq2):
    """Calculates Kimura 2-Parameter distance between two sequences."""
    transitions = 0
    transversions = 0
    valid_sites = 0
    
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    
    for a, b in zip(seq1, seq2):
        if a == b:
            valid_sites += 1
        elif a in 'AGCT' and b in 'AGCT':
            valid_sites += 1
            if (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
                transitions += 1
            else:
                transversions += 1
                
    if valid_sites == 0: return 0
    
    P = transitions / valid_sites
    Q = transversions / valid_sites
    
    try:
        # K2P formula for substitutions per site
        d = -0.5 * math.log(1 - 2*P - Q) - 0.25 * math.log(1 - 2*Q)
        # Scale back to total expected differences across the sequence length
        return d * valid_sites 
    except ValueError:
        # Fallback to standard Hamming distance if math domain error occurs (extreme divergence)
        return transitions + transversions

def tajimas_d_analysis(sequences):
    seqs = list(sequences.values())
    n = len(seqs) 
    if n < 3: 
        print("Need at least 3 sequences")
        return
    
    # TRUNCATE to the shortest sequence to ensure column iteration works without IndexErrors
    min_len = min(len(s) for s in seqs)
    seqs = [s[:min_len] for s in seqs]
    seq_length = min_len
    
    # Calculate Segregating Sites (S)
    S = 0
    for i in range(seq_length):
        column = [s[i] for s in seqs]
        if len(set(column)) > 1: 
            S += 1
            
    # Tajima's Coefficients
    a1 = sum(1.0/i for i in range(1, n))
    a2 = sum(1.0/(i**2) for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))
    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1**2))
    e1 = c1 / a1
    e2 = c2 / ((a1**2) + a2)
    
    # Watterson's Estimator (M)
    M = S / a1
    
    # Pairwise Differences (Standard & Modified)
    std_pair_diffs = 0
    k2p_pair_diffs = 0
    num_pairs = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            # Standard Hamming Distance
            diffs = sum(1 for a, b in zip(seqs[i], seqs[j]) if a != b)
            std_pair_diffs += diffs
            
            # Modified K2P Distance
            k2p_diffs = calculate_k2p_distance(seqs[i], seqs[j])
            k2p_pair_diffs += k2p_diffs
            
            num_pairs += 1
            
    std_avg_polymorphism = std_pair_diffs / num_pairs
    k2p_avg_polymorphism = k2p_pair_diffs / num_pairs
    
    # Variances & Final Calculations
    variance = (e1 * S) + (e2 * S * (S - 1))
    
    if variance == 0:
        print("Variance is zero; Tajima's D cannot be calculated.")
        return
        
    std_d_val = std_avg_polymorphism - M
    std_tajimas_d = std_d_val / math.sqrt(variance)
    
    k2p_d_val = k2p_avg_polymorphism - M
    k2p_tajimas_d = k2p_d_val / math.sqrt(variance)
    
    #Output
    print("\n--- Results ---")
    print(f"Segregating Sites (S): {S}")
    print(f"Watterson's Estimator (M): {M:.4f}\n")
    
    print("STANDARD TAJIMA's D:")
    print(f"Average Pairwise Differences: {std_avg_polymorphism:.4f}")
    print(f"Standard Tajima's D: {std_tajimas_d:.4f}\n")
    
    print("K2P MODIFIED TAJIMA's D:")
    print(f"K2P Average Pairwise Differences: {k2p_avg_polymorphism:.4f}")
    print(f"Modified Tajima's D: {k2p_tajimas_d:.4f}\n")
    
    print("--- Comparison ---")
    if k2p_tajimas_d > std_tajimas_d:
        print("The K2P modified Tajima's D is slightly higher/more positive than the standard value.")
        print("This implies the standard method was underestimating the true number of historical mutations (due to saturation/multiple hits at the same site).")
    else:
        print("The K2P value is similar or slightly altered, indicating sequence divergence is low enough that saturation is not heavily skewing the results.")


input_seqs = parse_fasta("simulated_population.fasta")
tajimas_d_analysis(input_seqs)
