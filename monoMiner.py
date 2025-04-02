import readline  # Import readline to enable autocomplete
import os
import glob
import subprocess
import re
import concurrent.futures
from functools import partial
import shutil

# Enable TAB autocomplete
readline.parse_and_bind("tab: complete")

def find_motifs(library_file, reference_sequence, similarity_threshold=0.6):
    motifs = []

    # Load the reference sequence
    with open(reference_sequence, 'r') as ref_file:
        ref_file.readline()
        reference_seq = ref_file.read().replace('\n', '')

    # Calculate maximum allowed divergence
    max_divergence = int(len(reference_seq) * (1 - similarity_threshold))

    # Process the library
    with open(library_file, 'r') as lib_file:
        current_sequence = ''
        for line in lib_file:
            if line.startswith('>'):
                if current_sequence:
                    i = 0
                    while i <= len(current_sequence) - len(reference_seq):
                        subsequence = current_sequence[i:i+len(reference_seq)]
                        divergence = sum(1 for j in range(len(subsequence)) if subsequence[j] != reference_seq[j])
                        if divergence <= max_divergence:
                            motifs.append(subsequence)
                            i += len(reference_seq)
                        else:
                            i += 1
                current_sequence = ''
            else:
                current_sequence += line.strip()

        # Check the last sequence
        if current_sequence:
            i = 0
            while i <= len(current_sequence) - len(reference_seq):
                subsequence = current_sequence[i:i+len(reference_seq)]
                divergence = sum(1 for j in range(len(subsequence)) if subsequence[j] != reference_seq[j])
                if divergence <= max_divergence:
                    motifs.append(subsequence)
                    i += len(reference_seq)
                else:
                    i += 1

    return motifs

def save_motifs_to_fasta(motifs, reference_name, library_name):
    library_base = os.path.splitext(os.path.basename(library_name))[0]
    library_prefix = library_base[:3]  # First 3 letters of library
    first_number = re.search(r'\d', library_base).group()
    
    output_file = f"{reference_name}_motifs_{library_base}.fasta"
    valid_nucleotides = {'A', 'T', 'C', 'G'}

    with open(output_file, 'w') as out_fasta:
        for i, motif in enumerate(motifs, start=1):
            if set(motif.upper()).issubset(valid_nucleotides):
                read_name = f"{library_prefix}{first_number}{i}"
                out_fasta.write(f">{read_name}\n")
                out_fasta.write(f"{motif}\n")
            else:
                print(f"Motif {i} from library {library_base} was discarded (invalid nucleotides).")

    return output_file

def concatenate_motifs(reference_name):
    motif_files = glob.glob(f"{reference_name}_motifs_*.fasta")
    if not motif_files:
        print("No motif files found to concatenate.")
        return False

    with open("final.fasta", 'w') as final_fasta:
        for motif_file in motif_files:
            with open(motif_file, 'r') as mf:
                final_fasta.write(mf.read())
    print("All motif files have been concatenated into 'final.fasta'.")
    return True

def run_cd_hit_filter_size(min_copies):
    try:
        cd_hit_script = shutil.which("cd_hit_filter_size.py")
        if not cd_hit_script:
            print("Error: 'cd_hit_filter_size.py' not found in PATH.")
            return False

        print(f"Executing '{cd_hit_script} final.fasta {min_copies} 1'...")
        subprocess.run([cd_hit_script, "final.fasta", str(min_copies), "1"], check=True)
        print("cd_hit_filter_size.py executed successfully.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing cd_hit_filter_size.py: {e}")
        return False

def completer(text, state):
    if not text:
        matches = glob.glob('*')
    else:
        matches = glob.glob(text + '*')
    try:
        return matches[state]
    except IndexError:
        return None

def setup_autocomplete():
    readline.set_completer(completer)
    readline.parse_and_bind('tab: complete')

def get_user_input(prompt, file_only=True):
    while True:
        try:
            user_input = input(prompt).strip()
            if file_only:
                user_input = os.path.expanduser(user_input)
                if os.path.isfile(user_input):
                    return user_input
                else:
                    print(f"File '{user_input}' not found. Please try again.")
            else:
                return user_input
        except KeyboardInterrupt:
            print("\nOperation canceled by user.")
            exit()
        except EOFError:
            print("\nOperation terminated.")
            exit()

def read_mapping_file(mapping_file):
    mapping = {}
    with open(mapping_file, 'r', encoding='utf-8') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                key, value = parts[0], parts[1]
                mapping[key] = value
    return mapping

def extract_species_code(seq_id):
    match = re.match(r'([a-zA-Z]+)', seq_id)
    if match:
        return match.group(1)
    else:
        return None

def process_fasta(fasta_file, mapping):
    sequences = []
    with open(fasta_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip()
                species_code = extract_species_code(seq_id)
                if species_code and species_code in mapping:
                    species_name = mapping[species_code]
                else:
                    species_name = "Unknown"
                sequences.append((seq_id, species_name))
    return sequences

def write_output(output_file, sequences):
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("Sequence\tSpecies\n")
        for seq, species in sequences:
            f.write(f"{seq}\t{species}\n")
    print(f"\nOutput file '{output_file}' created successfully.")

def process_final_fasta(mapping, min_copies):
    final_fasta = f"final.fasta.nr0.{min_copies}.sel.fasta"
    if not os.path.isfile(final_fasta):
        print(f"Error: File {final_fasta} not found. Ensure 'cd_hit_filter_size.py' ran correctly.")
        return

    sequences = process_fasta(final_fasta, mapping)
    output_file = "output.tsv"
    write_output(output_file, sequences)

def main():
    print("=== NO MORE MANUAL PROCESSING! AUTOMATED MOTIF ANALYSIS PIPELINE ===\n")

    # Set up autocomplete
    setup_autocomplete()

    mapping_file = get_user_input("Enter mapping file name (e.g., amex for Astyanax mexicanus): ")
    mapping = read_mapping_file(mapping_file)
    print("File loaded successfully.\n")

    reference_sequence = get_user_input("Enter reference sequence file path (FASTA format): ")

    if not os.path.isfile(reference_sequence):
        print(f"Error: File {reference_sequence} does not exist.")
        return

    min_copies = get_user_input("For CD-hit, remove sequences with how many copies or fewer? ", file_only=False)

    reference_name = os.path.basename(reference_sequence).split('.')[0]

    library_files = glob.glob("*.fq")
    if not library_files:
        print("No .fq files found in current directory.")
        return

    print(f"\nFound {len(library_files)} .fq files to process.\n")

    process_func = partial(find_motifs, reference_sequence=reference_sequence)
    motifs_files = []
    max_threads = 50

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_threads) as executor:
        futures = {executor.submit(process_func, lib_file): lib_file for lib_file in library_files}
        for future in concurrent.futures.as_completed(futures):
            lib_file = futures[future]
            try:
                motifs = future.result()
                if motifs:
                    print(f"{len(motifs)} motifs found in {lib_file}.")
                    saved_file = save_motifs_to_fasta(motifs, reference_name, lib_file)
                    motifs_files.append(saved_file)
                    print(f"Motifs saved to {saved_file}")
                else:
                    print(f"No motifs found in {lib_file}.")
            except Exception as exc:
                print(f"File {lib_file} generated an exception: {exc}")

        print("\nConcatenating all motif files...")

    if concatenate_motifs(reference_name):
        if run_cd_hit_filter_size(min_copies):
            process_final_fasta(mapping, min_copies)
        else:
            print("Error running cd_hit_filter_size.py. Verify script availability.")
    else:
        print("No motif files concatenated. Execution halted.")

if __name__ == "__main__":
    main()
