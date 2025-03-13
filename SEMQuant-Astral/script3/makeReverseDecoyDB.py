import sys

def reverse_fasta_sequences(input_file, output_file):
    with open(input_file, 'r') as fasta, open(output_file, 'w') as output:
        sequence = ''
        header = ''
        for line in fasta:
            if line.startswith('>'):
                if header:
                    output.write(header + '\n' + sequence[::-1] + '\n')
                header = '>Decoy_' + line[1:].strip()
                sequence = ''
            else:
                sequence += line.strip()
        if header:  # write the last sequence
            output.write(header + '\n' + sequence[::-1] + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)
    reverse_fasta_sequences(sys.argv[1], sys.argv[2])