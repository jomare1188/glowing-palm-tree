# Define file paths
input_file_path = 'matrix_full_cane.txt'
output_file_path = 'triplets_p60_cane_fancy.tsv'
temp_dir = 'temp_outputs'  # Directory to store temporary output files

# Set the threshold
threshold = 0.6

def process_chunk(args):
    start_line, end_line, chunk_index, gene_names = args
    temp_output_file = os.path.join(temp_dir, f'temp_output_{chunk_index}.tsv')

    with open(input_file_path, 'r') as infile, open(temp_output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Skip the header and lines before the start_line
        if start_line > 0:
            for _ in range(start_line + 1):  # +1 to skip the header
                next(reader)
        else:
            next(reader)  # Skip the header

        for row_index, row in enumerate(reader, start=start_line):
            if row_index >= end_line:
                break
            row_gene = row[0]  # Gene name for the current row
            values = row[1:]   # Values for the current row

            # Iterate over columns up to the current row index to cover the lower triangle
            for col_index in range(row_index):
                value_str = values[col_index]
                try:
                    value = float(value_str)
                    if abs(value) >= threshold:
                        source_gene = row_gene
                        target_gene = gene_names[col_index]
                        writer.writerow([source_gene, target_gene, value])
                except ValueError:
                    # Handle non-numeric values if necessary
                    continue  # Skip invalid entries

def merge_files(num_chunks):
    with open(output_file_path, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')

        for chunk_index in range(num_chunks):
            temp_output_file = os.path.join(temp_dir, f'temp_output_{chunk_index}.tsv')
            with open(temp_output_file, 'r') as infile:
                reader = csv.reader(infile, delimiter='\t')
                for row in reader:
                    writer.writerow(row)

def main():
    # Read the header to get gene names
    with open(input_file_path, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        header = next(reader)
        gene_names = header[1:]  # Exclude the first entry which is assumed to be an empty placeholder

    # Determine the number of lines in the file
    with open(input_file_path, 'r') as infile:
        num_lines = sum(1 for _ in infile) - 1  # Exclude header

    # Set up parallel processing
    num_chunks = 50
    chunk_size = num_lines // num_chunks
    chunk_ranges = [(i * chunk_size, (i + 1) * chunk_size if i < num_chunks - 1 else num_lines, i, gene_names) for i in range(num_chunks)]

    # Create temp directory if it doesn't exist
    os.makedirs(temp_dir, exist_ok=True)

    # Process chunks in parallel
    with Pool(processes=num_chunks) as pool:
        pool.map(process_chunk, chunk_ranges)

    # Merge temporary files into the final output file
    merge_files(num_chunks)

    # Clean up temporary files
    for chunk_index in range(num_chunks):
        os.remove(os.path.join(temp_dir, f'temp_output_{chunk_index}.tsv'))
    os.rmdir(temp_dir)

if __name__ == "__main__":
    main()
