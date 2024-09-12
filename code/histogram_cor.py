import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count

# Define file paths
input_file_path = 'matrix_full_sorghum.txt'
temp_dir = 'temp_outputs'  # Directory to store temporary output files

def process_chunk(args):
    start_line, end_line, chunk_index, gene_names = args
    lower_tri_values = []

    with open(input_file_path, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')

        # Skip the header and lines before the start_line
        if start_line > 0:
            for _ in range(start_line + 1):  # +1 to skip the header
                next(reader)
        else:
            next(reader)  # Skip the header

        for row_index, row in enumerate(reader, start=start_line):
            if row_index >= end_line:
                break
            values = row[1:]   # Values for the current row

            # Iterate over columns up to the current row index to cover the lower triangle
            for col_index in range(min(row_index, len(values))):
                value_str = values[col_index]
                try:
                    value = float(value_str)
                    lower_tri_values.append(value)
                except ValueError:
                    continue  # Skip invalid entries

    return lower_tri_values

def merge_values(results):
    all_values = []
    for result in results:
        all_values.extend(result)
    return np.array(all_values)

def plot_histogram(values):
    plt.hist(values, bins=20, color='blue')
    plt.title("")
    plt.xlabel("Values")
    plt.ylabel("Frequency")
    plt.savefig("histogram_sorghum_cor_.png")
    plt.close()

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
    num_chunks = 100
    chunk_size = num_lines // num_chunks
    chunk_ranges = [(i * chunk_size, (i + 1) * chunk_size if i < num_chunks - 1 else num_lines, i, gene_names) for i in range(num_chunks)]

    # Create temp directory if it doesn't exist
    os.makedirs(temp_dir, exist_ok=True)

    # Process chunks in parallel
    with Pool(processes=num_chunks) as pool:
        results = pool.map(process_chunk, chunk_ranges)

    # Merge all the collected lower triangular values
    lower_tri_values = merge_values(results)

    # Plot the histogram of all lower triangular values
    plot_histogram(lower_tri_values)

    # Clean up temporary files
    os.rmdir(temp_dir)

if __name__ == "__main__":
    main()

