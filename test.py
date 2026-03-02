# script to check if tsv files are the same

import pandas as pd

def compare_tsv_files(file1, file2):
    '''Compares two TSV files and prints out any differences.'''
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t')

    if df1.equals(df2):
        print("The TSV files are the same.")
    else:
        print("The TSV files are different.")
        # Optionally, you can print out the differences
        diff = pd.concat([df1, df2]).drop_duplicates(keep=False)
        print("Differences:")
        print(diff)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Compare two TSV files for differences.')
    parser.add_argument('file1', type=str, help='Path to the first TSV file.')
    parser.add_argument('file2', type=str, help='Path to the second TSV file.')
    args = parser.parse_args()

    compare_tsv_files(args.file1, args.file2)
    