#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def confusion_to_ortho(table_in, table_out):
    hog_table = pd.read_csv(table_in, delimiter="\t", cols=["Gene", "Species"] )

    with np.errstate(divide='ignore', invalid='ignore'):
        hog_table["log_E"] = np.log10(hog_table["E-value"].astype(float))

    hog_table = hog_table[hog_table["log_E"] <= -23]

    gene_groups = hog_table.groupby(['HmmHogHit', 'Species'])['Gene'].apply(lambda x: ','.join(x))

    pivot_df = gene_groups.unstack(fill_value='')

    pivot_df.reset_index(inplace=True)

    pivot_df.to_csv(table_out, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert hmmscore confusion matrix to orthofinder format.')
    parser.add_argument('input', help='Path to the input file')
    parser.add_argument('output', help='Path to the output file')

    args = parser.parse_args()
    confusion_to_ortho(args.input, args.output)