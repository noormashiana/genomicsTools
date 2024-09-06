#!/usr/bin/env python

import pandas as pd
import argparse

def ortho_to_confusion(table_in, table_out):
    ortho_table = pd.read_csv(table_in, delimiter="\t")

    def parse_row(row, f_out):
        orthogroup = row['Orthogroup']
        
        for column_name in ortho_table.columns[1:]:
            species_name = column_name.strip()

            if pd.isna(row[species_name]):
                continue
            
            genes = row[species_name].split(',')
            
            for gene in genes:
                line = f"{gene.strip()}\t{orthogroup}\t{species_name}\n"
                f_out.write(line)

    with open(table_out, 'w') as f_out:
        f_out.write("Gene\tHmmHogHit\tSpecies\n")
        ortho_table.apply(lambda row: parse_row(row, f_out), axis=1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert orthogroup table to confusion matrix format.')
    parser.add_argument('input', help='Path to the input file')
    parser.add_argument('output', help='Path to the output file')

    args = parser.parse_args()
    ortho_to_confusion(args.input, args.output)