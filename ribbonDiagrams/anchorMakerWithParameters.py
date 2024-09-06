#!/usr/bin/env python

import os
import sys
import logging
import numpy as np
import pandas as pd
import argparse
import yaml

from OrthofinderOutputToConfusionMatrix import ortho_to_confusion

'''
TODO: 
    1) Make an order_by flag so pick a species's chromosomes to order by, use like idk figure something out
    2) Make a color_by so pick a species "Species", take Species.bed, get number of chromosomes as n, find n unique colors to make up color_dict
    3) Maybe combine 1 and 2

'''

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def make_seqids(bed_path, species_list_in_order, min_genes_in_chromosome, cluster_table, output_directory):
    logging.info("Creating seqids file...")
    # Right now, largest to smallest ordered, need to implement min crossing algorithm

    seqids = []

    for species in species_list_in_order:
        bed_file = os.path.join(bed_path, f"{species}.bed")
        bed = pd.read_csv(bed_file, delimiter="\t", header=None, usecols=[0, 3], names=["Chrom", "Gene"])

        #Filter chromosomes based on the group size
        chrom_counts = bed["Chrom"].value_counts()
        filtered_chroms = chrom_counts[chrom_counts >= min_genes_in_chromosome].index.tolist()

        #Merge HOG table with bed data to get the Chrom column
        gene_connections = cluster_table[cluster_table["Species"] == species]
        gene_connections = pd.merge(gene_connections, bed, on="Gene", how="left")

        # #Check if 'Chrom' column exists after merge
        # if 'Chrom' not in gene_connections.columns:
        #     logging.warning(f"'Chrom' column not found for species {species}")
        #     continue

        # Order chromosomes based on gene connections in the HOG table
        chrom_order = gene_connections.groupby("Chrom_x").size().reindex(filtered_chroms).sort_values(ascending=False).index.tolist()

        seqids.append(",".join(chrom_order))

    seqids_file = os.path.join(output_directory, "seqids")
    with open(seqids_file, 'w') as f:
        for seqid in seqids:
            f.write(f"{seqid}\n")

    logging.info(f"Created {seqids_file}")

def make_seqids2(species_list_in_order, cluster_table, output_directory, min_genes_in_chromosome):
    def order_chromosomes(merged_data, min_genes_in_chromosome):
        chrom_order = merged_data.groupby(['Chrom_one', 'Chrom_two'])['Counts'].sum().reset_index()
        chrom_order = chrom_order[chrom_order['Counts'] >= min_genes_in_chromosome]
        chrom_order = chrom_order.sort_values(by='Counts', ascending=False)
        chrom_one_order = chrom_order['Chrom_one'].unique().tolist()
        chrom_two_order = chrom_order['Chrom_two'].unique().tolist()
        return chrom_one_order, chrom_two_order

    chromosome_orders = []
    for i in range(len(species_list_in_order) - 1):
        species_one = species_list_in_order[i]
        species_two = species_list_in_order[i+1]
        species_one_data = cluster_table[cluster_table["Species"] == species_one]
        species_two_data = cluster_table[cluster_table["Species"] == species_two]
        merged_data = pd.merge(species_one_data, species_two_data, on="HmmHogHit", suffixes=('_one', '_two'))
        merged_data['Counts'] = merged_data.groupby(['Chrom_one', 'Chrom_two'])['HmmHogHit'].transform('size')
        chrom_one_order, chrom_two_order = order_chromosomes(merged_data, min_genes_in_chromosome)
        
        # Optimize ordering based on previous pair if available
        if i > 0:
            prev_chrom_two_order = chromosome_orders[-1][1]
            # Reorder chrom_one_order based on the previous chrom_two_order
            chrom_one_order = sorted(chrom_one_order, key=lambda x: prev_chrom_two_order.index(x) if x in prev_chrom_two_order else len(prev_chrom_two_order))
        
        chromosome_orders.append([chrom_one_order, chrom_two_order])

    seqids = []
    for i, species in enumerate(species_list_in_order):
        if i == 0:
            chrom_order = chromosome_orders[0][0]  # For the first species
        else:
            chrom_order = chromosome_orders[i-1][1]  # For subsequent species
        seqids.append(",".join(chrom_order))

    seqids_file = os.path.join(output_directory, "seqids")
    with open(seqids_file, 'w') as f:
        for seqid in seqids:
            f.write(f"{seqid}\n")
    logging.info(f"Created {seqids_file}")

def make_anchors(color_dict, species_list_in_order, cluster_table, bed_path, e_cutoff, gene_count_max, synteny_file, 
                 post_prob, strand_count_min, output_directory, min_genes_in_chromosome):
    
    logging.info("Creating anchors...")
    cluster_table = pd.read_csv(cluster_table, delimiter="\t")

    if e_cutoff:
        with np.errstate(divide='ignore', invalid='ignore'):
            cluster_table["log_E"] = np.log10(cluster_table["E-value"].astype(float))
        cluster_table = cluster_table[cluster_table["log_E"] <= e_cutoff] # PARAMETER: E-23 Cutoff

    # add chromosome of gene to cluster_table
    concatenated_bed = pd.DataFrame()
    for species in species_list_in_order:
        bed_file = os.path.join(bed_path, f"{species}.bed")
        bed = pd.read_csv(bed_file, delimiter="\t", usecols=[0, 3], names=["Chrom", "Gene"])
        concatenated_bed = pd.concat([concatenated_bed, bed])
    cluster_table = pd.merge(cluster_table, concatenated_bed, on="Gene", how="left")

    # Min 20 genes in hog to get drawn - try getting rid of this now that there is the chrom-chrom req
    if "Species" not in cluster_table.columns:
        cluster_table["Species"] = cluster_table["Gene"].apply(lambda gene: gene.split("|")[0][:3])
    gene_counts = cluster_table.groupby("HmmHogHit")["Species"].count()
    #print(gene_counts.sort_values(ascending=False).head())
    if gene_count_max:
        gene_counts = gene_counts[gene_counts <= gene_count_max] # PARAMETER: gene count cutoff
    cluster_table = cluster_table[cluster_table["HmmHogHit"].isin(gene_counts.index)]

    # only grassfams in synteny with PostProb >= .99 will be drawn
    synteny = pd.read_csv(synteny_file, delimiter="\t")
    synteny.dropna(axis=0, inplace=True)
    if post_prob:
        synteny = synteny[synteny["PostProb"] >= post_prob] # PARAMETER: Posterior Probability
    else:
        synteny.columns = ["HmmHogHit", "PostID"]
    synteny_dict = synteny.set_index("HmmHogHit")["PostID"].to_dict()
    cluster_table = cluster_table[cluster_table["HmmHogHit"].isin(synteny["HmmHogHit"])]

    #make_seqids2(species_list_in_order, cluster_table, output_directory, min_genes_in_chromosome)
    make_seqids(bed_path, species_list_in_order, min_genes_in_chromosome, cluster_table, output_directory)

    # here can set minimum strand count between chromosomes
    def addAnchor(row, f_out):
        if row['Counts'] >= strand_count_min: # PARAMETER: Cutoff for Chrom-Chrom strand min/max
            block = synteny_dict.get(row["HmmHogHit"])
            color = color_dict[block]
            
            line = f"{color}*{row['Gene_one']}\t{row['Gene_one']}\t{row['Gene_two']}\t{row['Gene_two']}\t1\t+\n"
            f_out.write(line)
            return 1
        return 0
        
    layout_file = os.path.join(output_directory, "layout")
    layout_species = 0
    
    with open(layout_file, 'w') as layout_out:
        layout_out.write("# y, xstart, xend, rotation, color, label, va,  bed\n")
        y_positions = np.linspace(0.95, 0.05, len(species_list_in_order))
        
        for i, species in enumerate(species_list_in_order):
            va = "bottom" if i == len(species_list_in_order) - 1 else "top"
            bed_out = os.path.join(bed_path, f"{species}.bed")
            layout_out.write(f"{y_positions[i]:.2f}, .15, .95, 0, black, {species}, {va}, {bed_out}\n")
            layout_species += 1
        
        layout_out.write("# edges\n")

        #make n-1 anchors.simple files
        for i in range(len(species_list_in_order)-1):
            species_one = species_list_in_order[i]
            species_two = species_list_in_order[i+1]
            anchorfile = os.path.join(output_directory, f"{species_one}.{species_two}.anchors.simple")

            species_one_data = cluster_table[cluster_table["Species"] == species_one]
            species_two_data = cluster_table[cluster_table["Species"] == species_two]
            
            merged_data = pd.merge(species_one_data, species_two_data, on="HmmHogHit", suffixes=('_one', '_two'))
            merged_data['Counts'] = merged_data.groupby(['Chrom_one', 'Chrom_two'])['HmmHogHit'].transform('size')

            with open(anchorfile, 'w') as f_out:
                written = merged_data.apply(lambda row: addAnchor(row, f_out), axis=1).sum()
            logging.info(f"Created {anchorfile} with {written} entries.") # this logic can be used to generate layout file
            
            layout_out.write(f"e, {i}, {i + 1}, {anchorfile}\n")

        logging.info(f"Created {layout_file} with {layout_species} species.")

def main(config_path):
    try:
        with open(config_path, "r") as file:
            config = yaml.safe_load(file)
    except Exception as e:
        logging.error(f"Error reading config file: {e}")
        sys.exit(1)
    
    output_directory = config.get("output_directory")
    if not output_directory:
        logging.error("Output directory is not specified in the config file.")
        sys.exit(1)
    
    os.makedirs(output_directory, exist_ok=True)
    logging.info(f"Output directory: {output_directory}")

    cluster_format = config.get("cluster_table_format")
    if cluster_format == "Orthofinder":
        confusion = os.path.join(output_directory, "confusion_matrix.tsv")
        ortho_to_confusion(config["cluster_table"], confusion)
        #make_seqids(config["bed_path"], config["species_list_in_order"], config["min_genes_in_chromosome"], confusion, output_directory)
        make_anchors(
            color_dict=config["color_dict"],
            species_list_in_order=config["species_list_in_order"],
            cluster_table=confusion,
            bed_path=config["bed_path"],
            e_cutoff=config["e_cutoff"],
            gene_count_max=config["gene_count_max"],
            synteny_file=config["synteny_file"],
            post_prob=config["post_prob"],
            strand_count_min=config["strand_count_min"],
            output_directory=output_directory,
            min_genes_in_chromosome=config["min_genes_in_chromosome"]
        )
    elif cluster_format == "Confusion_Matrix":
        #make_seqids(config["bed_path"], config["species_list_in_order"], config["min_genes_in_chromosome"], config["cluster_table"], output_directory)
        make_anchors(
            color_dict=config["color_dict"],
            species_list_in_order=config["species_list_in_order"],
            cluster_table=config["cluster_table"],
            bed_path=config["bed_path"],
            e_cutoff=config["e_cutoff"],
            gene_count_max=config["gene_count_max"],
            synteny_file=config["synteny_file"],
            post_prob=config["post_prob"],
            strand_count_min=config["strand_count_min"],
            output_directory=output_directory,
            min_genes_in_chromosome=config["min_genes_in_chromosome"]
        )
    else:
        logging.error("Invalid cluster_table_format. Must be 'Orthofinder' or 'Confusion_Matrix'")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate anchors for synteny analysis")
    parser.add_argument("config", type=str, help="Path to the configuration file yaml format")
    args = parser.parse_args()

    main(args.config)