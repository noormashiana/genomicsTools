import argparse
import numpy as np
import pandas as pd
from kmodes.kmodes import KModes
import matplotlib.pyplot as plt

# Takes in  Jessen .all. file

def parse_arguments():
    parser = argparse.ArgumentParser(description='KModes clustering to identify groups of chromosome vectors. Takes in sorted Jessen .all. file')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input TSV file.')
    parser.add_argument('-s', '--min_species', type=int, required=True, help='Minimum number of species.')
    parser.add_argument('-c', '--min_count', type=int, required=True, help='Minimum count.')
    parser.add_argument('-k', '--num_clusters', type=int, required=True, help='Number of clusters for KModes.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to output the binned TSV.')
    parser.add_argument('-E', '--do_elbow', action='store_true', help='Flag to perform elbow method for determining optimal clusters.')

    return parser.parse_args()

def load_and_filter_data(input_file, min_species, min_count):
    data = pd.read_csv(input_file, delimiter="\t")
    data = data.fillna("")
    data["species_count"] = data[data.columns[2:]].apply(lambda x: (x != "").sum(), axis=1)
    data = data[(data["species_count"] >= min_species) & (data["Count"] >= min_count)]
    return data

def expand_data(data):
    data_expanded = data.loc[data.index.repeat(data['Count'])]
    data_expanded = data_expanded.reset_index(drop=True).set_index('Cluster')
    data_expanded = data_expanded.drop(columns=['Count', 'species_count'])
    return data_expanded.apply(lambda row: row.replace('', row.name), axis=1)

def kmodes_clustering(data, num_clusters):
    km = KModes(n_clusters=num_clusters, init='Cao', n_init=1, random_state=1, verbose=1)
    return km.fit_predict(data)

def kmodes_cost(data, num_clusters):
    return KModes(n_clusters=num_clusters, init='Cao', n_init=1, random_state=1, verbose=0).fit(data).cost_


def merge_clusters(data, data_expanded, clusters):
    cluster = pd.DataFrame(clusters, columns=['Group'], index=data_expanded.index).reset_index()
    cluster = cluster.drop_duplicates(keep='first')
    return data.merge(cluster, on='Cluster', how='left')

def write_output(data, output_file):
    with open(output_file, 'w') as file:

        grouped = data.groupby('Group')
        first_group = True

        for name, group in grouped:
            group = group.drop(columns=["species_count", "Group"])
            group = group.sort_values(by='Count', ascending=False)

            if first_group:
                header = '\t'.join(group.columns)
                file.write(header + '\n')
                first_group = False

            file.write(f'##group={name}\n')

            group.to_csv(file, sep='\t', index=False, header=False)

def main():
    args = parse_arguments()

    data = load_and_filter_data(args.input, args.min_species, args.min_count)

    # expand to rows on count replace missing with ClusterID so that missing isn't counted as match
    data_expanded = expand_data(data)

    if args.do_elbow:
        # if I feel like later can add function to plot costs for elbow visualization
        costs = []
        for k in range(1, args.num_clusters + 1):
            costs.append(kmodes_cost(data_expanded, k))
            print(f"K = {k}: Cost = {costs[k - 1]}, dCost = {costs[k - 2] - costs[k - 1] if k > 1 else 'N/A'}")
        print("\nElbow Plot (pick K as lowest int of parrelel segment):")
        for i, cost in enumerate(costs):
            print(f"K = {i + 1}: {'*' * int(cost / max(costs) * 50)}")

    else:
        clusters = kmodes_clustering(data_expanded, args.num_clusters)

        data = merge_clusters(data, data_expanded, clusters)

        write_output(data, args.output)


if __name__ == "__main__":
    main()