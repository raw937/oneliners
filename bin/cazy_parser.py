#!/usr/bin/env python3

#AUTHOR: Dr. Richard Allen White III
#REVISED: April 2022

# import libraries
import argparse
import csv
import os
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Aggregate Cazy18 annotations.')
    parser.add_argument('cazy18_file_paths', metavar='CAZY18', type=str, nargs='+', help='a Cazy18 file')
    parser.add_argument('-o', '--output', dest='output', default='.', metavar='DIR', type=str, help='the output directory')

    args = parser.parse_args()

    keys_by_cazy18_file_path_and_entry = {}
    keys_by_cazy18_file_path_and_gene_id = {}
    keys_by_cazy18_file_path_and_protein_id = {}

    for cazy18_file_path in args.cazy18_file_paths:
        keys_by_cazy18_file_path_and_entry[cazy18_file_path] = {}
        keys_by_cazy18_file_path_and_gene_id[cazy18_file_path] = {}
        keys_by_cazy18_file_path_and_protein_id[cazy18_file_path] = {}

        with open(cazy18_file_path, mode='r') as cazy18_file:
            for line in cazy18_file:
                key, *entries = re.split('\\s+', line.rstrip())

                for entry in entries:
                    if entry not in keys_by_cazy18_file_path_and_entry[cazy18_file_path]:
                        keys_by_cazy18_file_path_and_entry[cazy18_file_path][entry] = []

                    keys_by_cazy18_file_path_and_entry[cazy18_file_path][entry].append(key)

                    protein_id, *gene_ids = re.split(re.escape('|'), entry)

                    if protein_id not in keys_by_cazy18_file_path_and_protein_id[cazy18_file_path]:
                        keys_by_cazy18_file_path_and_protein_id[cazy18_file_path][protein_id] = []

                    keys_by_cazy18_file_path_and_protein_id[cazy18_file_path][protein_id].append(key)

                    for gene_id in gene_ids:
                        if gene_id not in keys_by_cazy18_file_path_and_gene_id[cazy18_file_path]:
                            keys_by_cazy18_file_path_and_gene_id[cazy18_file_path][gene_id] = []

                        keys_by_cazy18_file_path_and_gene_id[cazy18_file_path][gene_id].append(key)

    entry_out = {}

    for cazy18_file_path, keys_by_entry in keys_by_cazy18_file_path_and_entry.items():
        for entry, keys in keys_by_entry.items():
            if entry not in entry_out:
                entry_out[entry] = {}

            if cazy18_file_path not in entry_out[entry]:
                entry_out[entry][cazy18_file_path] = 0

            entry_out[entry][cazy18_file_path] += len(keys)

    gene_id_out = {}

    for cazy18_file_path, keys_by_gene_id in keys_by_cazy18_file_path_and_gene_id.items():
        for gene_id, keys in keys_by_gene_id.items():
            if gene_id not in gene_id_out:
                gene_id_out[gene_id] = {}

            if cazy18_file_path not in gene_id_out[gene_id]:
                gene_id_out[gene_id][cazy18_file_path] = 0

            gene_id_out[gene_id][cazy18_file_path] += len(keys)

    protein_id_out = {}

    for cazy18_file_path, keys_by_protein_id in keys_by_cazy18_file_path_and_protein_id.items():
        for protein_id, keys in keys_by_protein_id.items():
            if protein_id not in protein_id_out:
                protein_id_out[protein_id] = {}

            if cazy18_file_path not in protein_id_out[protein_id]:
                protein_id_out[protein_id][cazy18_file_path] = 0

            protein_id_out[protein_id][cazy18_file_path] += len(keys)

    def write_csv(name, data, header=''):
        with open(os.path.join(args.output, '{0}.csv'.format(name)), mode='w') as file:
            writer = csv.DictWriter(file, [header] + list(args.cazy18_file_paths), quoting=csv.QUOTE_NONNUMERIC)

            writer.writeheader()

            for key, orig_row in sorted(data.items(), key=lambda x: x[0]):
                new_row = orig_row.copy()
                new_row[header] = key

                for cazy18_file_path in args.cazy18_file_paths:
                    if cazy18_file_path not in new_row:
                        new_row[cazy18_file_path] = 0

                writer.writerow(new_row)

        return

    # def write_json(name, data):
    #     with open(os.path.join(args.output, '{0}.json'.format(name)), mode='w') as file:
    #         json.dump(data, file)
    #
    #     return

    for name, data in {
        'entry_counts': entry_out,
        'gene_id_counts': gene_id_out,
        'protein_id_counts': protein_id_out,
    }.items():
        write_csv(name, data)
        # write_json(name, data)

    pass
