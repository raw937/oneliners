#!/usr/bin/env python3

import argparse
import csv
# import json
import os
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Aggregate KO classifications in KEGG annotations.')
    parser.add_argument('classification_file_paths', metavar='CLASSIFICATION', type=str, nargs=1, help='a KO classification file')
    parser.add_argument('annotation_file_paths', metavar='ANNOTATION', type=str, nargs='+', help='a KEGG annotation file')
    parser.add_argument('-o', '--output', dest='output', default='.', metavar='DIR', type=str, help='the output directory')

    args = parser.parse_args()

    keys_by_annotation_file_path_and_entry = {}

    for annotation_file_path in args.annotation_file_paths:
        keys_by_annotation_file_path_and_entry[annotation_file_path] = {}

        with open(annotation_file_path, mode='r') as annotation_file:
            for line in annotation_file:
                key, *entries = re.split('\\s+', line.strip())

                for entry in entries:
                    if entry not in keys_by_annotation_file_path_and_entry[annotation_file_path]:
                        keys_by_annotation_file_path_and_entry[annotation_file_path][entry] = []

                    keys_by_annotation_file_path_and_entry[annotation_file_path][entry].append(key)

    data = {}

    previous_depth = None
    previous_parts = []

    for classification_file_path in args.classification_file_paths:
        with open(classification_file_path, mode='r') as classification_file:
            for line in classification_file:
                orig_parts = re.split('\\t', line.rstrip())

                depth = 0

                for part in orig_parts:
                    if part == '':
                        depth += 1
                    else:
                        break

                if depth == 2:
                    orig_parts = orig_parts[0:2:] + ['\t'.join(orig_parts[2::])]

                orig_parts = (previous_parts[0:depth:] + orig_parts[depth::])[0:depth + 1:]

                my_data = data

                for index, part in enumerate(orig_parts):
                    if index == 3:
                        counts = {}

                        for annotation_file_path, keys_by_entry in keys_by_annotation_file_path_and_entry.items():
                            counts[annotation_file_path] = len(keys_by_entry[part]) if (part in keys_by_entry) else 0

                        my_data[part] = counts

                        break
                    else:
                        if part not in my_data:
                            my_data[part] = {}

                        my_data = my_data[part]

                previous_parts = orig_parts

    level_1_out = {}

    for level_1_key, level_1_data in data.items():
        if level_1_key not in level_1_out:
            level_1_out[level_1_key] = {}

            for annotation_file_path in args.annotation_file_paths:
                level_1_out[level_1_key][annotation_file_path] = 0

        for level_2_key, level_2_data in level_1_data.items():
            for level_3_key, level_3_data in level_2_data.items():
                for level_4_key, level_4_data in level_3_data.items():
                    for level_5_key, level_5_data in level_4_data.items():
                        level_1_out[level_1_key][level_5_key] += level_5_data

    level_2_out = {}

    for level_1_key, level_1_data in data.items():
        for level_2_key, level_2_data in level_1_data.items():
            if level_2_key not in level_2_out:
                level_2_out[level_2_key] = {}

                for annotation_file_path in args.annotation_file_paths:
                    level_2_out[level_2_key][annotation_file_path] = 0

            for level_3_key, level_3_data in level_2_data.items():
                for level_4_key, level_4_data in level_3_data.items():
                    for level_5_key, level_5_data in level_4_data.items():
                        level_2_out[level_2_key][level_5_key] += level_5_data

    level_3_out = {}

    for level_1_key, level_1_data in data.items():
        for level_2_key, level_2_data in level_1_data.items():
            for level_3_key, level_3_data in level_2_data.items():
                if level_3_key not in level_3_out:
                    level_3_out[level_3_key] = {}

                    for annotation_file_path in args.annotation_file_paths:
                        level_3_out[level_3_key][annotation_file_path] = 0

                for level_4_key, level_4_data in level_3_data.items():
                    for level_5_key, level_5_data in level_4_data.items():
                        level_3_out[level_3_key][level_5_key] += level_5_data

    level_4_out = {}

    for level_1_key, level_1_data in data.items():
        for level_2_key, level_2_data in level_1_data.items():
            for level_3_key, level_3_data in level_2_data.items():
                for level_4_key, level_4_data in level_3_data.items():
                    if level_4_key not in level_4_out:
                        level_4_out[level_4_key] = {}

                        for annotation_file_path in args.annotation_file_paths:
                            level_4_out[level_4_key][annotation_file_path] = 0

                    for level_5_key, level_5_data in level_4_data.items():
                        level_4_out[level_4_key][level_5_key] += level_5_data

    def write_csv(name, data, header=''):
        with open(os.path.join(args.output, '{0}.csv'.format(name)), mode='w') as file:
            writer = csv.DictWriter(file, [header] + list(args.annotation_file_paths), quoting=csv.QUOTE_NONNUMERIC)

            writer.writeheader()

            for key, orig_row in sorted(data.items(), key=lambda x: x[0]):
                new_row = orig_row.copy()
                new_row[header] = key
                writer.writerow(new_row)

        return

    # def write_json(name, data):
    #     with open(os.path.join(args.output, '{0}.json'.format(name)), mode='w') as file:
    #         json.dump(data, file)
    #
    #     return

    for name, data in {
        'level_1': level_1_out,
        'level_2': level_2_out,
        'level_3': level_3_out,
        'level_4': level_4_out,
    }.items():
        write_csv(name, data)
        # write_json(name, data)

    pass
