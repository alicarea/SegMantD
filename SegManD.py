#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Jun 23 2016 

"""
DESCRIPTION OF PROGRAM
"""
import re


def prep_data(input_file):
    # read in Chromosome 22 and create pseudo dataframe
    with open(input_file, "r") as ifile:
        file_content = ifile.read()
        file_content = re.findall(r"#+ (Chromosome [0-9XY]+) #+\n([^#]*)", file_content, flags=re.DOTALL)

    output = {}
    for chrom, content in file_content:
        if chrom != "Chromosome 22":
            continue
        output[chrom] = []
        content = content.split("\n")[1:-1]
        for indx, line in enumerate(content):
            line = line.split("\t")
            line[0] = tuple(map(int, line[0].split("_")[1:3]))
            line[1] = line[1].split(";")
            line[1] = [tuple(map(int, loc.split("_"))) for loc in line[1]]
            for segment in line[1]:
                output[chrom].append((line[0], segment))
            content[indx] = line

        output[chrom] = sorted(output[chrom], key=lambda x: x[1][1], reverse=True)
        output[chrom] = sorted(output[chrom], key=lambda x: x[1][0])
        output[chrom] = sorted(output[chrom], key=lambda x: x[0][1], reverse=True)
        output[chrom] = sorted(output[chrom], key=lambda x: x[0][0])

    return output


def strip_internals(input_data):
    """
    Loop through dataframe and remove any cases that are completely contained by another
    :param input_data:
    :return:
    """
    output_data = []
    while input_data:
        output_data.append(input_data[0])
        p_bac, p_chrom = input_data[0]
        del input_data[0]
        delete_indicies = []
        for indx, line in enumerate(input_data):
            n_bac, n_chrom = line
            # first ensure upcoming line is a subset of the previous line (within the BAC)
            if n_bac[0] >= p_bac[0] and n_bac[1] <= p_bac[1]:
                # then check within the chromosomal sequences.
                if n_chrom[0] >= p_chrom[0] and n_chrom[1] <= p_chrom[1]:
                    delete_indicies.append(indx)
                    break
            elif n_bac[0] < p_bac[1]:
                continue
            else:
                # The postitions are sorted, so no need to keep checking once passed the previous subset
                break
        if delete_indicies:
            delete_indicies = sorted(delete_indicies, reverse=True)
            for indx in delete_indicies:
                del input_data[indx]
    return output_data


def extend_singles(input_data):
    output_data = []
    while input_data:
        output_data.append(input_data[0])
        p_bac, p_chrom = input_data[0]
        del input_data[0]
        delete_indicies = []
        for indx, line in enumerate(input_data):
            n_bac, n_chrom = line
            if n_bac[0] > p_bac[1]:
                break
            elif n_chrom[0] > p_chrom[1] or n_chrom[1] < p_chrom[0]:
                continue
            else:
                p_bac = (p_bac[0], max([p_bac[1], n_bac[1]]))  # NOTE: No need to check p_bac[0] because of sort order
                p_chrom = (min([p_chrom[0], n_chrom[0]]), p_chrom[1])
                p_chrom = (p_chrom[0], max([p_chrom[1], n_chrom[1]]))
                delete_indicies.append(indx)

        output_data[-1] = (p_bac, p_chrom)
        if delete_indicies:
            delete_indicies = sorted(delete_indicies, reverse=True)
            for indx in delete_indicies:
                del input_data[indx]
    output_data = sorted(output_data, key=lambda x: (x[0][0], x[1][0]))
    return output_data

if __name__ == '__main__':
    data = prep_data("parsed_chrAL_segDup_5kb.txt")
    for chromosome in data:
        cleaned_singles = strip_internals(data[chromosome])
        cleaned_singles = extend_singles(cleaned_singles)  # I'm not certain this will always complete in one round, need to think more
        out_string = "#### %s ####\n" % chromosome
        for data_line in cleaned_singles:
            out_string += "al591856_%s_%s_bp\t%s_%s\n" % (data_line[0][0], data_line[0][1],
                                                          data_line[1][0], data_line[1][1])
        print(out_string)

