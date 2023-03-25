#! /usr/bin/env python
__author__ = 'Jianjun Jin'

import os
import sys
from platform import system
if system() == "Windows":
    line_br = "\r\n"
elif system() == "Darwin":
    line_br = "\r"
else:
    line_br = "\n"
try:
    from Bio import SeqIO, SeqFeature
except ImportError:
    sys.stdout.write("Python package biopython not found!" + line_br +
                     "You could use \"pip install biopython\" to install it." + line_br)
    sys.exit()
from optparse import OptionParser
from glob import glob


"""Convert Genbank format file to tbl format (https://www.ncbi.nlm.nih.gov/Sequin/table.html)"""
# example1 from Seq2 of https://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html


def get_options(description=""):
    usage = "This is a GetOrganelle script for converting genbank format to tbl format, \n" \
            "which can be further used to submit through Banklt (https://www.ncbi.nlm.nih.gov/WebSub/)\n" \
            "Usage: gb_to_tbl.py gb_files"
    parser = OptionParser(description=description, usage=usage)
    parser.add_option("-o", dest="output",
                      help="Output directory. Default: along with the original file.")
    parser.add_option("-t", dest="gene_types", default="CDS,tRNA,rRNA,gene,repeat_region,source",
                      help="Annotation type taken as gene. Set 'all' to report all types. Default: %default")
    parser.add_option("-q", dest="qualifiers", default="gene,product,note,standard_name,rpt_type",
                      help="The qualifiers to record. Set 'all' to report all qualifiers. Default: %default.")
    # parser.add_option("--ignore-format-error", dest="ignore_format_error", default=False, action="store_true",
    #                   help="Skip the Error: key \"*\" not found in annotation. Not suggested.")
    options, argv = parser.parse_args()
    if not len(argv):
        parser.print_help()
        sys.exit()
    if system() == "Windows":
        new_argv = []
        for input_fn_pattern in argv:
            new_argv.extend(glob(input_fn_pattern))
        argv = new_argv
    return options, argv


def parse_bio_gb_locations(location_feature):
    if type(location_feature) == SeqFeature.CompoundLocation:
        return [parse_bio_gb_locations(location)[0] for location in location_feature.parts]
    elif type(location_feature) == SeqFeature.FeatureLocation:
        return [(int(location_feature.start), int(location_feature.end), location_feature.strand)]
    else:
        raise ValueError(str(type(location_feature)))


def location_feature_to_str(location_feature, feature_type):
    locations = parse_bio_gb_locations(location_feature=location_feature)
    # switch location[0][0] and locations[0][1] if strand/location[0][2] says reverse/-1
    # locations[0][0] + 1 because Location records indices in Biopython
    # example1: 2626  2590    tRNA
    lines = ["\t".join([str(locations[0][0] + 1), str(locations[0][1])][::locations[0][2]]) +
             "\t" + feature_type + line_br]
    # add more parts if location_feature is a CompoundLocation
    # example1: 2570  2535
    for loc in locations[1:]:
        lines.append("\t".join([str(loc[0] + 1), str(loc[1])][::loc[2]]) + line_br)
    return "".join(lines)


def genbank_to_tbl(genbank_f, out_base, accepted_type_dict, qualifiers_dict):
    this_records = list(SeqIO.parse(genbank_f, "genbank"))
    if accepted_type_dict is None:
        record_all_types = True
        accepted_type_dict = {}
    else:
        record_all_types = False
    if qualifiers_dict is None:
        record_all_qualifiers = True
        qualifiers_dict = {}
    else:
        record_all_qualifiers = False
    with open(out_base + ".fasta", "w") as output_fs:
        pass
    with open(out_base + ".tbl", "w") as output_tbl:
        pass
    for go_record, seq_record in enumerate(this_records):
        this_seq_name = (str(go_record + 1) + "--") * int(bool(len(this_records) > 1)) + \
                        seq_record.name * int(bool(seq_record.name))
        if not this_seq_name:
            raise Exception("Sequence name not found in the " + str(go_record + 1) + " record of " + genbank_f)
        else:
            with open(out_base + ".fasta", "a") as output_fs:
                output_fs.write(">" + this_seq_name + line_br + str(seq_record.seq) + line_br)
            with open(out_base + ".tbl", "a") as output_tbl:
                output_tbl.write(">Feature " + this_seq_name + line_br)
                for feature in seq_record.features:
                    if record_all_types or feature.type.lower() in accepted_type_dict:
                        this_type = accepted_type_dict.get(feature.type.lower(), feature.type)
                        try:
                            # some locations are compound locations
                            location_str = location_feature_to_str(feature.location, feature_type=this_type)
                        except ValueError as e:
                            sys.stdout.write("Warning: abnormal location " + str(e) +
                                             line_br + str(feature) +
                                             "\nin the " +
                                             str(go_record + 1) + " record of " + genbank_f + " .. skipped!" + line_br)
                            continue
                        else:
                            output_tbl.write(location_str)
                            for qualifier_k in feature.qualifiers:
                                if record_all_qualifiers or qualifier_k.lower() in qualifiers_dict:
                                    this_qualifier = qualifiers_dict.get(qualifier_k.lower(), qualifier_k)
                                    this_values = feature.qualifiers[qualifier_k]
                                    for this_val in this_values:
                                        output_tbl.write("\t\t\t" + this_qualifier + "\t" + this_val + "\n")


def main():
    options, argv = get_options(
        "Convert Genbank format file to tbl format (https://www.ncbi.nlm.nih.gov/Sequin/table.html)" + line_br +
        "By jinjianjun@mail.kib.ac.cn")
    if options.gene_types != "all":
        accepted_types = {this_type.lower(): this_type for this_type in options.gene_types.split(",")}
    else:
        accepted_types = None
    if options.qualifiers != "all":
        qualifiers = {this_q.lower(): this_q for this_q in options.qualifiers.split(",")}
    else:
        qualifiers = None

    # check output file names
    output_base_names = []
    if options.output:
        if not os.path.exists(options.output):
            os.mkdir(options.output)
        elif os.path.isfile(options.output):
            raise FileExistsError(options.output + " is a file!")
    else:
        options.output = ""
    for gb_file in argv:
        this_out_f = os.path.join(options.output, os.path.basename(gb_file))
        if this_out_f.endswith(".gb"):
            this_out_f = this_out_f[:-3]
        elif this_out_f.endswith(".gbk"):
            this_out_f = this_out_f[:-4]
        elif this_out_f.endswith(".genbank"):
            this_out_f = this_out_f[:-8]
        output_base_names.append(this_out_f)

    # converting
    for go_gb, gb_file in enumerate(argv):
        if os.path.isfile(gb_file):
            genbank_to_tbl(gb_file, output_base_names[go_gb],
                           accepted_type_dict=accepted_types, qualifiers_dict=qualifiers)
        else:
            sys.stdout.write("Error: " + gb_file + " not found!" + line_br)


if __name__ == '__main__':
    main()
