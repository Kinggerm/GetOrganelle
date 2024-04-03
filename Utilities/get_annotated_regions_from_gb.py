#!/usr/bin/env python

import os
import time
import sys
try:
    from Bio import SeqIO, SeqFeature
except ImportError:
    sys.stdout.write("Python package biopython not found!\n"
                     "You could use \"pip install biopython\" to install it.\n")
    sys.exit()
from optparse import OptionParser
from platform import system
from glob import glob
# moved from https://github.com/Kinggerm/PersonalUtilities

# Copyright(C) 2022 Jianjun Jin


major_version, minor_version = sys.version_info[:2]
if major_version == 2 and minor_version >= 7:
    python_version = "2.7+"
elif major_version == 3 and minor_version >= 5:
    python_version = "3.5+"
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)


def get_options():
    usage = "Usage: get_annotated_regions_from_gb.py gb_files -o out_dir"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", dest="out_put",
                      help="Output.")
    parser.add_option("-t", dest="gene_types", default="CDS,tRNA,rRNA",
                      help="Annotation type taken as gene. Default: CDS,tRNA,rRNA")
    parser.add_option("--separate-copy", dest="one_copy", default=True, action="store_false",
                      help="By default, only keep one copy (see '--copy-mode' for more) "
                           "if there are several regions with the same name. "
                           "Exception: if there are one copy with intron(s) and another copy without intron, "
                           "they would be both kept. This exception was specially made for the convenience of "
                           "commonly-incorrectly-annotated rps12 gene of plastome.")
    parser.add_option("--copy-mode", dest="copy_mode", default="leastN_longest",
                      help="first|longest|leastN|leastN_longest (default).")
    parser.add_option("--separate-exon", dest="combine_exon", default=True, action="store_false",
                      help="By default, combining exons.")
    parser.add_option("--keys", dest="gene_keys", default="gene,label,product,note",
                      help="The key to the gene name: gene, label, product or other keys in the qualifiers region."
                           "Default: %default.")
    parser.add_option("--mix", dest="mix", default=False, action="store_true",
                      help="Mix different genes into a single fasta file. "
                           "In this mode, the sequence header will be >gene_name - gb_info")
    parser.add_option("--case-mode", dest="case_treatment", default="first",
                      help="first: Gene name case-non-sensitive. Consistent to the first appearance. \n"
                           "lower: Gene name case-non-sensitive. All gene name set to lower case. \n"
                           "upper: Gene name case-non-sensitive. All gene name set to Upper case. \n"
                           "raw: Gene name case-sensitive. ")
    parser.add_option("--ignore-format-error", dest="ignore_format_error", default=False, action="store_true",
                      help="Skip the Error: key \"*\" not found in annotation. Not suggested.")
    parser.add_option("--translate-to-product", dest="product_to_gene", default=True, action="store_false",
                      help="Translate the tRNA gene name to the form of their product. Default: False")
    parser.add_option("--overwrite", dest="overwrite", default=False, action="store_true",
                      help="Choose to overwrite previous result.")
    options, argv = parser.parse_args()
    if not (options.out_put and bool(len(argv))):
        parser.print_help()
        sys.exit()
    if system() == "Windows":
        new_argv = []
        for input_fn_pattern in argv:
            new_argv.extend(glob(input_fn_pattern))
        argv = new_argv
    if options.copy_mode not in {"longest", "first", "leastN", "leastN_longest"}:
        sys.stdout.write("Error: invalid value " + options.copy_mode + " for '--copy-mode'!\n")
        sys.exit()
    if options.case_treatment not in {"first", "upper", "lower", "raw"}:
        sys.stdout.write("Error: invalid value " + options.case_treatment + " for '--case-mode'!\n")
        sys.exit()
    return options, argv


if python_version == "2.7+":
    # python2
    import string
    translator = string.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")


    def complementary_seq(input_seq):
        return string.translate(input_seq, translator)[::-1]

else:
    # python3
    translator = str.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")


    def complementary_seq(input_seq):
        return str.translate(input_seq, translator)[::-1]


def complementary_seqs(input_seq_iter):
    return tuple([complementary_seq(seq) for seq in input_seq_iter])


missing_base = {"N", "?", "n"}
head_ = "_+_"
tail_ = "_-_"


def count_n(seq):
    return len([base for base in seq if base in missing_base])


def parse_bio_gb_locations(location_feature):
    if type(location_feature) == SeqFeature.CompoundLocation:
        return [parse_bio_gb_locations(location)[0] for location in location_feature.parts]
    elif type(location_feature) == SeqFeature.FeatureLocation:
        return [(int(location_feature.start), int(location_feature.end), location_feature.strand)]
    else:
        raise ValueError(str(type(location_feature)))


def embed_in(candidate_small, candidate_large):
    small_start, small_end = candidate_small
    large_start, large_end = candidate_large
    # both circular
    if small_start >= small_end and large_start >= large_end:
        return small_end <= large_end
    elif small_start >= small_end:
        return False
    elif large_start >= large_end:
        return True
    else:
        return small_end <= large_end


trna_translate_table = {"Ala": "A",
                        "Arg": "R",
                        "Asn": "N",
                        "Asp": "D",
                        "Cys": "C",
                        "Gln": "Q",
                        "Glu": "E",
                        "Gly": "G",
                        "His": "H",
                        "Ile": "I",
                        "Leu": "L",
                        "Lys": "K",
                        "Met": "M",
                        "fMet": "fM",
                        "Phe": "F",
                        "Pro": "P",
                        "Ser": "S",
                        "Thr": "T",
                        "Trp": "W",
                        "Tyr": "Y",
                        "Val": "V"}


gene_name_lower_to_first = {}


def modify_gene_name(product_name, translate_product_to_gene_name, case_mode):
    # treat case
    if case_mode == "upper":
        product_name = product_name.upper()
    elif case_mode == "lower":
        product_name = product_name.lower()
    elif case_mode == "first":
        if product_name.lower() in gene_name_lower_to_first:
            product_name = gene_name_lower_to_first[product_name.lower()]
        else:
            gene_name_lower_to_first[product_name.lower()] = product_name
    elif case_mode == "raw":
        pass
    # treat product name
    if translate_product_to_gene_name:
        if product_name.startswith("tRNA-") or product_name.startswith("trna-"):
            short_name = product_name.replace("tRNA-", "").replace("trna-", "")
            if short_name[:3] in trna_translate_table:
                return "trn" + trna_translate_table[short_name[:3]] + "-" + short_name[3:].replace("(", "").replace(")", "")
            elif short_name[:4] in trna_translate_table:
                return "trn" + trna_translate_table[short_name[:4]] + "-" + short_name[4:].replace("(", "").replace(")", "")
            else:
                return product_name
        elif "rrna" in product_name or "rRNA" in product_name:
            res_name = "rrn" + product_name.replace("rrna", "").replace("rRNA", "").replace(" ", "").replace("_", "")
            if "S" in res_name:
                if res_name[res_name.index("S") - 1].isdigit():
                    return res_name.replace("S", "")
            return res_name
        else:
            return product_name
    else:
        return product_name


def get_seqs(seq_record, accepted_types, gene_keys,
             ignore_format_error=False, trans_product_to_gene=True, case_m="first"):
    original_seq = str(seq_record.seq)

    def get_seq_with_gb_loc(in_location):
        in_start, in_end, in_strand = in_location
        if in_start >= in_end:
            in_seq = original_seq[in_start:] + original_seq[:in_end]
        else:
            in_seq = original_seq[in_start: in_end]
        if in_strand == 1:
            return in_seq
        else:
            return complementary_seq(in_seq)
    gene_regions = []
    name_counter = {}
    taken_loc = set()
    for feature in seq_record.features:
        if feature.type in accepted_types:
            this_key_found = False
            location_error = False
            for gene_key in gene_keys:
                if gene_key in feature.qualifiers:
                    try:
                        locations = parse_bio_gb_locations(feature.location)
                    except ValueError as e:
                        location_error = True
                        sys.stdout.write("Warning: " + str(e) + "\n")
                        break
                    else:
                        this_name = [modify_gene_name(feature.qualifiers[gene_key][0], trans_product_to_gene, case_m),
                                     "", ""]
                        if this_name[0] not in name_counter:
                            name_counter[this_name[0]] = 1
                        else:
                            name_counter[this_name[0]] += 1
                            this_name[1] = "__copy" + str(name_counter[this_name[0]])
                        if len(locations) > 1:
                            for i, loc in enumerate(locations):
                                this_name[2] = "__exon" + str(i + 1)
                                if loc not in taken_loc:
                                    gene_regions.append(
                                        [tuple(this_name)] + list(loc) + [get_seq_with_gb_loc(loc), feature.type])
                                    taken_loc.add(loc)
                        else:
                            gene_regions.append(
                                [tuple(this_name)] + list(locations[0]) + [get_seq_with_gb_loc(locations[0]), feature.type])
                        this_key_found = True
                        break
            if not location_error and not this_key_found and not ignore_format_error:
                sys.stdout.write("\nError: ")
                sys.stdout.write("Present key(s) \"" + ",".join(gene_keys) + "\" not found in annotation:\n")
                sys.stdout.write(str(feature))
                raise NotImplementedError
    gene_regions.sort(key=lambda x: (x[1], -x[2], x[0]))
    intergenic_regions = []
    end_of_last_region = 0
    if len(gene_regions) == 1:
        if gene_regions[0][1] == gene_regions[0][2]:
            pass
        else:
            anchor1 = [gene_regions[0][0][0], gene_regions[0][0][2], tail_ if gene_regions[0][3] == 1 else head_]
            anchor2 = [gene_regions[0][0][0], gene_regions[0][0][2], head_ if gene_regions[0][3] == 1 else tail_]
            this_name = sorted([tuple(anchor1), tuple(anchor2)]) + [""]
            if tuple(this_name[:2]) not in name_counter:
                name_counter[tuple(this_name[:2])] = 1
            else:
                name_counter[tuple(this_name[:2])] += 1
                this_name[2] = "__copy" + str(name_counter[tuple(this_name[:2])])
            this_loc = [gene_regions[0][2], gene_regions[0][1], 1 * int(2 * ((anchor1 <= anchor2) - 0.5))]
            intergenic_regions.append([tuple(this_name)] + this_loc + [get_seq_with_gb_loc(this_loc), "noncoding"])
    elif len(gene_regions) > 1:
        first_region = gene_regions[0]
        circular_regions = [in_region for in_region in gene_regions if in_region[1] >= in_region[2]]
        if circular_regions:
            last_region = sorted(circular_regions, key=lambda x: (-x[2], x[1], x[0]))[0]
            end_of_last_region = last_region[2]
        else:
            last_region = gene_regions[-1]
        # if both of the terminal annotations across the ends (circular), they apparently overlapped
        if first_region[1] >= first_region[2] and last_region[1] >= last_region[2]:
            pass
        # elif embedded
        elif first_region[1] >= first_region[2]:
            pass
        elif last_region[1] >= last_region[2]:
            if last_region[2] >= first_region[1]:
                pass
            else:
                anchor1 = [last_region[0][0], last_region[0][2], tail_ if last_region[3] == 1 else head_]
                anchor2 = [first_region[0][0], first_region[0][2], head_ if first_region[3] == 1 else tail_]
                this_name = sorted([tuple(anchor1), tuple(anchor2)]) + [""]
                if tuple(this_name[:2]) not in name_counter:
                    name_counter[tuple(this_name[:2])] = 1
                else:
                    name_counter[tuple(this_name[:2])] += 1
                    this_name[2] = "__copy" + str(name_counter[tuple(this_name[:2])])
                this_loc = [last_region[2], first_region[1], 1 * int(2 * ((anchor1 <= anchor2) - 0.5))]
                intergenic_regions.append([tuple(this_name)] + this_loc + [get_seq_with_gb_loc(this_loc), "noncoding"])
        else:
            anchor1 = [last_region[0][0], last_region[0][2], tail_ if last_region[3] == 1 else head_]
            anchor2 = [first_region[0][0], first_region[0][2], head_ if first_region[3] == 1 else tail_]
            this_name = sorted([tuple(anchor1), tuple(anchor2)]) + [""]
            if tuple(this_name[:2]) not in name_counter:
                name_counter[tuple(this_name[:2])] = 1
            else:
                name_counter[tuple(this_name[:2])] += 1
                this_name[2] = "__copy" + str(name_counter[tuple(this_name[:2])])
            this_loc = [last_region[2], first_region[1], 1 * int(2 * ((anchor1 <= anchor2) - 0.5))]
            intergenic_regions.append([tuple(this_name)] + this_loc + [get_seq_with_gb_loc(this_loc), "noncoding"])
    go2 = 0
    while go2 < len(gene_regions) - 1:
        go_add = 1
        while go2 + go_add < len(gene_regions) and embed_in(gene_regions[go2 + go_add][1:3], gene_regions[go2][1:3]):
            go_add += 1
        if go2 + go_add == len(gene_regions):
            break
        this_region, next_region = gene_regions[go2], gene_regions[go2 + go_add]
        if this_region[1] >= this_region[2] and next_region[1] >= next_region[2]:
            pass
        elif this_region[2] < next_region[1] and end_of_last_region < next_region[1]:
            anchor1 = [this_region[0][0], this_region[0][2], tail_ if this_region[3] == 1 else head_]
            anchor2 = [next_region[0][0], next_region[0][2], head_ if next_region[3] == 1 else tail_]
            this_loc = [this_region[2], next_region[1], 1 * int(2 * ((anchor1 <= anchor2) - 0.5))]
            this_name = sorted([tuple(anchor1), tuple(anchor2)]) + [""]
            if tuple(this_name[:2]) not in name_counter:
                name_counter[tuple(this_name[:2])] = 1
            else:
                name_counter[tuple(this_name[:2])] += 1
                this_name[2] = "__copy" + str(name_counter[tuple(this_name[:2])])
            intergenic_regions.append([tuple(this_name)] + this_loc + [get_seq_with_gb_loc(this_loc), "noncoding"])
        go2 += go_add
    return gene_regions, intergenic_regions


def write_fasta(out_file, seq_dict, overwrite=False):
    names = sorted(list(seq_dict))
    if not overwrite and os.path.exists(out_file):
        existed_f = [x for x in os.listdir(os.path.split(out_file)[0])
                     if x.lower() == os.path.split(out_file)[-1].lower()]
        sys.stdout.write("Warning: running on a case-non-sensitive disk. "
                         "Cannot create " + out_file + " while " + existed_f[0] + " exists! "
                         "Appending seqs to " + existed_f[0] + " ... \n")
        with open(out_file, "a") as out_put_handler:
            for name in names:
                out_put_handler.write(">" + name + "\n" + seq_dict[name] + "\n\n")
    else:
        with open(out_file, "w") as out_put_handler:
            for name in names:
                out_put_handler.write(">" + name + "\n" + seq_dict[name] + "\n\n")


def write_statistics(out_file, base_name_list, gene_dict, intergene_dict):
    gene_names = sorted(list(gene_dict))
    str_gene_names = ["".join(n).replace(" ", "_") for n in gene_names]
    inter_names = sorted(list(intergene_dict))
    str_inter_names = ["--".join(["".join(x) for x in n[:2]]).replace(" ", "_") + n[2] for n in inter_names]
    with open(out_file, "w") as out_put_handler:
        out_put_handler.write("\t".join(["gb_name"] + str_gene_names + str_inter_names) + "\n")
        for gb_name in base_name_list:
            out_put_handler.write(gb_name)
            for loci_name in gene_names:
                if gb_name in gene_dict[loci_name]:
                    out_put_handler.write("\t" + str(len(gene_dict[loci_name][gb_name])))
                else:
                    out_put_handler.write("\t-")
            for loci_name in inter_names:
                if gb_name in intergene_dict[loci_name]:
                    out_put_handler.write("\t" + str(len(intergene_dict[loci_name][gb_name])))
                else:
                    out_put_handler.write("\t-")
            out_put_handler.write("\n")


def main():
    time0 = time.time()

    options, argv = get_options()
    options.gene_keys = options.gene_keys.split(",")
    gene_dir = os.path.join(options.out_put, "gene")
    intergenic_dir = os.path.join(options.out_put, "intergene")
    if not os.path.exists(options.out_put):
        os.mkdir(options.out_put)
        os.mkdir(gene_dir)
        os.mkdir(intergenic_dir)
    else:
        if options.overwrite:
            if not os.path.exists(gene_dir):
                os.mkdir(gene_dir)
            if not os.path.exists(intergenic_dir):
                os.mkdir(intergenic_dir)
        else:
            # raise FileExistsError(options.out_put + " exists!")
            raise IOError(options.out_put + " exists!")

    types = set()
    for this_t in options.gene_types.split(","):
        types.add(this_t)
        types.add(this_t.capitalize())
        types.add(this_t.lower())
        types.add(this_t.upper())

    out_gene_dict = {}
    out_intergenic_dict = {}
    base_name_list = []
    region_name_to_type = {}
    for this_gb in argv:
        if os.path.exists(this_gb):
            gb_base_name = os.path.basename(this_gb).replace(".gb", "").replace(".genbank", "")
            # base_name_list.append(gb_base_name)
            try:
                this_records = list(SeqIO.parse(this_gb, "genbank"))
            except ValueError as e:
                sys.stdout.write("Err loc: file " + this_gb + "\n")
                sys.stdout.write(str(e) + "\n")
                sys.exit()
            for go_record, seq_record in enumerate(this_records):
                try:
                    this_description = seq_record.description.replace("\n", " ").replace("\t", " ").strip()
                    this_seq_name = gb_base_name + \
                                    ("--" + str(go_record + 1)) * int(bool(len(this_records) > 1)) + \
                                    ("--" + this_description) * int(bool(this_description))
                    gene_regions, intergenic_regions = get_seqs(seq_record, types, options.gene_keys,
                                                                options.ignore_format_error,
                                                                options.product_to_gene,
                                                                options.case_treatment)
                    if options.one_copy:
                        temp_gene_dict = {}
                        for region_name, start, end, strand, this_seq, region_type in gene_regions:
                            # if region_name not in temp_gene_dict:
                            #     temp_gene_dict[region_name] = {}
                            temp_gene_dict[region_name] = this_seq
                            region_name_to_type[region_name] = region_type
                        temp_inter_dict = {}
                        for region_name, start, end, strand, this_seq, region_type in intergenic_regions:
                            # if region_name not in temp_inter_dict:
                            #     temp_inter_dict[region_name] = {}
                            region_name_to_type[region_name] = region_type
                            temp_inter_dict[region_name] = this_seq
                        # processing gene
                        go_to = 0
                        sorted_region_names = sorted(list(temp_gene_dict), key=lambda x: (x[0], x[2], x[1]))
                        while go_to < len(sorted_region_names):
                            region_name = sorted_region_names[go_to]
                            go_plus = 1
                            # if bool(next_loci[1]) == True, multiple copies exist.
                            while go_to + go_plus < len(sorted_region_names):
                                next_loci = sorted_region_names[go_to + go_plus]
                                if (next_loci[0], next_loci[2]) != (region_name[0], region_name[2]):
                                    if next_loci[1]:
                                        # if next_loci[0], next_loci[2])!=(region_name[0], region_name[2]
                                        # then next_loci seemed to be the first copy of a new region
                                        # but the first copy of a new region should not have __copy label (next_loci[1])
                                        # so next_loci is not a new region, but the same region with different exon
                                        sys.stdout.write("Warning: cannot find " + "".join(
                                            [next_loci[0], next_loci[2]]) +
                                             " while there's " + "".join(next_loci) + " in " + this_seq_name + "\n")
                                    break
                                else:
                                    go_plus += 1
                            if go_plus > 1:
                                this_seqs = []
                                for go_candidate in range(go_to, go_to + go_plus):
                                    this_seqs.append(temp_gene_dict[sorted_region_names[go_candidate]])
                                if len(set(this_seqs)) > 1:
                                    sys.stdout.write("Warning: distinct copies of " + "".join(
                                        region_name) + " in " + this_seq_name + "\n")
                                if options.copy_mode == "longest":
                                    temp_gene_dict[region_name] = sorted(this_seqs, key=lambda x: -len(x))[0]
                                elif options.copy_mode == "leastN":
                                    temp_gene_dict[region_name] = sorted(this_seqs, key=lambda x: count_n(x))[0]
                                elif options.copy_mode == "leastN_longest":
                                    temp_gene_dict[region_name] = sorted(
                                        this_seqs, key=lambda x: (count_n(x), -len(x)))[0]
                                for go_del in range(go_to + 1, go_to + go_plus):
                                    del temp_gene_dict[sorted_region_names[go_del]]
                            go_to += go_plus
                        # processing intergene
                        go_to = 0
                        sorted_inter_names = sorted(list(temp_inter_dict), key=lambda x: (x[:2], x[2]))
                        while go_to < len(sorted_inter_names):
                            inter_name = sorted_inter_names[go_to]
                            go_plus = 1
                            while go_to + go_plus < len(sorted_inter_names):
                                next_inter = sorted_inter_names[go_to + go_plus]
                                if inter_name[:2] != next_inter[:2]:
                                    if next_inter[2]:
                                        sys.stdout.write("Warning: cannot find " + "".join(next_inter[0]) + "--" +
                                            "".join(next_inter[1]) + " while there's " + "".join(next_inter[0]) +
                                            "--" + "".join(next_inter[1]) + next_inter[2] + " in " + this_seq_name + "\n")
                                    break
                                else:
                                    go_plus += 1
                            if go_plus > 1:
                                this_seqs = []
                                for go_candidate in range(go_to, go_to + go_plus):
                                    this_seqs.append(temp_inter_dict[sorted_inter_names[go_candidate]])
                                if len(set(this_seqs)) > 1:
                                    sys.stdout.write(
                                        "Warning: distinct copies of " + "".join(inter_name[0]) + "--" +
                                        "".join(inter_name[1]) + " in " + this_seq_name + "\n")
                                if options.copy_mode == "longest":
                                    temp_inter_dict[inter_name] = sorted(this_seqs, key=lambda x: -len(x))[0]
                                elif options.copy_mode == "leastN":
                                    temp_inter_dict[inter_name] = sorted(this_seqs, key=lambda x: count_n(x))[0]
                                elif options.copy_mode == "leastN_longest":
                                    temp_inter_dict[inter_name] = sorted(
                                        this_seqs, key=lambda x: (count_n(x), -len(x)))[0]
                                for go_del in range(go_to + 1, go_to + go_plus):
                                    del temp_inter_dict[sorted_inter_names[go_del]]
                            go_to += go_plus
                        # transfer temp to main dict
                        for region_name in temp_gene_dict:
                            if region_name not in out_gene_dict:
                                out_gene_dict[region_name] = {}
                            out_gene_dict[region_name][this_seq_name] = temp_gene_dict[region_name]
                        for region_name in temp_inter_dict:
                            if region_name not in out_intergenic_dict:
                                out_intergenic_dict[region_name] = {}
                            out_intergenic_dict[region_name][this_seq_name] = temp_inter_dict[region_name]
                    else:
                        for region_name, start, end, strand, this_seq in gene_regions:
                            if region_name not in out_gene_dict:
                                out_gene_dict[region_name] = {}
                            out_gene_dict[region_name][this_seq_name] = this_seq
                        for region_name, start, end, strand, this_seq in intergenic_regions:
                            if region_name not in out_intergenic_dict:
                                out_intergenic_dict[region_name] = {}
                            out_intergenic_dict[region_name][this_seq_name] = this_seq
                    base_name_list.append(this_seq_name)
                except NotImplementedError as e:
                    sys.stdout.write("Err loc: " + str(go_record + 1) + "th record in file " + this_gb + "\n")
                    sys.stdout.write("\nSolutions: "
                                     "\n1. Add suitable key (in above qualifiers part) to \"--keys\".")
                    sys.stdout.write("\n2. Add suitable present key and its value to the problematic annotation record.")
                    sys.stdout.write("\n3. Use \"--ignore-format-error\" to skip this annotation record.\n")
                    sys.exit()
        else:
            sys.stdout.write("")

    if options.combine_exon:
        regions_with_exon = [x for x in list(out_gene_dict) if x[2]]
        region_set_dict = {}
        region_set_types = {}
        for region_name in regions_with_exon:
            region_set_name = region_name[:2]
            exon_num = int(region_name[2].replace("__exon", ""))
            if region_set_name not in region_set_dict:
                region_set_dict[region_set_name] = []
                region_set_types[region_set_name] = set()
            region_set_dict[region_set_name].append(exon_num)
            region_set_types[region_set_name].add(region_name_to_type[region_name])
        for region_set_name in region_set_dict:
            region_set_dict[region_set_name].sort()
            seq_names = set()
            for exon_num in region_set_dict[region_set_name]:
                for gb_name in out_gene_dict[tuple(list(region_set_name) + ["__exon" + str(exon_num)])]:
                    seq_names.add(gb_name)
            new_name = tuple(list(region_set_name) + [""])
            if new_name not in out_gene_dict:
                out_gene_dict[new_name] = {}
                region_name_to_type[new_name] = "_".join(sorted(region_set_types[region_set_name]))
            for gb_name in seq_names:
                out_gene_dict[new_name][gb_name] = ""
                for exon_num in region_set_dict[region_set_name]:
                    out_gene_dict[new_name][gb_name] += \
                        out_gene_dict[tuple(list(region_set_name) + ["__exon" + str(exon_num)])].get(gb_name, "")
            for exon_num in region_set_dict[region_set_name]:
                del out_gene_dict[tuple(list(region_set_name) + ["__exon" + str(exon_num)])]

    if options.mix:
        gene_mixture = dict()
        for region_name in out_gene_dict:
            region_name_str = "".join(region_name).replace(" ", "_")
            for this_seq_name_ in out_gene_dict[region_name]:
                new_seq_name_ = region_name_str + " " + region_name_to_type.get(region_name, "region") + \
                                " - " + this_seq_name_.replace(" ", "_").replace(",", "_")
                gene_mixture[new_seq_name_] = out_gene_dict[region_name][this_seq_name_]
        write_fasta(os.path.join(gene_dir, "gene.fasta"), gene_mixture, overwrite=options.overwrite)
        intergene_mixture = dict()
        for region_name in out_intergenic_dict:
            region_name_str = "--".join(["".join(x) for x in region_name[:2]]).replace(" ", "_") + region_name[2]
            for this_seq_name_ in out_intergenic_dict[region_name]:
                new_seq_name_ = region_name_str + " " + region_name_to_type.get(region_name, "region") + \
                                " - " + this_seq_name_.replace(" ", "_").replace(",", "_")
                intergene_mixture[new_seq_name_] = out_intergenic_dict[region_name][this_seq_name_]
        write_fasta(os.path.join(intergenic_dir, "intergene.fasta"), intergene_mixture, overwrite=options.overwrite)
    else:
        for region_name in out_gene_dict:
            write_fasta(os.path.join(gene_dir, "".join(region_name).replace(" ", "_") + ".fasta"),
                        out_gene_dict[region_name],
                        overwrite=options.overwrite)
        for region_name in out_intergenic_dict:
            write_fasta(os.path.join(intergenic_dir,
                                     "--".join(["".join(x) for x in region_name[:2]]).replace(" ", "_") +
                                     region_name[2] + ".fasta"), out_intergenic_dict[region_name],
                        overwrite=options.overwrite)
    write_statistics(os.path.join(options.out_put, "statistics.txt"), base_name_list, out_gene_dict, out_intergenic_dict)

    sys.stdout.write("Time cost: "+str(time.time() - time0) + "\n")


if __name__ == '__main__':
    main()
