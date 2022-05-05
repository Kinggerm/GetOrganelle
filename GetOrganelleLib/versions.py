

def get_versions():
    return versions[0]["number"]


versions = [
    {
        "number": "1.7.6.0",
        "features": [
            "1. improve the target component recognition on non-circular cases (discussion 138 & issue 141)",
            "2. assembly_graph.py & statistical_func.py: specify scipy error (issue 132)",
            "3. compatible with the newly released GetOrganelleDB v0.0.1.minima",
            "4. get_organelle_config.py: fix a bug when there was not directory made, --config-dir is invalid",
        ],
        "time": "2022-05-05 16:20 UTC-5"
    },
    {
        "number": "1.7.5.3",
        "features": [
            "1. get_organelle_from_reads.py: fix a bug when no qualified reads found (issue 123)",
            "2. get_organelle_from_reads.py: pass --spades-options to pre-assembly for Mac M1 (issue 127)",
            "3. get_organelle_from_reads.py:slim_spades_result: avoid ValueError",
            "4. Update README.md: input read qc; issues->discussions",
        ],
        "time": "2022-01-21 11:20 UTC-5"
    },
    {
        "number": "1.7.5.2",
        "features": [
            "1. ask the questions publicly",
            "2. statistical_func.py: weighted_gmm_with_em_aic(): fix a bug that will be triggered by "
            "   graph produced by join_spades_fastg_by_blast.py (reported by Mergi Dinka); "
            "   also fix a hidden mis-indexing issue there",
            "3. Utilities/join_spades_fastg_by_blast.py: update on a v name issue "
            "   (issues 119)",
        ],
        "time": "2021-12-15 02:35 UTC-5"
    },
    {
        "number": "1.7.5.1",
        "features": [
            "1. make_batch_for_get_organelle.py: usage updated (reported by Fei Zhao@KIB)",
            "2. space in the output path should be forbidden (reported by Manuela Sann@UFreiburg)",
        ],
        "time": "2021-05-13 16:50 UTC+8"
    },
    {
        "number": "1.7.5",
        "features": [
            "1. assembly_parser.py: fix a bug in estimation of the multiplicity of self-loop vertex, "
            "   which was falsely forced to be at least 2. (detected in a case of Yan Zhong@SCNU)",
            "2. get_organelle_from_reads.py: typo in instruction for --max-extending-len corrected",
            "3. redundant Assembly.merging_history() removed",
            "4. fix a small bug in Assembly.is_sequential_repeat(), "
            "   isolate path_without_leakage as find_pair_closing_the_path",
            "5. get_organelle_from_reads.py: fix an inconformity between the document and the default for -R, ",
            "6. pipe_control_func.py: detect bowtie2 version improved",
        ],
        "time": "2021-05-13 16:50 UTC+8"
    },
    {
        "number": "1.7.4.1",
        "features": [
            "1. get_organelle_config.py: provide guidance for old code and new database incompatibility "
            "   (reported by Wenxiang Liu@SWFU)",
            "2. assembly_parser.py: fix a bug after scaffolding with SPAdes path "
            "   (introduced in 1.7.4 feature 5; reported by Robin van Velzen@WUR)",
            "3. update README.md with improved instruction",
        ],
        "time": "2021-04-16 14:46 UTC+8"
    },
    {
        "number": "1.7.4",
        "features": [
            "1. SequenceList: self.__indexed updated",
            "2. get_static_html_context: requests.exceptions.ReadTimeout added"
            "3. get_organelle_from_reads.py: '--overwrite' added; catch shutil.rmtree() errors; no python-lib -> error."
            "4. get_organelle_from_*.py: turning off bandage if the result is not circular",
            "5. assembly_parser.py: recording every overlap value rather than using a universal value for hifiasm "
            "   (in response to Christopher Benson@PSU)",
            "6. optparse -> argparse (in response to Matthias Bernt@UFZ)",
            "7. get_organelle_from_assembly.py: fix a bug with the malfunction of --continue when the input is gfa ",
            "8. get_organelle_from_reads.py/disentangle_organelle_assembly.py: correct typos",
            "9. pipe_control_func.py: map_with_bowtie2: warn reads integrity; build_bowtie2_db: rm small index",
            "10. get_organelle_config.py: verbose log for bowtie2 and blast",
            "11. update README.md with a reframed instruction",
        ],
        "time": "2021-04-14 17:52 UTC+8"
    },
    {
        "number": "1.7.4-pre2",
        "features": [
            "1. README.md: updated",
            "2. partial fix: subprocess may 'fail' if error was in the directory name due to detecting error from log",
            "3. scaffolding failures: try except AssertError to skip abnormal SPAdes paths (reported by Jinjing Jian@FDU)",
        ],
        "time": "2021-03-04 23:55 UTC+8"
    },
    {
        "number": "1.7.4-pre",
        "features": [
            "1. setup.py: modify scripts with utf-8",
            "2. filtered -> extended, to clarify the process",
            "3. get_organelle_from_reads.py: fix illustration of --which-blast, --which-bowtie2, --which-spades",
            "4. fix a bug of pre_assembly_mapped_reads_for_base_cov introduced in 1.7.3.5b",
        ],
        "time": "2021-02-26 11:30 UTC+8"
    },
    {
        "number": "1.7.3.5b",
        "features": [
            "1. pipe_control_func: executable modified",
            "2. realtime monitoring SPAdes log (solving the stuck in the segmentation fault)",
        ],
        "time": "2021-02-24 18:11 UTC+8"
    },
    {
        "number": "1.7.3.5a",
        "features": [
            "1. setup.py: fix a bug of codec for in-situ installation: invalid attempt to modify ._*py files",
            "2. setup.py: fix a bug introduced by 1.7.3 while relocating GetOrganelle databases",
        ],
        "time": "2021-02-24 12:30 UTC+8"
    },
    {
        "number": "1.7.3.5",
        "features": [
            "1. README.md updated with embplant_mt notes",
            "2. remove invalid --genes file check",
        ],
        "time": "2021-02-23 01:11 UTC+8"
    },
    {
        "number": "1.7.3.4",
        "features": [
            "1. fix bugs with '--config-dir'",
            "2. customized databases (--genes/--ex-genes) passed to slim during pre-assembly & depth estimation",
            "3. fix a bug with get_graph_coverages_range_simple() when no contigs received"
        ],
        "time": "2021-02-18 19:00 UTC+8"
    },
    {
        "number": "1.7.3.4-pre",
        "features": [
            "1. fix bugs with '--config-dir'",
            "2. customized databases (--genes/--ex-genes) passed to slim during pre-assembly & depth estimation",
        ],
        "time": "2021-02-12 02:00 UTC+8"
    },
    {
        "number": "1.7.3.3",
        "features": [
            "1. early termination on invalid path characters for spades",
            "2. fix a bug introduced by '--max-reads inf'",
            "3. get_organelle_config.py: fix a bug if a new organelle types was added and '--use-local' was used",
        ],
        "time": "2021-02-11 01:00 UTC+8"
    },
    {
        "number": "1.7.3.2",
        "features": [
            "1. improve support for standard gfa format: E/L with optional fields",
        ],
        "time": "2021-02-03 13:30 UTC+8"
    },
    {
        "number": "1.7.3.1",
        "features": [
            "1. bug fixed: a bug of sorting consensus vertex names using smart_trans_for_sort",
            "2. option --max-reads inf configured",
        ],
        "time": "2021-01-25 13:00 UTC+8"
    },
    {
        "number": "1.7.3",
        "features": [
            "1. fungus_nr added",
            "2. relocate default GetOrganelle databases using GETORG_PATH",
            "3. log platform info",
        ],
        "time": "2021-01-20 12:30 UTC+8"
    },
    {
        "number": "1.7.2b",
        "features": [
            "1. --reverse-lsc malfunction issue solved",
        ],
        "time": "2020-12-19 14:00 UTC+8"
    },
    {
     "number": "1.7.2a",
     "features": [
         "1. slim_graph.py: --evalue added",
         "2. trans_word_cov: using log() to avoid large int converting",
         "3. get_organelle_from_reads.py: add url to FAQ for incomplete result",
     ],
     "time": "2020-12-16 23:50 UTC+8"
    },
    {
     "number": "1.7.2",
     "features": [
         "1. assembly_parser.py: fix bugs in scaffolding",
     ],
     "time": "2020-12-10 23:50 UTC+8"
    },
    {
     "number": "1.7.2beta2",
     "features": [
         "1. get_organelle_from_reads.py: fix a bug in parsing options",
     ],
     "time": "2020-12-09 UTC+8"
    },
    {
     "number": "1.7.2beta",
     "features": [
         "1. slim_graph.py: fix a bug in parsing options (reported by Fei Zhao @ KIB)",
         "2. automatically adding --phred-offset to avoid hammer failures",
         "3. MergingHistory and ConsensusHistory added, in replace of using the names to detect merging history",
         "4. Compatible with flye gfa format",
         "5. Scaffolding function bug fixed & improved",
         "6. get_organelle_config.py: example typo fixed",
     ],
     "time": "2020-12-03 UTC+8"
    },
    {
     "number": "1.7.1a",
     "features": [
         "output index bug fixed",
     ],
     "time": "2020-07-28 01:30 UTC+8"
    },
    {
     "number": "1.7.1",
     "features": [
         "get_organelle_from_assembly.py: do not convert gfa to fastg; ",
         "Assembly.get_all_circular_paths(): optimized for plastome order",
         "Assembly.reduce_to_subgraph: limit_offset_current_vertex -> bait_offsets; safer slim_graph.py performance",
         "Assembly.get_all*_paths(): more detailed log info",
     ],
     "time": "2020-07-25 02:50 UTC+8"
    },
    {
     "number": "1.7.0c",
     "features": [
         "Utilities/slim_graph.py: fix a bug with anonym mode",
         "README.md: updated"
     ],
     "time": "2020-07-21 10:55 UTC+8"
    },
    {
     "number": "1.7.0b",
     "features": [
         "1. get_organelle_from_reads.py: "
         "  1) --ignore-k work for small k disentanglement "
         "  2) fix a bug when input reads are very few ",
         "2. better log info",
     ],
     "time": "2020-07-08 22:30 UTC+8"
    },
    {
     "number": "1.7.0",
     "features": [
         "1. get_organelle_from*.py: reorganize some importing code, fix minor issues",
         "2. get_organelle_from_reads.py: rm output/seed/embplant_pt.initial.fq.spades by default",
     ],
     "time": "2020-06-28 00:00 UTC+8"
    },
    {
     "number": "1.7.0-beta7",
     "features": [
         "1. get_organelle_from_reads.py: "
         "   1) fix option --max-extending-len typo "
         "   2) --disentangle-time-limit 600 => 1800"
         "   3) fix parameter estimation, pre-assembly minor bug",
         "2. get_organelle_from_assembly.py: echo scaffolding.",
         "3. assembly_parser.Assembly: "
         "   1) --keep-temp generate more intermediate results "
         "   2) Vertex: connections[end]: set() -> OrderedDict() "
         "   3) fix a bug of multiplicity estimation on assembly graph with self-loop contigs",
         "4. assembly_parser.SPAdesScaffolds: improved with more situations",
         "5. README.md: updated",
     ],
     "time": "2020-06-26 03:00 GMT-6"},
    {
     "number": "1.7.0-beta6",
     "features": [
         "1. get_organelle_from_reads.py: make pre-assembly and --ignore-k work for small read length",
     ],
     "time": "2020-06-07 01:40 GMT-6"},
    {
     "number": "1.7.0-beta5",
     "features": [
         "1. get_organelle_config.py: alternative repository (gitee.com/jinjianjun/GetOrganelleDB) added",
         "2. setup.py: dependent python lib requests added",
     ],
     "time": "2020-05-28 18:40 GMT-6"},
    {"number": "1.7.0-beta4",
     "features": [
         "1. Utilities/slim_fastg -> Utilities/slim_graph: 1) reorganized; 2) added support for gfa format graph file;"
             "3) --max-slim-extending-len added",
         "2. get_organelle_config.py added with *Database removed",
         "3. get_organelle_from_reads.py: "
             "1) use SPAdes generated scaffolds.paths to create gap containing scaffolds (assembly_parser.py)"
             "2) rm --gradient-k "
             "3) output fasta name modified"
             "4) log Database version",
         "4. get_organelle_from_assembly.py: "
             "1) use SPAdes generated scaffolds.paths to create gap containing scaffolds (assembly_parser.py)"
             "2) output fasta name modified"
             "3) log Database version",
         "5. assembly_parser.py: 1) merge_all_possible_vertices: fix a bug for overlap=0; 2) reduce_to_subgraph added; "
             "3) processing_polymorphism: fix a bug; "
             "4) class SPAdesScaffolds: use SPAdes generated scaffolds.paths to create gap containing scaffolds",
         "6. Utilities/reconstruct_graph_from_fasta.py: fix a bug",
         "7. Bandage generate temp file",
         "8. README.md: updated",
     ],
     "time": "2020-05-26 21:00 GMT-6"},
    {"number": "1.6.4",
     "features": [
         "1. log plastome structure",
         "2. --max-paths-num added for get_organelle*.py & disentangle*.py",
         "3. reorganize codes: class SimpleAssembly & detect_plastome_architecture()",
         "4. evaluate_assembly_using_mapping.py: --stat-mode, --bowtie2-options, --plot-font added",
         "5. isolate GetOrganelleDep again",
         "6. README.md: updated with conda installation",
     ],
     "time": "2020-02-27 17:14 GMT-6"},
    {"number": "1.6.3a",
     "features": [
         "1. Minor bugs fixes",
     ],
     "time": "2020-02-27 17:14 GMT-6"},
    {"number": "1.6.3-beta",
     "features": [
         "1. get_organelle_from_assembly.py & disentangle_organelle_assembly.py: --max-multiplicity added",
         "2. Assembly.estimate_copy_and_depth_precisely() modified: constraint_max_function() for --max-multiplicity",
         "3. Assembly.tag_in_between() modified",
         "4. Assembly.estimate_copy_and_depth_by_cov() modified: min average coverage limit",
         "5. Assembly.processing_polymorphism():"
         "   fix a bug when kmer-len repeats shared by two contigs; fix a bug that cause RuntimeError",
         "6. Assembly: too many results due to palindromic repeats, problem solved",
         "7. Utilities/reconstruct_graph_from_fasta.py & NaiveKmerNodeGraph added",
         "8. Utilities/gfa_to_fasta.py, Utilities/fastg_to_gfa.py: description corrected",
         "9. Assembly.parse_gfa(): compatibility increased",
         "10. Utilities/gfa2fastg.py: compatibility increased",
         "11. Assembly.estimate_copy_and_depth_precisely(): fix a bug on a rare case that multiplicities res are 4,8,4",
         "12. README.md: updated",
     ],
     "time": "2020-02-22 02:40 GMT-6"},
    {"number": "1.6.2e",
     "features": [
         "1. seq_parser.py: fix a bug for fastq format: @*****#/1",
         "2. get_organelle_from_reads.py: separate_fq_by_pair(), fix a bug when detecting pair info failed; ",
         "3. evaluate_assembly_using_mapping.py: fix a bug for --plot-transparent",
         "4. GetOrganelleLib.__init__.py: __version__",
         "5. README.md: updated",
     ]},
    {"number": "1.6.2d",
     "features": [
         "1. get_organelle_from_reads.py: fix a bug with '-F anonym'",
     ]},
    {"number": "1.6.2c",
     "features": [
         "1. GetOrganelleLib/assembly_parser.py: SSC direction set according to orf",
         "2. disentangle: --reverse-lsc option added; fix a bug for disentangling contigs with no overlaps",
         "3. Utilities/plastome_arch_info.py: GC content added",
         "4. get_organelle_from_reads.py: fix a bug for --flush-step inf"
     ]},
    {"number": "1.6.2b",
     "features": [
         "1. fix a minor bug when raising ProcessingGraphFailed with # tags",
         "2. setup.py install modified",
         "3. open() modified",
     ]},
    {"number": "1.6.2a",
     "features": [
         "1. the bug with option \"--genes\" fixed",
         "2. the bug with \"Mixing iteration and read methods\" introduced by 1.6.2 fixed",
     ]},
    {"number": "1.6.2",
     "features": [
         "1. get_organelle_from_reads.py: --reduce-reads-for-cov/estimate_maximum_n_reads_using_mapping() added; "
         "   problem with pre_assembly_mapped_reads_for_base_cov() fixed; "
         "   better target-hitting base coverage estimation",
         "2. get_organelle_from_assembly.py: fix a bug on parsing gfa file with long seq head names; "
         "   --keep-temp fixed; fix a bug with '-h'; ",
         "3. Utilities/slim_fastg.py: --no-merge -> --merge; disable merge by default",
         "4. GetOrganelleLib/assembly_parser.py: fix a bug with generating new vertices, "
         "   as well as merge_all_possible_contigs; export plastome-LSC direction according to convention based on "
         "   accumulated orf lengths (the conventional reverse direction has more accumulated orf lengths), which "
         "   makes users easier to use; remove processing_polymorphism() before filter_by_coverage() to better "
         "   cluster organelle contigs by coverages",
         "5. blastn GLIBC_2.14 not found problem fixed",
         "6. '-F embplant_pt' does not remove embplant_mt-hitting contigs, which makes more accurate clustering",
     ]},
    {"number": "1.6.1a",
     "features": [
         "1. GetOrganelleLib/SeedDatabase: embplant_pt updated. ",
     ]},
    {"number": "1.6.1",
     "features": [
         "1. GetOrganelleLib/SeedDatabase: updated with repeats removed. "
         " save computational time and generate better target-hitting base coverage estimation.",
     ]},
    {"number": "1.6.0",
     "features": [
         "1. setup.py added with new installation way.",
         "2. GetOrganelleDep added for easier dependencies installation",
         "3. get_organelle_reads.py -> get_organelle_from_reads.py;"
         " --max-extending-len, --ex-genes, --which-blast, --which-bowtie2, --which-spades, --zip-files added; "
         "multi-organelle mode supported; add support for fq head @digits; "
         "--safe,-r,--auto-wss,--soft-max-n-words etc removed; --flush-step modified (background mode); "
         "4. get_organelle_from_assembly.py (basically slim + disentangle) added",
         "5. Library/SeqReference -> GetOrganelleLib/SeedDatabase",
         "6. Library/NotationReference -> GetOrganelleLib/LabelDatabase",
         "7. plant_mt -> embplant_mt; plant_nr -> embplant_nr; plant_pt -> embplant_pt; other_pt added;",
         "8. assembly_parser.py: keep terminal contigs if --linear; fix a bug when --acyclic-allowed; "
         "optimized for self-loop contigs",
         "9. Utilities/evaluate_assembly_using_mapping.py: log statistics without --draw",
         "10. Utilities/disentangle_organelle_assembly.py: --acyclic-allowed -> --linear",
         "11. Utilities/slim_fastg.py: --no-merge added",
     ]},
    {"number": "1.5.2a",
     "features": [
         "1. more descriptive log",
     ]},
    {"number": "1.5.2",
     "features": [
         "1. get_organelle_reads.py: Bowtie2 index files written to sample output directory rather than "
         " to GetOrganelle/Library/SeqReference",
         "2. get_organelle_reads.py: more descriptive log",
         "3. seq_parser.py: re_linear_circular_seqs: fix a bug for evaluating the assembly result of DR plastomes",
         "4. evaluate_assembly_using_mapping.py: bowtie2-build --seed",
     ]},
    {"number": "1.5.1c",
     "features": [
         "1. --random-seed added with default value 12345",
         "2. wider suggested/default k-mer values",
         "3. get_organelle_reads.py: exit after illegal -F",
         "4. evaluate_assembly_using_mapping.py: customized error rate info added",
         "5. evaluate_assembly_using_mapping.py: robust to illegitimate usage of duplicated seq names in fasta",
         "6. evaluate_assembly_using_mapping.py: fix a bug when no aligned bases found",
         "7. sam_parser.py: keep redundant cigar chars",
         "8. README.md: toolkit",
         "9. evaluate_assembly_using_mapping.py: re-linearize circular sequence before mapping",
         "10. plastome_arch_info.py: added",
     ]},
    {"number": "1.5.1b",
     "features": [
         "1. get_organelle_reads.py: value of mesh size should have effect on --out-per-round (fix a bug since 1.4.2)",
     ]},
    {"number": "1.5.1a",
     "features": [
         "1. get_organelle_reads.py: from math import inf is not compatible with Python2; -R default set to 1000",
         "2. pipe_control_func.py: MEM_TRANS, influence summary_get_organelle_output.py",
     ]},
    {"number": "1.5.1",
     "features": [
         "1. fix a bug in get_organelle_reads.py: pre_grouping(): generate_forward_and_reverse()",
         "2. fix a bug in evaluate_assembly_using_mapping.py: --debug",
     ]},
    {"number": "1.5.0h",
     "features": [
         "1. re-organize importing codes",
         "2. minimum of -R: 2 -> 1",
         "3. slim_fastg.py: remove default -F",
         "4. round_statistics.py: increase significant digits",
     ]},
    {"number": "1.5.0g",
     "features": [
         "1.get_organelle_reads.py: fix a bug in --out-per-round & --min-quality-score, chop_seqs -> chop_seq_list",
         "2.get_organelle_reads.py: expand user-defined word size scope, 49 -> 29 (AUTO_MIN_WS, GLOBAL_MIN_WS)",
         "3.README.md: updated",
     ]},
    {"number": "1.5.0f",
     "features": [
         "1.disentangle: more instructive log.",
         "2.Set default logging level of round_statistics.py and evaluate_assembly_using_mapping.py to INFO",
         "3.round_statistics.py: set larger value to max_cov_tick",
     ]},
    {"number": "1.5.0e",
     "features": [
         "1.get_organelle_reads.py: --continue skip disentangling when *.gfa & *.fasta",
     ]},
    {"number": "1.5.0d",
     "features": [
         "1.disentangle: parallel contigs remained; --contamination-depth",
         "2.seq_parser.py: find_string_dif adjusted",
     ]},
    {"number": "1.5.0c",
     "features": [
         "1.evaluate_assembly_using_mapping.py & sam_parser.py: mapped reads counted; echo bug correctly",
     ]},
    {"number": "1.5.0b",
     "features": [
         "1.evaluate_assembly_using_mapping.py: --debug added",
         "2.slim_fastg.py: compatible with older python version"
     ]},
    {"number": "1.5.0a",
     "features": [
         "1.evaluate_assembly_using_mapping.py: modifying the layout; plot options added",
     ]},
    {"number": "1.5.0",
     "features": [
         "1.evaluate_assembly_using_mapping.py added",
     ]},
    {"number": "1.5.0-pre2",
     "features": [
         "-F anonym added",
         "-F fungus_mt added",
         "-F animal_mt added but not activated",
         "re-estimate base coverage by counting seed word frequencies if the result (directly from sam) < 200",
         "fix a bug for logging base-coverage when no kmer detected from graph",
         "fix a bug of --continue",
     ]},
    {"number": "1.5.0-pre",
     "features": [
         "cp -> plant_pt; mt -> plant_mt; nr -> plant_nr; for adding animals",
         "Comparison url (https://github.com/Kinggerm/GetOrganelleComparison) added",
     ]},
    {"number": "1.4.4b",
     "features": [
         "1.assembly_parser.py: fix a bug for disentangling single-contig graph; remove redundant 'repeat_pattern';",
     ]},
    {"number": "1.4.4a",
     "features": [
         "time limit works only for disentangling graph as a circular genome",
         "more informative log info for disentangling",
     ]},
    {"number": "1.4.4",
     "features": [
         "1.get_organelle_reads.py: fix a bug with --continue & --prefix when LogInfo() added; ",
         "2.assembly_parser.py & statistical_func.py: "
         "if single copy vertex percentage is < 50%, continue dropping suspicious vertices",
         "3.pip_control_func.py: for --prefix",
     ]},
    {"number": "1.4.3a",
     "features": [
         "1.get_organelle_reads.py: check_kmers() modified; ",
         "2.pipe_control_func.py: LogInfo() modified",
     ]},
    {"number": "1.4.3",
     "features": [
         "1.get_organelle_reads.py: output renamed; fix a bug of logging",
         "2.summary_get_organelle_output.py: added",
     ]},
    {"number": "1.4.3-beta",
     "features": [
         "1.get_organelle_reads.py: a bug in logging seed reads; moving re-setting kmers after extending;",
         "2.disentangle_organelle_assembly.py & assembly_parser.py: "
         "2a.'--acyclic-allowed' activated; "
         "2b.'--continue' added; "
         "2c.better output for polymorphyism-contained graph, default degenerate similarity threshold increased, "
         "   print warning when degenerate base used; "
         "2d.find_target_graph(broken_graph_allowed)",
         "4.pip_control_func.py: LogInfo added",
         "5.NotationReference: cp updated"
     ]},
    {"number": "1.4.2",
     "features": [
         "1.get_organelle_reads.py: better seed fastq file log; "
         " increase the default values of jump_step and mesh_size",
         "2.assembly_parser.py: fix a bug in filter_by_coverage().",
         "3.disentangle_organelle_assembly.py: '--acyclic-allowed' added (not activated yet)",
         "4.statistical_func.py: fix a bug in assign_cluster_labels",
         "5.join_spades_fastg_by_blast.py: signs of gap and overlap",
         "6.SeqReference/cp.fasta: source id specified",
     ]},
    {"number": "1.4.1a",
     "features": [
         "1.get_organelle_reads.py: '--no-pre-reading' added",
     ]},
    {"number": "1.4.1",
     "features": [
         "1.assembly_parser.py: Assembly.export_path() and Assembly.merge_all_possible_vertices():"
         " name of merged vertices optimized",
         "2.README: PATH configuration",
         "3.mk_batch_for_iteratively_mapping_assembling.py: -t threads",
         "4.get_organelle_reads.py --fast mode modified",
         "5.sam_parser.py added: for 1.5.0",
     ]},
    {"number": "1.4.0j",
     "features": [
         "1.default values (--max-n-words & --auto-wss) set to make GetOrganelle perform like older versions",
     ]},
    {"number": "1.4.0i",
     "features": [
         "1.empirically reduce maximum word size",
         "2.report SPAdes failed when not output folder exist.",
     ]},
    {"number": "1.4.0h",
     "features": [
         "1.fix a bug: calling slim_fastg.py failed.",
     ]},
    {"number": "1.4.0g",
     "features": [
         "1.slim_fastg.py: fix the import error when using python2.*",
         "2.README.md: in case of HTTP request failed",
     ]},
    {"number": "1.4.0f",
     "features": [
         "1.parse_gfa() added to Library/assembly_parser.py",
         "2.get_organelle_reads.py -h",
     ]},
    {"number": "1.4.0e",
     "features": [
         "1.print python version",
         "2.gfa2fastg.py modified, gfa2fasta.py added, fastg2gfa.py added",
     ]},
    {"number": "1.4.0d",
     "features": [
         "1.some default values adjusted.",
         "2.slim_fastg.py: '--depth-threshold' -> '--min-depth'&'--max-depth'",
         "3.print python version"
         "4.gfa2fastg.py modified, gfa2fasta.py added, fastg2gfa.py added"
     ]},
    {"number": "1.4.0c",
     "features": [
         "1.'--pre-w' added mainly for reproducing results when word size changes during reads extending process.",
     ]},
    {"number": "1.4.0b",
     "features": [
         "1.--max-reads also works for mapping now, which would do better in target coverage estimation.",
     ]},
    {"number": "1.4.0a",
     "features": [
         "1.default reference seq set as Library/SeqReference/*.fasta.",
     ]},
    {"number": "1.4.0",
     "features": [
         "1.estimate_word_size() added.",
         "2.auto_word_size_step (--auto-wss, --soft-max-words, -r) added.",
         "3.mean_error_rate added.",
         "4.options re-organized and introductions optimized: '-h' versus '--help'.",
         "5.Utilities/mk_get_organelle.py recognize *.gz files.",
         "6.change the default setting of --max-n-words"
     ]},
    {"number": "1.3.1",
     "features": [
         "1.'--max-discard-percent' added to prevent discarding too much data.",
         "2.fix the bug in get_low_quality_char_pattern, which causes misidentification for quality encoding format.",
         "3.--flush-frequency added",
         "4.better log info",
         "5.'--overwrite' -> '--no-overwrite' in slim_fastg.py"
     ]},
    {"number": "1.3.0d",
     "features": [
         "1.continue to process assembly results based on available kmers if SPAdes failed at one of the planned kmers",
     ]},
    {"number": "1.3.0c",
     "features": [
         "1.fix a bug: compatibility with '--continue' option of SPAdes",
     ]},
    {"number": "1.3.0b",
     "features": [
         "1.fix a bug for exporting organelle",
     ]},
    {"number": "1.3.0a",
     "features": [
         "1.automatically discard improper input kmers",
     ]},
    {"number": "1.3.0",
     "features": [
         "1.Read quality control (--min-quality-score) added.",
         "2.--trim option removed.",
         "3.fix a bug on word size estimation",
     ]},
    {"number": "1.2.0d",
     "features": [
         "1.Add --max-words.",
     ]},
    {"number": "1.2.0c",
     "features": [
         "1.Go over assemblies based on all kmer values, from large to small, until the solvable assembly is found.",
         "2.overwrite option added for slim_fastg.py",
         "3.Optimize the log",
         "4.multiprocessing function added (planed, not utilized yet)",
     ]},
    {"number": "1.2.0b",
     "features": [
         "1.Assembly.parse_fastg(): (more robust) Add connection information to both of the related vertices"
         " even it is only mentioned once;",
         "2.Assembly.is_sequential_repeat(): fix a bug that leak in the reverse direction;",
         "3.add depth_factor to the main script;",
         "4.remove unnecessary warning when #read equals maximum#reads setting, show warning only when overrunning;",
     ]},
    {"number": "1.2.0a",
     "features": [
         "1.set time limit for disentangling",
     ]},
    {"number": "1.2.0",
     "features": [
         "1.more robust and precise in disentangling graph: ",
         "2.report contamination; ",
         "3.detect parallel contigs and generate consensus;",
         "4.estimate chloroplast coverage distribution pattern using weighted GMM with EM and BIC;",
         "5.re-organize codes",
         "6.update NotationReference",
     ]},
    {"number": "1.1.0d",
     "features": [
         "1.more precise in disentangling graph.",
         "2.--prefix added",
     ]},
    {"number": "1.1.0c",
     "features": [
         "1.--max-reads: default=10,000,000 for cp and nr, default=50,000,000 for mt",
     ]},
    {"number": "1.1.0b",
     "features": [
         "1.'-w' could be followed with a ratio of word_size to average_read_length now",
     ]},
    {"number": "1.1.0a",
     "features": [
         "1.Add options-maximum_number_of_reads with default=10,000,000 to avoid unnecessary overloading",
     ]},
    {"number": "1.1.0",
     "features": [
         "1.automatically exporting final sequence(s)",
         "2.adding disentangle_organelle_assembly.py",
         "3.adding assembly_parser.py",
         "4.re-organize codes",
         "5.revise README.md",
     ]},
    {"number": "1.0.5",
     "features": [
         "1.re-organize codes",
     ]},
    {"number": "1.0.4",
     "features": [
         "1.support fq head with @XXXX.XXXX.X",
         "2.automatically skip empty fq files for spades",
     ]},
    {"number": "1.0.3a",
     "features": [
         "1.gunzip",
     ]},
    {"number": "1.0.3",
     "features": [
         "1.accept .gz/.zip files as input",
         "2.logging errors as utf8",
     ]},
    {"number": "1.0.2a",
     "features": [
         "1.prompt failure in Utilities/slim_fastg.py",
     ]},
    {"number": "1.0.2",
     "features": [
         "1.removing duplicates become a parameter to control the memory usage.",
     ]},
    {"number": "1.0.1a",
     "features": [
         "1.Fix the bug of running spades.py with --continue when no output file exists.",
     ]},
    {"number": "1.0.1",
     "features": [
         "1.Add default reference (Library/SeqReference).",
     ]}
]
