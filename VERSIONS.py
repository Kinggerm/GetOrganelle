# /usr/bin/env python


def get_versions():
    return versions[0]["number"]


versions = [
    {"number": "1.5.1c",
     "features": [
         "1. get_organelle_reads.py & disentangle_organelle_assembly.py: --random-seed added with default value 12345",
         "2. evaluate_assembly_using_mapping.py: customized error rate info added",
         "3. get_organelle_reads.py: change default & suggested kmers",
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
         "cp -> plant_cp; mt -> plant_mt; nr -> plant_nr; for adding animals",
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
