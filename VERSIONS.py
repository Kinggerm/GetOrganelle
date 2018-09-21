# /usr/bin/env python


def get_versions():
    return versions[0]["number"]


versions = [
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
