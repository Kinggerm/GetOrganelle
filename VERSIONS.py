# /usr/bin/env python


def get_versions():
    return versions[0]["number"]


versions = [
    {"number": "1.2.0",
     "features": [
         "more robust and precise in disentangling graph: "
         "report contamination; "
         "detect parallel contigs and generate consensus;"
         "estimate chloroplast coverage distribution pattern using weighted GMM with EM and BIC;",
         "re-organize codes",
     ]},
    {"number": "1.1.0d",
     "features": [
         "more precise in disentangling graph.",
         "--prefix added"
     ]},
    {"number": "1.1.0c",
     "features": [
         "--max-reads: default=10,000,000 for cp and nr, default=50,000,000 for mt",
     ]},
    {"number": "1.1.0b",
     "features": [
         "'-w' could be followed with a ratio of word_size to average_read_length now",
     ]},
    {"number": "1.1.0a",
     "features": [
         "Add options-maximum_number_of_reads with default=10,000,000 to avoid unnecessary overloading",
     ]},
    {"number": "1.1.0",
     "features": [
         "automatically exporting final sequence(s)",
         "adding disentangle_organelle_assembly.py",
         "adding assembly_parser.py"
         "re-organize codes",
         "revise README.md",
     ]},
    {"number": "1.0.5",
     "features": [
         "re-organize codes",
     ]},
    {"number": "1.0.4",
     "features": [
         "support fq head with @XXXX.XXXX.X",
         "automatically skip empty fq files for spades",
     ]},
    {"number": "1.0.3a",
     "features": [
         "gunzip",
     ]},
    {"number": "1.0.3",
     "features": [
         "accept .gz/.zip files as input",
         "logging errors as utf8",
     ]},
    {"number": "1.0.2a",
     "features": [
         "prompt failure in Utilities/slim_fastg.py",
     ]},
    {"number": "1.0.2",
     "features": [
         "removing duplicates become a parameter to control the memory usage.",
     ]},
    {"number": "1.0.1a",
     "features": [
         "Fix the bug of running spades.py with --continue when no output file exists.",
     ]},
    {"number": "1.0.1",
     "features": [
         "Add default reference (Library/SeqReference).",
     ]}
]
