# /usr/bin/env python


def get_versions():
    return versions[0]["number"]


versions = [
    {"number": "1.0.1a",
     "features": [
         "Fix the bug of running spades.py with --continue when no output file exists.",
     ]},
    {"number": "1.0.1",
     "features": [
         "Add default reference (Library/SeqReference).",
     ]}
]
