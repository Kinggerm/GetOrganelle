#!/usr/bin/env python
"""This script converts a fastg file into a fasta file"""
import sys
import os


def main():
    if len(sys.argv) > 1:
        for i in sys.argv:
            if '-h' in i or 'help' in i:
                print("Usage: fastg2gfa.py *.fastg")
                break
        else:
            PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
            sys.path.insert(0, os.path.join(PATH_OF_THIS_SCRIPT, ".."))
            from GetOrganelleLib.assembly_parser import Assembly
            for fastg in sys.argv[1:]:
                this_assembly = Assembly(fastg)
                this_assembly.write_to_gfa(fastg + ".gfa")
    else:
        if type(2/1) == float:
            fastg = input('Please input gfa file:').strip()
        else:
            fastg = raw_input('Please input gfa file:').strip() # type: ignore
        PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
        sys.path.insert(0, os.path.join(PATH_OF_THIS_SCRIPT, ".."))
        from GetOrganelleLib.assembly_parser import Assembly
        if fastg.strip():
            this_assembly = Assembly(fastg)
            this_assembly.write_to_gfa(fastg + ".gfa")


if __name__ == '__main__':
    main()


"""Copyright 2016 Jianjun Jin

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License."""