#!/usr/bin/env python
"""This script converts a gfa (Graphical Fragment Assembly) file into a fastg file"""
import sys
import os



def main():
    if len(sys.argv) > 1:
        for i in sys.argv:
            if '-h' in i or 'help' in i:
                print("Usage: gfa2fastg.py *.gfa")
                break
        else:
            path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
            sys.path.insert(0, os.path.join(path_of_this_script, ".."))
            from GetOrganelleLib.assembly_parser import Assembly
            for gfa_file in sys.argv[1:]:
                Assembly(gfa_file).write_to_fastg(gfa_file + '.fastg', rename_if_needed=True, echo_rename_warning=True)
    else:
        if type(2/1) == float:
            gfa_file = input('Please input gfa file:').strip()
        else:
            gfa_file = raw_input('Please input gfa file:').strip()
        path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
        sys.path.insert(0, os.path.join(path_of_this_script, ".."))
        from GetOrganelleLib.assembly_parser import Assembly
        if gfa_file.strip():
            Assembly(gfa_file).write_to_fastg(gfa_file + '.fastg', rename_if_needed=True, echo_rename_warning=True)


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