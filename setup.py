import setuptools
from setuptools import setup
from GetOrganelleLib.versions import get_versions
from GetOrganelleLib.pipe_control_func \
    import remove_db_postfix, check_fasta_seq_names, make_blast_db, build_bowtie2_db, executable
import platform
import sys
import os


# system info
SYSTEM_NAME = ""
if platform.system() == "Linux":
    SYSTEM_NAME = "linux"
elif platform.system() == "Darwin":
    SYSTEM_NAME = "macOS"
else:
    sys.stdout.write("Error: currently GetOrganelle is not supported for " + platform.system() + "! ")
    exit()
# python version
MAJOR_VERSION, MINOR_VERSION = sys.version_info[:2]
if MAJOR_VERSION == 2 and MINOR_VERSION >= 7:
    pass
elif MAJOR_VERSION == 3 and MINOR_VERSION >= 5:
    pass
else:
    sys.stdout.write("Python version have to be 2.7+ or 3.5+")
    sys.exit(0)

sys.stdout.write("Python " + str(sys.version).replace("\n", " ") + "\n")
sys.stdout.write("Using setuptools " + str(setuptools.__version__) + "\n")

# python libs
install_dependencies = []
try:
    import numpy
except ImportError:
    if MAJOR_VERSION == 3:
        install_dependencies.append("numpy>=1.16.4")
    else:
        install_dependencies.append("numpy==1.16.4")
else:
    sys.stdout.write("Existed module numpy " + str(numpy.__version__) + "\n")
try:
    import scipy
except ImportError:
    if MAJOR_VERSION == 3:
        install_dependencies.append("scipy>=1.3.0")
    else:
        # higher version not compatible with python2
        install_dependencies.append("scipy==1.2.1")
else:
    sys.stdout.write("Existed module numpy " + str(scipy.__version__) + "\n")
try:
    import sympy
except ImportError:
    if MAJOR_VERSION == 3:
        install_dependencies.append("sympy>=1.4")
    else:
        install_dependencies.append("sympy==1.4")
else:
    sys.stdout.write("Existed module sympy " + str(sympy.__version__) + "\n")

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
LIB_NAME = "GetOrganelleLib"
LIB_DIR = os.path.join(PATH_OF_THIS_SCRIPT, LIB_NAME)
DEP_NAME = "GetOrganelleDep"
DEP_DIR = os.path.join(PATH_OF_THIS_SCRIPT, DEP_NAME)
NOT_NAME = "LabelDatabase"
NOT_DIR = os.path.join(PATH_OF_THIS_SCRIPT, LIB_NAME, NOT_NAME)
SEQ_NAME = "SeedDatabase"
SEQ_DIR = os.path.join(PATH_OF_THIS_SCRIPT, LIB_NAME, SEQ_NAME)
if "--continue" in sys.argv:
    RESUME = True
    sys.argv.remove("--continue")
else:
    RESUME = False
if "--in-situ" in sys.argv:
    in_situ = True
    sys.argv.remove("--in-situ")
else:
    in_situ = False
if "--keep-index" in sys.argv:
    keep_index = True
    sys.argv.remove("--keep-index")
else:
    keep_index = False


def get_recursive_files(target_dir, start_from="", exclude_files=None):
    if exclude_files is None:
        exclude_files = set()
    assert target_dir.startswith(start_from), "target_dir should be starting with start_from!"
    omit_len = len(start_from.rstrip("/") + "/") if start_from else 0
    for f_dir, sub_dirs, files in os.walk(target_dir):
        for i_file in files:
            if not i_file.startswith(".") and os.path.join(f_dir, i_file)[omit_len:] not in exclude_files:
                yield os.path.join(f_dir, i_file)[omit_len:]


EXCLUDE_SHARE_SPADES_PATHS = set()
# do not exclude other version in case some users install it at root with python2 but use it with python3 ...
# if os.path.exists(DEP_DIR):
#     if MAJOR_VERSION == 2:
#         for ex_f in get_recursive_files(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes/share/spades/joblib3"), DEP_DIR):
#             EXCLUDE_SHARE_SPADES_PATHS.add(ex_f)
#         for ex_f in get_recursive_files(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes/share/spades/pyyaml3"), DEP_DIR):
#             EXCLUDE_SHARE_SPADES_PATHS.add(ex_f)
#     elif MAJOR_VERSION == 3:
#         for ex_f in get_recursive_files(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes/share/spades/joblib2"), DEP_DIR):
#             EXCLUDE_SHARE_SPADES_PATHS.add(ex_f)
#         for ex_f in get_recursive_files(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes/share/spades/pyyaml2"), DEP_DIR):
#             EXCLUDE_SHARE_SPADES_PATHS.add(ex_f)


scripts_to_install = ["get_organelle_from_reads.py",
                      "get_organelle_from_assembly.py",
                      "Utilities/check_annotations.py",
                      "Utilities/cook_coding_for_blast.py",
                      "Utilities/disentangle_organelle_assembly.py",
                      "Utilities/evaluate_assembly_using_mapping.py",
                      "Utilities/fastg_to_gfa.py",
                      "Utilities/get_pair_reads.py",
                      "Utilities/gfa_to_fastg.py",
                      "Utilities/gfa_to_fasta.py",
                      "Utilities/join_spades_fastg_by_blast.py",
                      "Utilities/make_alignment_with_selected_path_when_overlaps.py",
                      "Utilities/make_batch_for_iteratively_mapping_assembling.py",
                      "Utilities/make_batch_for_get_organelle.py",
                      "Utilities/plastome_arch_info.py",
                      "Utilities/rm_low_coverage_duplicated_contigs.py",
                      "Utilities/round_statistics.py",
                      "Utilities/slim_fastg.py",
                      "Utilities/summary_get_organelle_output.py"]
dep_scripts_to_change = []
if os.path.isdir(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "bin")):
    for spades_script in os.listdir(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "bin")):
        if spades_script.endswith(".py"):
            dep_scripts_to_change.append(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "bin", spades_script))
if os.path.exists(os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2", "bowtie2-build")):
    dep_scripts_to_change.append(os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2", "bowtie2-build"))
# rename execution program if not python
if os.path.basename(sys.executable) != "python":
    for rename_py_script in scripts_to_install + dep_scripts_to_change:
        original_lines = open(rename_py_script).readlines()
        original_lines[0] = "#!" + sys.executable + "\n"
        open(rename_py_script, "w").writelines(original_lines)


def initialize_notation_database(which_blast, overwrite=False):
    for fasta_f in os.listdir(NOT_DIR):
        if fasta_f.endswith(".fasta") and fasta_f[:-6] in ("embplant_pt", "other_pt", "embplant_mt", "embplant_nr",
                                                           "animal_mt", "fungus_mt"):
            fasta_f = os.path.join(NOT_DIR, fasta_f)
            output_base = remove_db_postfix(fasta_f)
            sys.stdout.write("makeblastdb " + output_base + " ... ")
            sys.stdout.flush()
            if overwrite or sum([os.path.exists(output_base + postfix) for postfix in (".nhr", ".nin", ".nsq")]) != 3:
                make_blast_db(input_file=fasta_f,
                              output_base=output_base,
                              which_blast=which_blast)
                sys.stdout.write("finished\n")
            else:
                sys.stdout.write("skipped\n")


def initialize_seed_database(which_bowtie2, overwrite=False):
    for fasta_f in os.listdir(SEQ_DIR):
        if fasta_f.endswith(".fasta") and fasta_f[:-6] in ("embplant_pt", "other_pt", "embplant_mt", "embplant_nr",
                                                           "animal_mt", "fungus_mt"):
            fasta_f = os.path.join(SEQ_DIR, fasta_f)
            new_seed_file = fasta_f + ".modified"
            changed = check_fasta_seq_names(fasta_f, new_seed_file)
            if changed:
                seed_file = new_seed_file
            else:
                seed_file = fasta_f
            output_base = remove_db_postfix(fasta_f) + ".index"
            sys.stdout.write("bowtie2-build " + output_base + " ... ")
            sys.stdout.flush()
            if overwrite or sum([os.path.exists(output_base + postfix)
                                 for postfix in
                                 (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")]) != 6:
                build_bowtie2_db(seed_file=seed_file, seed_index_base=output_base, which_bowtie2=which_bowtie2,
                                 overwrite=overwrite, random_seed=12345, silent=True)
                sys.stdout.write("finished\n")
            else:
                sys.stdout.write("skipped\n")
            if changed:
                os.remove(seed_file)


# check BLAST and make blast database
if os.path.exists(os.path.join(DEP_DIR, SYSTEM_NAME, "ncbi-blast")):
    files_to_check = ["blastn", "makeblastdb"]
    for check_f in files_to_check:
        check_file_path = os.path.join(DEP_DIR, SYSTEM_NAME, "ncbi-blast", check_f)
        if not os.path.exists(check_file_path):
            raise EnvironmentError(check_file_path + " not exists!")
        os.chmod(check_file_path, 0o755)
# TODO set overwrite=True
if executable(os.path.join(DEP_DIR, SYSTEM_NAME, "ncbi-blast", "makeblastdb")):
    initialize_notation_database(which_blast=os.path.join(DEP_DIR, SYSTEM_NAME, "ncbi-blast"), overwrite=not RESUME)
elif executable("makeblastdb"):
    initialize_notation_database(which_blast="", overwrite=not RESUME)
else:
    raise EnvironmentError("makeblastdb not found in the $PATH nor in " +
                           os.path.join(DEP_DIR, SYSTEM_NAME, "ncbi-blast") + "!\n"
                           "change directory to GetOrganelle and git clone git://github.com/Kinggerm/GetOrganelleDep\n"
                           "or get BLAST via http://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/LATEST")

# check Bowtie2 and build seed index
if os.path.exists(os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2")):
    files_to_check = ["bowtie2", "bowtie2-align-l", "bowtie2-build", "bowtie2-build-l"]
    for check_f in files_to_check:
        check_file_path = os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2", check_f)
        if not os.path.exists(check_file_path):
            raise EnvironmentError(check_file_path + " not exists!")
        os.chmod(check_file_path, 0o755)
# TODO set overwrite=True
if executable(os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2", "bowtie2-build")):
    initialize_seed_database(os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2"), overwrite=not RESUME)
elif executable("bowtie2-build"):
    initialize_seed_database("", overwrite=not RESUME)
else:
    raise EnvironmentError("bowtie2-build not found in the $PATH nor in " +
                           os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2") + "!\n"
                           "change directory to GetOrganelle and git clone git://github.com/Kinggerm/GetOrganelleDep\n"
                           "or get Bowtie2 via http://bowtie-bio.sourceforge.net/bowtie2/index.shtml")

# check SPAdes
if os.path.exists(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes")):
    files_to_check = ["metaspades.py", "spades-bwa", "spades-gbuilder", "spades-ionhammer", "spades.py",
                      "plasmidspades.py", "spades-core", "spades-gmapper", "spades-kmercount", "spades_init.py",
                      "rnaspades.py", "spades-corrector-core", "spades-hammer", "spades-truseq-scfcorrection",
                      "truspades.py"]
    for check_f in files_to_check:
        check_file_path = os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "bin", check_f)
        if not os.path.exists(check_file_path):
            raise EnvironmentError(check_file_path + " not exists!")
        os.chmod(check_file_path, 0o755)
    files_to_check = ["configs", "spades_pipeline"]
    for check_f in files_to_check:
        check_file_path = os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "share", "spades", check_f)
        if not os.path.exists(check_file_path):
            raise EnvironmentError(check_file_path + " not exists!")
if executable(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "bin", "spades.py")):
    pass
elif executable("spades.py"):
    pass
else:
    raise EnvironmentError("spades.py not found in the $PATH nor in " +
                           os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "bin") + "!\n"
                           "change directory to GetOrganelle and git clone git://github.com/Kinggerm/GetOrganelleDep\n"
                           "or get SPAdes via http://cab.spbu.ru/software/spades/")


if not in_situ:
    setup(
        name="GetOrganelle",
        version=get_versions(),
        description="a fast and versatile toolkit for accurate de novo assembly of organelle genomes.",
        author="Jian-Jun Jin",
        author_email="jinjianjun@mail.kib.ac.cn",
        url="http://github.org/Kinggerm/GetOrganelle",
        license="GNU General Public License, version 3",
        packages=[LIB_NAME, DEP_NAME],
        platforms="linux/MacOS",
        scripts=scripts_to_install,
        # relative path to each package
        package_data={LIB_NAME: [os.path.join(NOT_NAME, "*.n*"),
                                 os.path.join(SEQ_NAME, "*.bt2l"),
                                 os.path.join(SEQ_NAME, "*.fasta")],
                      DEP_NAME: [this_file
                                 for this_file in
                                 get_recursive_files(target_dir=os.path.join(DEP_DIR, SYSTEM_NAME),
                                                     start_from=DEP_DIR, exclude_files=EXCLUDE_SHARE_SPADES_PATHS)]},
        install_requires=install_dependencies,
        zip_safe=False
        )
    if keep_index:
        for temp_dir_or_files in ("build", "dist", "*.pyc", "*.tgz", "*.egg-info"):
            os.system("rm -vrf " + str(os.path.join(PATH_OF_THIS_SCRIPT, temp_dir_or_files)))
    else:
        for temp_dir_or_files in ("build", "dist", "*.pyc", "*.tgz", "*.egg-info",
                                  os.path.join(LIB_NAME, NOT_NAME, "*.n*"), os.path.join(LIB_NAME, SEQ_NAME, "*.bt2l")):
            os.system("rm -vrf " + str(os.path.join(PATH_OF_THIS_SCRIPT, temp_dir_or_files)))
else:
    for script_chmod in scripts_to_install:
        os.chmod(script_chmod, 0o755)
