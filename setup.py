# coding utf8
import setuptools
from setuptools import setup
from GetOrganelleLib.versions import get_versions
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
sys.stdout.write("PLATFORM: " + " ".join(platform.uname()) + "\n")
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
try:
    import requests
except ImportError:
    install_dependencies.append("requests[security]")
else:
    sys.stdout.write("Existed module requests " + str(requests.__version__) + "\n")

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
LIB_NAME = "GetOrganelleLib"
# LIB_DIR = os.path.join(PATH_OF_THIS_SCRIPT, LIB_NAME)
DEP_NAME = "GetOrganelleDep"
DEP_DIR = os.path.join(PATH_OF_THIS_SCRIPT, DEP_NAME)
LBL_NAME = "LabelDatabase"
# LBL_DIR = os.path.join(PATH_OF_THIS_SCRIPT, LIB_NAME, LBL_NAME)
SEQ_NAME = "SeedDatabase"
# SEQ_DIR = os.path.join(PATH_OF_THIS_SCRIPT, LIB_NAME, SEQ_NAME)
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
if "--keep-temp" in sys.argv:
    keep_temp = True
    sys.argv.remove("--keep-temp")
else:
    keep_temp = False


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

scripts_to_install = ["get_organelle_from_reads.py",
                      "get_organelle_from_assembly.py",
                      "Utilities/check_annotations.py",
                      "Utilities/cook_coding_for_blast.py",
                      "Utilities/disentangle_organelle_assembly.py",
                      "Utilities/evaluate_assembly_using_mapping.py",
                      "Utilities/fastg_to_gfa.py",
                      "Utilities/get_organelle_config.py",
                      "Utilities/get_pair_reads.py",
                      "Utilities/gfa_to_fastg.py",
                      "Utilities/gfa_to_fasta.py",
                      "Utilities/join_spades_fastg_by_blast.py",
                      "Utilities/make_batch_for_iteratively_mapping_assembling.py",
                      "Utilities/make_batch_for_get_organelle.py",
                      "Utilities/plastome_arch_info.py",
                      "Utilities/rm_low_coverage_duplicated_contigs.py",
                      "Utilities/round_statistics.py",
                      "Utilities/slim_graph.py",
                      "Utilities/summary_get_organelle_output.py",
                      "Utilities/reconstruct_graph_from_fasta.py"]
# rename execution program if not python
dep_scripts_to_change = []
if os.path.isdir(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "bin")):
    for spades_script in os.listdir(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "bin")):
        if spades_script.endswith(".py") and not spades_script.startswith("."):
            dep_scripts_to_change.append(os.path.join(DEP_DIR, SYSTEM_NAME, "SPAdes", "bin", spades_script))
if os.path.exists(os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2", "bowtie2-build")):
    dep_scripts_to_change.append(os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2", "bowtie2-build"))
if os.path.basename(sys.executable) != "python":
    for rename_py_script in scripts_to_install + dep_scripts_to_change:
        original_lines = open(rename_py_script, encoding="utf-8").readlines()
        original_lines[0] = "#!" + sys.executable + "\n"
        open(rename_py_script, "w", encoding="utf-8").writelines(original_lines)


# check local BLAST
if os.path.exists(os.path.join(DEP_DIR, SYSTEM_NAME, "ncbi-blast")):
    files_to_check = ["blastn", "makeblastdb"]
    for check_f in files_to_check:
        check_file_path = os.path.join(DEP_DIR, SYSTEM_NAME, "ncbi-blast", check_f)
        if not os.path.exists(check_file_path):
            raise EnvironmentError(check_file_path + " not exists!")
        os.chmod(check_file_path, 0o755)

# check local Bowtie2
if os.path.exists(os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2")):
    files_to_check = ["bowtie2", "bowtie2-align-l", "bowtie2-build", "bowtie2-build-l"]
    for check_f in files_to_check:
        check_file_path = os.path.join(DEP_DIR, SYSTEM_NAME, "bowtie2", check_f)
        if not os.path.exists(check_file_path):
            raise EnvironmentError(check_file_path + " not exists!")
        os.chmod(check_file_path, 0o755)

# check local SPAdes
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


PACKAGES = [LIB_NAME]
PACKAGE_DATA = {}
# PACKAGE_DATA = {LIB_NAME: [os.path.join(LBL_NAME, "VERSION"),
#                            os.path.join(SEQ_NAME, "VERSION")]}
if os.path.isdir(DEP_DIR) and os.path.isfile(os.path.join(DEP_DIR, "__init__.py")):
    PACKAGES.append(DEP_NAME)
    PACKAGE_DATA[DEP_NAME] = [this_file
                              for this_file in
                              get_recursive_files(target_dir=os.path.join(DEP_DIR, SYSTEM_NAME),
                                                  start_from=DEP_DIR, exclude_files=EXCLUDE_SHARE_SPADES_PATHS)]


if not in_situ:
    setup(
        name="GetOrganelle",
        version=get_versions(),
        description="a fast and versatile toolkit for accurate de novo assembly of organelle genomes.",
        author="Jian-Jun Jin",
        author_email="jinjianjun@mail.kib.ac.cn",
        url="http://github.com/Kinggerm/GetOrganelle",
        license="GNU General Public License, version 3",
        packages=PACKAGES,
        platforms="linux/MacOS",
        scripts=scripts_to_install,
        # relative path to each package
        package_data=PACKAGE_DATA,
        install_requires=install_dependencies,
        zip_safe=False
        )
    if keep_temp:
        for temp_dir_or_files in ("build", "dist", "*.pyc", "*.tgz", "*.egg-info"):
            os.system("rm -vrf " + str(os.path.join(PATH_OF_THIS_SCRIPT, temp_dir_or_files)))
    else:
        for temp_dir_or_files in ("build", "dist", "*.pyc", "*.tgz", "*.egg-info",
                                  os.path.join(LIB_NAME, LBL_NAME, "*.n*"), os.path.join(LIB_NAME, SEQ_NAME, "*.bt2l")):
            os.system("rm -vrf " + str(os.path.join(PATH_OF_THIS_SCRIPT, temp_dir_or_files)))
else:
    for script_chmod in scripts_to_install:
        os.chmod(script_chmod, 0o755)
