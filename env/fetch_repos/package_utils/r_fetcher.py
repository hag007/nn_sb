import sys
import subprocess
import constants
import os
sys.path.insert(0, '../../../')


def fetch_r_package(project_name):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    print subprocess.Popen("pwd".format(project_name), shell=True, stdout=subprocess.PIPE, cwd=dir_path).stdout.read()
    print subprocess.Popen("Rscript ../{}.r".format(project_name), shell=True, stdout=subprocess.PIPE, cwd=dir_path).stdout.read()

