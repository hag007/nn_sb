import shutil
import os
import json
import sys
sys.path.insert(0, '../../')
import constants

dir_path = os.path.dirname(os.path.realpath(__file__))

def sync_from_project():
    dest = constants.BASE_PROFILE
    source = dir_path
    replace_folder(source, dest)


def sync_to_project():
    dest = dir_path
    source = constants.BASE_PROFILE
    replace_folder(source, dest)

def replace_folder(source, dest):
    if os.path.exists(dest):
        print "delete {}".format(os.path.join(dest, "dictionaries"))
        shutil.rmtree(os.path.join(dest, "dictionaries"))
        print "delete {}".format(os.path.join(dest, "list"))
        shutil.rmtree(os.path.join(dest, "list"))
    else:
        print "create {}".format(dest)
        os.makedirs(dest)
        
    shutil.copytree(os.path.join(source,"dictionaries"), os.path.join(dest, "dictionaries"))
    shutil.copytree(os.path.join(source,"list"), os.path.join(dest, "list"))

