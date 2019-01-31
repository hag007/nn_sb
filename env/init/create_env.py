import sys
sys.path.insert(0, '../../')
import shutil
import os
import json
from utils import download_resources
from env.fetch_repos import runner

import constants

dir_path = os.path.dirname(os.path.realpath(__file__))

def main():
    dest = constants.BASE_PROFILE
    if not os.path.exists(os.path.join(dest, "output")):
        os.makedirs(os.path.join(dest, "output"))
    if not os.path.exists(os.path.join(dest, "repos")):
        os.makedirs(os.path.join(dest, "repos"))
    if not os.path.exists(constants.GO_DIR):
        os.makedirs(constants.GO_DIR)
    if not os.path.exists(constants.CACHE_GLOBAL_DIR):
        os.makedirs(constants.CACHE_GLOBAL_DIR)
    #runner.sync_repos()
    # download_resources.main()


if __name__ == "__main__":
    main()
