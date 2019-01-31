import constants
import hotnet2
import bionet
from package_utils.git_fetcher import git_sync_repo


def sync_repos():
    hotnet2.sync_package()
    bionet.sync_package()
