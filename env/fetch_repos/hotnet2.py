import constants
from package_utils.git_fetcher import git_sync_repo
GIT_REPO = "https://github.com/raphael-group/hotnet2"

def sync_package():
    git_sync_repo(GIT_REPO)

