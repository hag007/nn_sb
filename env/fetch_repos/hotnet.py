import constants
from package_utils.git_fetcher import git_sync_repo
GIT_REPO = "https://github.com/raphael-group/hotnet"

def sync_repo():
    git_sync_repo(GIT_REPO)

