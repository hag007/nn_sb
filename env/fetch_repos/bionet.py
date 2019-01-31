from package_utils.r_fetcher import fetch_r_package


PROJECT_NAME = "bionet"

def sync_package():
    fetch_r_package(PROJECT_NAME)

