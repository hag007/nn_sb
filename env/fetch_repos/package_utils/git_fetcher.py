import sys
import subprocess
import constants
import os
sys.path.insert(0, '../../../')


def git_clone(repository, branch="master"):
    project_name = os.path.basename(repository)
    print "{}_{}".format("git clone", repository)
    print subprocess.Popen("git clone {}".format(repository), shell=True, stdout=subprocess.PIPE, cwd=constants.REPOS_DIR).stdout.read()
    print subprocess.Popen(["git checkout {}".format(branch)], shell=True, stdout=subprocess.PIPE,
                           cwd=os.path.join(constants.REPOS_DIR, project_name)).stdout.read()

def git_pull(project_name, branch="master"):
    project_dir = os.path.join(constants.REPOS_DIR, project_name)
    print subprocess.Popen("git stash", shell=True, stdout=subprocess.PIPE,
                           cwd=project_dir).stdout.read()
    print subprocess.Popen("git pull origin {}".format(branch), shell=True, stdout=subprocess.PIPE, cwd=project_dir).stdout.read()

def git_sync_repo(repository, branch="master"):
    if os.path.exists(os.path.join(constants.REPOS_DIR, os.path.basename(repository))):
        print "about to pull git repo {}".format(os.path.basename(repository))
        git_pull(project_name=os.path.basename(repository), branch="master")
    else:
        print "about to clone git repo {}".format(os.path.basename(repository))
        git_clone(repository, branch="master")


