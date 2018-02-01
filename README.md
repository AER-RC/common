# Introduction

This library of Python modules was intended to provide utilities that could be used for many purposes rather than a specific project. It can be regarded as an SVN external.

There is enough utility in this library that it can be checked out by itself or as a submodule for another repo.

# Standalone Clone
To pull the library by itself:

```
git clone git@lex-gitlab.aer.com:RC/common_modules.git RC_common
```

This command will generate an `RC_common` subdirectory in the user's working directory. If the path to RC_common is included in the user's Python path (`sys.path.append(FULL_PATH_TO_RC_common)`), any module can then be imported. Currently there is only one branch -- master -- and all of the typical Git commands can be used:

```
git pull # pull any commits from remote repository (i.e., repo saved in Gitlab)
git add FILE # add file not in version control (FILE) to repository
git commit FILE -m 'log for changes to FILE' # commit to local repository
git push origin master # push any changes that were committed to local repository to remote repository
```

# Submodule Linking and Updating

There are a couple of ways to access this library as a submodule of another project, which we'll call `RC_Project`. 

## Adding as a submodule to an existing repository

Assuming that the user is already in the RC_Project directory, they can add this library with:

```
git submodule add git@lex-gitlab.aer.com:RC/common_modules.git externals/common
git submodule update --init --recursive
```

The directory `externals/common` is then added as a subdirectory, and a snapshot of the common_modules library is populated into it such that externals/common is associated *with a specific commit to common_modules*. It is important to note that the library is in a "detached HEAD" state, meaning that there is no branch association. If one makes changes to any file in `RC_Project/externals/common` and wants to commit those changes to the common_modules repository, they will have to:

```
git checkout master
```

in `externals/common` to get onto the master branch of common_modules, then they will be able to push as they normally would.

## updating local repository to what is in central repository
Because common_modules in this case is in a detached HEAD state, updates have to be explicitly requested:

```
cd externals/common
git pull
```

NOTE: this is not done automatically whenever a change to the common_modules.git remote repository is pushed (e.g., when one is in another RC_project2 and adds common_modules as externals/common, a simple "git pull" in RC_project2 will not update externals/common). There are justifications for this protocol, but the automatic update will be added in the future and may become the default: https://stackoverflow.com/questions/1899792/why-is-git-submodule-update-not-automatic-on-git-checkout

## committing local changes to central repository
For this library, there is only one branch -- master. To push any changes that one might think should be used outside of RC_Project:

```
cd externals/common
git checkout master
touch temp
git add temp
git commit temp -m 'temp added for testing'
git push origin master
git commit -a -m 'removed temp test file'
git push origin master
```

Should there be future branches, "master" in the previous code block should be replaced with the branch name.

## Cloning an existing repository with submodules

If common_modules has already been added to RC_Project, the clone of RC_Project can be done with one line:

```
git clone --recursive git@lex-gitlab.aer.com:RC/RC_Project.git RC_Project
```

This will generate an `RC_Project` directory that will contain common_modules in whatever subdirectory it was designated as when it was added as a submodule (e.g., `externals/common`). Again, common_modules will be in a detached HEAD state, so be sure to checkout the appropriate branch (probably master) if there is any intent on committing common_modules modifications.

# Submodule Revisions
Let us assume in the next couple of sections that a modification was made to `RC_File`, which already exists in the remote (Gitlab) repository.

## Working Copy vs. Local Repository

```
git diff RC_File
```

The output will be the standard output from a call to the Unix `diff` command.

I prefer the `tkdiff` approach, which is just a graphical "diff" that nicely visualizes (e.g., color codes) the differences between two files. Assuming this (or something similar) is installed on the user's system, they can set the difftool setting:

```
git config --global diff.tool tkdiff
git config --global merge.tool tkdiff
git config --global --add difftool.prompt false
git config --list | grep diff
  diff.tool=tkdiff
  difftool.prompt=false
```

Then we can use the command:

```
git difftool RC_File
```

If called like this verbatim, the HEAD version is in the left tkdiff panel, and the working copy is in the right panel.

## Local Repository vs. Central Repository

Should the simple change to the working copy of RC_File be committed into the user's local repository, the differences can be compared with what is in the central repository in Gitlab with:

```
% git difftool HEAD...origin/master
```

In tkdiff, the most recent commit in the central repo is on the right, and my local repo copy (HEAD) is on the left. these comparisons assume we just wanna compare the HEAD copies with the master branch -- origin/master could presumably just be changed to origin/whateverBranch if branches for the repo exist (which they don't, so i have not tested).

## Updating the Submodules

To update the submodule to the current state of the common_modules repository (assuming we are in `RC_Project`, which has had common_modules designated as a submodule):

```
% git submodule foreach git pull origin RC_Project
```

