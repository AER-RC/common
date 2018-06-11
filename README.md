# Introduction

Much like https://lex-gitlab.aer.com/RC/common_modules

# Adding as a submodule to an existing repository

Assuming that the user is already in the directory that hosts the repository in which these modules need to be added as submodules, they can add this library with:

```
git submodule add git@github.com:pernak18/common.git common
git submodule update --init --recursive
```

The directory `common` is then added as a subdirectory, and a snapshot of the common_modules library is populated into it such that `common` is associated with a specific commit to pernak18/common. It is important to note that the library is in a "detached HEAD" state, meaning that there is no branch association. If one makes changes to any file in `common` and wants to commit those changes to the common repository, they will have to:

```
git checkout master
```

in `common` to get onto the master branch of pernak18/common, then they will be able to push as they normally would.

## updating local repository to what is in central repository

Because `common` in this case is in a detached HEAD state, updates have to be explicitly requested:

```
cd common
git pull
```

NOTE: this is not done automatically whenever a change to the pernak18/common.git remote repository is pushed (e.g., when one is in another project directory and adds pernak18/common as common, a simple "git pull" in the other project directory will not update `common` in this project). There are justifications for this protocol, but the automatic update will be added in the future and may become the default: https://stackoverflow.com/questions/1899792/why-is-git-submodule-update-not-automatic-on-git-checkout.

