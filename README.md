# Introduction

This library of Python modules was intended to provide utilities that could be used for many purposes rather than a specific project. It can be regarded as an SVN external.


# Linking and Updating
## adding to local repository
% git submodule add git@lex-gitlab.aer.com:RC/common_modules.git externals/common

% git submodule update --init --recursive

## updating local repository to what is in central repository
% cd externals/common

% git pull

NOTE: this is not done automatically whenever a change to the common_modules.git is pushed (e.g., when one is in projectA and adds common_modules as externals/common, a simple "git pull" in projectA will not update externals/common). There are justifications for this protocol, but the automatic update will be added in the future and may become the default: https://stackoverflow.com/questions/1899792/why-is-git-submodule-update-not-automatic-on-git-checkout

## committing local changes to central repository
% cd externals/common

% touch temp

% git add temp

% git commit temp -m 'temp added for testing'

% git push origin master

% git commit -a -m 'removed temp test file'

% git push origin master
