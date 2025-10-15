# Git usage guidelines
This document provides a number of rules and guidelines for working with git
in the DREAM repository.


## Do not add all files
To avoid cluttering the DREAM repository with hidden and other unnecessary
files, please avoid commands such as ``git add -a`` and ``git add .``. Instead,
we recommend that you use
```bash
git add -u
```
which adds only modified files which have previously been tracked.
Alternatively, and in the case of brand new files, you may add specific files
individually.

## Work on own branch or fork
Developers who are granted write access to the repository are encouraged to
make their changes on a dedicated branch in the repository. Everyone with write
access can push such branches to the repository freely. If you do not have write
access to the repository but would still like to contribute, you are encouraged
to fork the repository.

Once you deem that your is mature and ready for others to use, you are
encouraged to make a pull request to the master branch on the official DREAM
repository. The pull request will be reviewed by members of the DREAM Developer
Council. Once everyone is satisfied with the branch, it will be merged into the
master branch by a maintainer.

## Do not commit large and/or binary files
Binary files are generally ill-suited for a git repository, as are large files
(> 1 MiB) of any kind. Exceptions can be made in specific cases, but please
reach out to the developers in case you feel that an exception to this rule
is warranted.

## Never store sensitive data in the DREAM directory
Data which should not be made public should never be stored inside the DREAM
directory, even if it is just for a temporary period and you have no intention
of committing it to the repository.

