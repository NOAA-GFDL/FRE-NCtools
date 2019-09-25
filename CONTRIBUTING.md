code># Contributing to FRE NCtools

Thank you for taking time to contribute.

FRE-NCtools is a collection of Fortran, C and POSIX shell scripts to query,
create and manipulate [netCDF](https://www.unidata.ucar.edu/software/netcdf/)
files for use with the [Geophysical Fluid Dynamics Laboratory](https://www.gfdl.noaa.gov)
(GFDL) Climate models run using the [Flexible Modeling System](https://www.gfdl.noaa.gov/fms)
(FMS) framework.  FRE-NCtools is developed internally at GFDL, but is made
available on GitHub in the [NOAA-GFDL](https://github.com/NOAA-GFDL)
organization.

What follows is set of guidelines and how-tos for contributing to FRE-NCtools.  
These are guidelines, not rules.  Use your best judgement and feel free to
propose changes to this document in a pull request.

## Table of Contents

* [Code of Conduct](#code-of-conduct)
* [Getting Support](#getting-support)
* [How Can I Contribute](#how-can-i-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Pull Requests](#pull-requests)
  * [Tests](#tests)
* [Styleguides](#styleguides)
* [Maintainer How-To](#maintainer-how-to)
  * [Migrate GitHub Pull Request](#migrate-github-pull-request)
  * [GitLab Push Mirror](#gitlab-push-mirror)

## Code of Conduct

This project and everyone participating in it is governed by the
[Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to
uphold this code. Please report unacceptable behavior to
[gfdl.climate.model.info@noaa.gov](mailto:gfdl.climate.model.info@noaa.gov).

## Getting Support

To get support, open a GitHub or GitLab issue using the support issue template.  
Please be as descriptive as possible.

The members of the [Modeling Systems](https://www.gfdl.noaa.gov/modeling-systems)
group at GFDL are the main developers and maintainers of the FRE-NCtools
package and applications.  Modeling Systems is a small group, and our ability to
offer support for these tools is limited.  We request your patience when
submitting a support issue.

## How Can I Contribute

### Reporting Bugs

This section guides you through submitting a bug report for FRE-NCtools.
Following these guidelines will help the maintainers and community understand
your report, reproduce the behavior and find related reports.

Before creating a new issue, please perform a [cursory search](https://github.com/search?q=+is%3Aissue+repo%3ANOAA-GFDL%2FFRE-NCtools&type=issue) to see if the
problem has already been reported.  If it has and *the issue is still open*, add
a comment to the existing issue instead of opening a new issue.

Bugs are tracked as [GitHub](https://help.github.com/en/articles/about-issues)
and [GitLab](https://docs.gitlab.com/ee/user/project/issues/) issues.  Where you
create an issue depends on if you have a GFDL account.  Those with GFDL accounts
can create issues in the [GFDL GitLab FRE-NCtools project](https://gitlab.gfdl.noaa.gov/fre/fre-nctools).
However, to better support the community, it may be better to report the issue
on [GitHub](https://github.com/NOAA-GFDL/FRE-NCtools).

Please use the bug report template when creating the issue.  Be sure to include
as much information as possible to describe the issue, and **use a clear and
descriptive** title.

### Suggesting Enhancements

This section guides you through submitting a enhancement suggestion for FRE-NCtools.
Following these guidelines will help the maintainers and community understand
your suggestion and find related suggestions.

Before creating an enhancement suggestion, please perform a
[cursory search](https://github.com/search?q=+is%3Aissue+repo%3ANOAA-GFDL%2FFRE-NCtools&type=issue)
as you may find out that you don't need to create one.  When creating an enhancement
suggetion, please include as many details as possible.  Fill in the feature
request template, and **use a clear and descriptive** title.

### Pull Requests

FRE-NCtools, on GitHub, follows the [Fork/Pull Request workflow](https://guides.github.com/activities/forking/).  When submitting a Pull Request, follow these steps:

1. Follow all instructions in the template
1. Follow the [styleguides](#styleguides)
1. After submitting your pull request, verify all [status checks](https://help.github.com/articles/about-status-checks/) are passing <details><summary>What if the status checks are
failing?</summary>If a status check is failing, and you believe it is unrelated
to your change, please leave a comment on the pull request why you believe the
failure is unrelated.  A maintainer will re-run the status check.</details>

### Tests

FRE-NCtools has a small number of tests that are run during the `make check`
step.  These tests are written using shell commands and the
[Bash Automated Testing System](https://github.com/sstephenson/bats)(BATS).  The
use of BATS is not a requirement when writing new tests, and the use of BATS
may be removed in the future.

All new features and bug fixes should, when possible, include a test.  If a
test is not included with a pull request, the pull request should explain the
reason for not including the test.

## Styleguides

FRE-NCtools does not have an official styleguide.  However, we request work in
a tool follows the same style already used.  We will accept pull requests that
fully change the style.  However, we request all source for that tool be updates
to use the same style, and the pull request *only* contain the style change, and
not other modifications.

Keep the following in mind when contributing to FRE-NCtools

### Git Commit Messages

* Use the present tense ("Add feature" not "Added feature")
* Limit the first line to 72 characters or less
* Describe the changes in the commit body
* When only changing documentation, include `[ci skip]` in the commit title

### General Language Styleguide

* Use meaningful names for variables, constants and functions
* Use good indentation
* Avoid the use of tab character for indentation
* Use whitespace to help break up the code into relevant blocks
* Use consistent style
* Avoid code duplication
* #define constants must be in uppercase
* Avoid source lines exceeding 80 characters
* Use comments to help describe code
* Clean trailing whitespace

### C Styleguide

* If variables have the same type, declare them on the same line if possible
* Leave one blank line between variable declarations and the first line of code in a function
* Opening brace should be on the same line as a conditional or function
* Avoid global variables where they are unnecessary

### Fortran Styleguide

* List only one module per `use <module>` statement
* Do not use tab characters (this is not supported in the Fortran standard)
* If using newer Fortran standard features, protect the feature using preprocessor macros, and include a test in `configure.ac` to check for that feature
* Do not use `;` to place multiple commands on a single line

### Shell Styleguide

* Use POSIX shell
* Avoid Bash specific extensions
* All error messages should go to STDERR
* Use `$(command)` instead of backticks
* Always check return values
* Favor builtin commands over external commands
* Executable commands should have no extension
* Libraries must have a `.sh` extension

### Perl Styleguide

* Follow [perlstyle](https://perldoc.perl.org/perlstyle.html) when possible

### Python Styleguide

* Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) when possible

## Maintainer How-To

The following are *How To* guides to help maintain the FRE-NCtools repository.

### Migrate GitHub Pull Request

When a pull request is submitted on GitHub and the GitHub status checks are
passing, the maintainer will need to migrate the pull request to the internal
[GFDL GitLab FRE-NCtools](https://gitlab.gfdl.noaa.gov/fre/fre-nctools)
repository.  The following steps can be taken to migrate the pull request.

1. <details><summary>Review the GitHub code to ensure no malicious code or scripts were introduced</summary>This is required as migrating the pull request to GitLab will automatically run the GitLab tests.</details>
1. <details><summary>Clone the GitLab repository</summary><code>git clone git@gitlab.gfdl.noaa.gov/fre/fre-nctools.git</code></details>
1. <details><summary>Add the GitHub repository as a remote</summary><code>git remote add github https://github.com/NOAA-GFDL/FRE-NCtools.git</code></details>
1. Find the ID number of the GitHub pull request.  This is the sequence of digits right after the pull request's title.
1. <details><summary>Fetch the reference pull request based on its ID number, and place it in a new branch</summary><code>git fetch github pull/<ID>/head:<BRANCH_NAME></code>  Replace <code>&lt;ID&gt;</code> with the pull request ID, and <code>&lt;BRANCH_NAME&gt;</code> with the new branch name.</details>
1. <details><summary>Switch to the new branch</summary><code>git checkout &lt;BRANCH_NAME&gt;</code></details>
1. <details><summary>Push the new branch to GitLab</summary><code>git push origin <BRANCH_NAME></code></details>

Once the GitHub pull request is migrated, open a GitLab merge request, and assign
a reviewer to review the code, and perform the merge once the GitLab.

### GitLab Push Mirror

The internal [GFDL GitLab FRE-NCtools](https://gitlab.gfdl.noaa.gov/fre/fre-nctools)
repository contains the *official* development trunk for FRE-NCtools.  To help
ensure the GitHub mirror remains up-to-date with development and releases, the
GitLab protected branches are automatically mirrored to GitHub using GitLab's
push [repository mirroring](https://docs.gitlab.com/ce/workflow/repository_mirroring.html#setting-up-a-push-mirror-from-gitlab-to-github-core).
