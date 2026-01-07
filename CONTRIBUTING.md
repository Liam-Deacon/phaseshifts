# Contributing to the phaseshifts project

A big welcome and thank you for considering contributing to this open source project! We all stand on the shoulders of giants - and with each new contribution we all stand a little taller.

Reading and following these guidelines will help us make the contribution process easy and effective for everyone involved. It also communicates that you agree to respect the time of the developers managing and developing these open source projects. In return, we will reciprocate that respect by addressing your issue, assessing changes, and helping you finalize your pull requests.

## Quick Links

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
  - [Issues](#issues)
  - [Pull Requests](#pull-requests)
- [Getting Help](#getting-help)

## Code of Conduct

We take our open-source community seriously and hold ourselves and other
contributors to high standards of communication. By participating and
contributing to this project, you agree to uphold our
[Code of Conduct][code-of-conduct].

## Getting Started

Contributions are made to this repo via Issues and Pull Requests (PRs). A few general guidelines that cover both:

- Search for existing Issues and PRs before creating your own.
- The maintainer(s) try to ensure issues are handled in a timely manner but, depending on the impact, it could take a while to investigate the root cause. A friendly ping in the comment thread to the submitter or a contributor can help draw attention if your issue is blocking.

### Issues

Issues should be used to report problems with the library, request a new feature, or to discuss potential changes before a PR is created. When you create a new Issue, a template will be loaded that will guide you through collecting and providing the information we need to investigate.

If you find an Issue that addresses the problem you're having, please add your own reproduction information to the existing issue rather than creating a new one. Adding a [reaction](https://github.blog/2016-03-10-add-reactions-to-pull-requests-issues-and-comments/) can also help be indicating to our maintainers that a particular problem is affecting more than just the reporter.

### Pull Requests

PRs are always welcome and can be a quick way to get your fix or improvement slated for the next release. In general, PRs should:

- Only fix/add the functionality in question **OR** address wide-spread whitespace/style issues, not both.
- Add unit or integration tests for fixed or changed functionality.
- Address a single concern in the least number of changed lines as possible.
- Include documentation in the repo by editing the `docs/` area.
- Be accompanied by a complete Pull Request template (loaded automatically when a PR is created).

For changes that address core functionality or would require breaking changes (e.g. a major release), it's best to open an Issue to discuss your proposal first. This is not required but can save time creating and reviewing changes.

In general, it is recommended to follow the ["fork-and-pull" Git workflow](https://github.com/susam/gitpr)

1. Fork the repository to your own Github account
2. Clone the project to your machine
3. Create a branch locally with a succinct but descriptive name
4. Commit changes to the branch
5. Following any formatting and testing guidelines specific to this repo
6. Push changes to your fork
7. Open a PR in our repository and follow the PR template so that we can efficiently review the changes.

## Contribution Requirements

Please make sure your contribution meets these requirements before submitting:

- Follow the code style guidelines:
  - PEP8 with a max line length of 120 characters.
  - Use absolute imports and group stdlib, third-party, and local imports.
  - Add type hints for public APIs where possible.
  - Write NumPy-style docstrings for functions/classes, following the
    [numpydoc documentation style][numpydoc-style].
  - The [Google Python style guide][google-pyguide] is recommended but not
    enforced.
  - For JavaScript in `wasm/src/**`, add [JSDoc][jsdoc] for exported/public
    functions.
  - For Fortran, preserve scientific comments; use `.f` for fixed-format
    and `.f90` for free-format.
- Run lint and tests as appropriate for your changes:
  - `flake8 . --max-line-length=120`
  - `pytest tests/ test/ --verbose`
  - If Fortran changes or tests require it: `make libphsh` (or
    `python setup.py build_ext --inplace`)
- Update documentation for user-visible changes (Sphinx sources live in
  `docs/`).

## Quality and CI

- Continuous integration runs on pull requests and pushes via
  [GitHub Actions][ci-actions].
- Static analysis includes [CodeQL][codeql-workflow] and other checks configured
  in CI; address findings before release.
- Linting uses flake8, black, and markdownlint; fix linter warnings before
  merging.
- When adding new functionality or fixing bugs, add or update automated tests.

## Getting Help

Please reach out to the [Q&A boards][qa-boards] and post your question there in
the correct category with a descriptive tag.

[code-of-conduct]:
  https://github.com/Liam-Deacon/phaseshifts/blob/master/CODE_OF_CONDUCT.md
[ci-actions]:
  https://github.com/Liam-Deacon/phaseshifts/actions
[codeql-workflow]:
  https://github.com/Liam-Deacon/phaseshifts/actions/workflows/codeql.yml
[google-pyguide]:
  https://google.github.io/styleguide/pyguide.html
[jsdoc]:
  https://jsdoc.app/
[numpydoc-style]:
  https://numpydoc.readthedocs.io/en/latest/format.html
[qa-boards]:
  https://github.com/Liam-Deacon/phaseshifts/discussions/categories/q-a
