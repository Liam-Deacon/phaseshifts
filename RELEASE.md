# Release and Versioning

## Versioning

This project follows Semantic Versioning (SemVer) for release versions. Each
release has a unique version identifier, and pre-releases may use suffixes
(e.g., `-rc`, `-dev`).

## Release tags

Releases are tagged in git using `vX.Y.Z` (for example, `v0.1.9`). Tags and
releases are published on GitHub:
https://github.com/Liam-Deacon/phaseshifts/releases

## Release notes

Release notes are maintained in ChangeLog and in the GitHub release entry:
https://github.com/Liam-Deacon/phaseshifts/blob/master/ChangeLog

If a release fixes a publicly known vulnerability (e.g., a CVE), the ChangeLog
and release notes will call it out explicitly.

## Distribution and integrity

Official releases are distributed via GitHub Releases and PyPI over HTTPS. We
avoid using insecure transport for release artifacts or checksums.
