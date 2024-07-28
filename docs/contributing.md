(contributing)=

```{include} ../CONTRIBUTING.md

```

In order to avoid gh actions linting errors, git commit message type must be one of [build, chore, ci, docs, feat, fix, perf, refactor, revert, style, test]

Please make sure to run pre-commit prior to committing to the main repo like so.

```
$ poetry run pre-commit run --all-files
```
