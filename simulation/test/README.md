# How to manage tests

## Generate reference files

Command from here: https://pypi.org/project/pytest-arraydiff/

```bash
py.test --arraydiff-generate-path=reference
```

## Run tests

From the `test/` directory run:

```bash
py.test --cov=simulation --cov-report html --cov-report term
```

the html report is in `htmlcov/index.html`