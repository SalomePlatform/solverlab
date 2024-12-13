# Test Salome - SOLVERLAB

These tests should check the use of SolverLab within Salome environment.

## Convert jupyter notebooks into python files

Include the notebooks in the tests to see if they are up to date.

```bash
cd ../notebooks
jupyter-nbconvert *.ipynb --to script
mv *.py ../test/notebooks
```
