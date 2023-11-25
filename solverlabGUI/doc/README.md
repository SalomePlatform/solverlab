# How to build Sphinx documentation (Linux)

The user and developer documentation is based on Sphinx. The build process is
automated with a Makefile. The file `unused_CMakeLists.txt` is broken, it must
not be used to build the documentation. 

## Dependencies

- python
- sphinx-build (`python3 -m pip install sphinx-build`) 

## Build

```sh
make html
```

## Open

```sh
firefox build/html/index.html &
```
