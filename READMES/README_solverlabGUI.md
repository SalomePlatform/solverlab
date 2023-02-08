
![logo](../images/logoSOLVERLAB.png)

## SOLVERLAB_GUI

### Introduction

- Is a new Graphical User Interface (GUI/IHM) developped for SOLVERLAB,
  this GUI uses the python-library PACKAGESPY.

- To use the GUI, set the cmake option compilation `-DSOLVERLAB_WITH_PACKAGESPY=ON`.

- if you have a local copy of the library
  [PACKAGESPY](https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git)
  (as in SALOME tuleap contex), give the link to solverlab using the cmake option
  `-DPACKAGESPY_ROOT_DIR=${PACKAGESPY_ROOT_DIR}`.

- Finally, to launch the Graphical User Interface of SOLVERLAB,
  you have to load the SOLVERLAB environment *in your terminal*, type:

  ```bash
  source .../SOLVERLAB_install/env_SOLVERLAB.sh
  ${SOLVERLAB_ROOT_DIR}/solverlabGUI --help # help CLI
  ${SOLVERLAB_ROOT_DIR}/solverlabGUI --gui  # lanch GUI
  ```


### PACKAGESPY (Tuleap) *outside* SALOME compilation/prerequisites

- PACKAGESPY is the *main* mandatory prerequisite for SOLVERLAB_GUI.

- The *reference* (SALOME-MATIX) base git for PACKAGESPY is
  `https://codev-tuleap.intra.cea.fr/plugins/git/matix/packagespy.git`

- NOTE: This web site is not reacheable for everybody, unhappily.

