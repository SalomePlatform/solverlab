#[==[.md
# paraview-config.cmake

This file is used by CMake when finding ParaView.

The following variables are provided by this module:

  * `ParaView_VERSION`: The version of ParaView found.
  * `ParaView_PREFIX_PATH`: Install prefix for ParaView.
  * `PARAVIEW_USE_QT`: If ParaView's Qt GUI is available.
  * `PARAVIEW_USE_MPI`: If ParaView is built with MPI support.
  * `PARAVIEW_USE_PYTHON`: If ParaView is built with Python support.
  * `PARAVIEW_PYTHONPATH`: Where ParaView's Python modules live under the
    install prefix. Unset if Python is not available.
  * `PARAVIEW_PLUGIN_SUBDIR`: The subdirectory under the library directory for
    plugins.
  * `ParaView_CLIENT_XML_FILES`: XML files for client applications to use to
    reproduce ParaView's menu items. Only provided if `PARAVIEW_USE_QT`
    is set.
  * `ParaView_LIBRARIES`: The list of modules specified by `COMPONENTS` and
    `OPTIONAL_COMPONENTS`. This may be used in `MODULES` arguments in the API
    (e.g., `vtk_module_autoinit`). All modules are also targets and may be
    linked to using `target_link_libraries`.
#]==]

set(${CMAKE_FIND_PACKAGE_NAME}_CMAKE_MODULE_PATH_save "${CMAKE_MODULE_PATH}")
list(INSERT CMAKE_MODULE_PATH 0
  "${CMAKE_CURRENT_LIST_DIR}")

set("${CMAKE_FIND_PACKAGE_NAME}_CMAKE_PREFIX_PATH_save" "${CMAKE_PREFIX_PATH}")
include("${CMAKE_CURRENT_LIST_DIR}/paraview-prefix.cmake")
set("${CMAKE_FIND_PACKAGE_NAME}_PREFIX_PATH"
  "${_vtk_module_import_prefix}")
unset(_vtk_module_import_prefix)
list(INSERT CMAKE_PREFIX_PATH 0
  "${${CMAKE_FIND_PACKAGE_NAME}_PREFIX_PATH}")

set("${CMAKE_FIND_PACKAGE_NAME}_VERSION" "@PARAVIEW_VERSION_FULL@")

unset("${CMAKE_FIND_PACKAGE_NAME}_FOUND")

set(_paraview_use_external_vtk "@PARAVIEW_USE_EXTERNAL_VTK@")
set(_paraview_find_package_args)
if (NOT _paraview_use_external_vtk)
  list(APPEND _paraview_find_package_args
    PATHS "${CMAKE_CURRENT_LIST_DIR}/vtk"
    NO_DEFAULT_PATH)
endif ()
if (${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
  list(APPEND _paraview_find_package_args
    QUIET)
endif ()
find_package(VTK REQUIRED
  ${_paraview_find_package_args})
if (NOT VTK_FOUND)
  set("${CMAKE_FIND_PACKAGE_NAME}_FOUND" 0)
endif ()
unset(_paraview_find_package_args)
unset(_paraview_use_external_vtk)

set(PARAVIEW_USE_QT "@PARAVIEW_USE_QT@")
set(PARAVIEW_USE_MPI "@PARAVIEW_USE_MPI@")
set(PARAVIEW_USE_PYTHON "@PARAVIEW_USE_PYTHON@")
set(PARAVIEW_PLUGIN_SUBDIR "paraview@paraview_version_suffix@/plugins")

if (PARAVIEW_USE_PYTHON)
  set(PARAVIEW_PYTHONPATH "@PARAVIEW_PYTHON_SITE_PACKAGES_SUFFIX@")
  include("${CMAKE_CURRENT_LIST_DIR}/ParaViewPython-targets.cmake")
  # Unset this for now; these targets will be defined later.
  unset("${CMAKE_FIND_PACKAGE_NAME}_FOUND")
  unset("${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE")
endif ()

include("${CMAKE_CURRENT_LIST_DIR}/${CMAKE_FIND_PACKAGE_NAME}-targets.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/${CMAKE_FIND_PACKAGE_NAME}-vtk-module-properties.cmake")

include("${CMAKE_CURRENT_LIST_DIR}/paraview-find-package-helpers.cmake" OPTIONAL)
include("${CMAKE_CURRENT_LIST_DIR}/${CMAKE_FIND_PACKAGE_NAME}-vtk-module-find-packages.cmake")

include("${CMAKE_CURRENT_LIST_DIR}/${CMAKE_FIND_PACKAGE_NAME}Plugins-paraview_plugins-targets-depends.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/${CMAKE_FIND_PACKAGE_NAME}Plugins-targets.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/${CMAKE_FIND_PACKAGE_NAME}Plugins-paraview-plugin-properties.cmake")

include("${CMAKE_CURRENT_LIST_DIR}/ParaViewTools-targets.cmake" OPTIONAL)

include("${CMAKE_CURRENT_LIST_DIR}/ParaViewClient.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ParaViewPlugin.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ParaViewServerManager.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ParaViewTesting.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/vtkModuleWrapClientServer.cmake")

if (@PARAVIEW_USE_PYTHON@) # PARAVIEW_USE_PYTHON
  include("${CMAKE_CURRENT_LIST_DIR}/paraview.modules-vtk-python-module-properties.cmake")
endif ()

if (@PARAVIEW_USE_QT@) # PARAVIEW_USE_QT
  include("${CMAKE_CURRENT_LIST_DIR}/${CMAKE_FIND_PACKAGE_NAME}Client-targets.cmake")
  include("${CMAKE_CURRENT_LIST_DIR}/ParaView-client-xml.cmake")
endif ()

set(_paraview_components_to_check)
foreach (_paraview_component IN LISTS "${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS")
  if (DEFINED "${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_FOUND")
    # It was already not-found (likely due to `find-package` failures).
  elseif (TARGET "${CMAKE_FIND_PACKAGE_NAME}::${_paraview_component}")
    list(APPEND _paraview_components_to_check
      "${_paraview_component}")
  else ()
    set("${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_FOUND" 0)
    list(APPEND "${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_NOT_FOUND_MESSAGE"
      "The ${_paraview_component} component is not available.")
  endif ()
endforeach ()
unset(_paraview_component)

set(_paraview_vtk_components)

while (_paraview_components_to_check)
  list(GET _paraview_components_to_check 0 _paraview_component)
  list(REMOVE_AT _paraview_components_to_check 0)
  if (DEFINED "${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_FOUND")
    # We've already made a determiniation.
    continue ()
  endif ()

  get_property(_paraview_dependencies
    TARGET    "${CMAKE_FIND_PACKAGE_NAME}::${_paraview_component}"
    PROPERTY  "INTERFACE_paraview_module_depends")
  string(REPLACE "${CMAKE_FIND_PACKAGE_NAME}::" "" _paraview_dependencies "${_paraview_dependencies}")
  set(_paraview_all_dependencies_checked TRUE)
  foreach (_paraview_dependency IN LISTS _paraview_dependencies)
    # Handle VTK module dependencies.
    string(FIND "${_paraview_component}" "VTK::" _paraview_vtk_idx)
    if (NOT _paraview_vtk_idx EQUAL -1)
      unset(_paraview_vtk_idx)
      if (NOT TARGET "${_paraview_dependency}")
        set("${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_FOUND" 0)
        list(APPEND "${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_NOT_FOUND_MESSAGE"
          "Failed to find the ${_paraview_dependency} module.")
      endif ()
      continue ()
    endif ()
    unset(_paraview_vtk_idx)

    if (DEFINED "${CMAKE_FIND_PACKAGE_NAME}_${_paraview_dependency}_FOUND")
      if (NOT ${CMAKE_FIND_PACKAGE_NAME}_${_paraview_dependency}_FOUND)
        set("${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_FOUND" 0)
        list(APPEND "${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_NOT_FOUND_MESSAGE"
          "Failed to find the ${_paraview_dependency} component.")
      endif ()
    else ()
      # Check its dependencies.
      list(APPEND _paraview_components_to_check
        "${_paraview_dependency}")
      set(_paraview_all_found FALSE)
    endif ()
  endforeach ()
  if (NOT DEFINED "${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_FOUND")
    if (_paraview_all_dependencies_checked)
      set("${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_FOUND" 1)
    else ()
      list(APPEND _paraview_components_to_check
        "${_paraview_component}")
    endif ()
  endif ()
  unset(_paraview_all_dependencies_checked)
  unset(_paraview_dependency)
  unset(_paraview_dependencies)
endwhile ()
unset(_paraview_component)
unset(_paraview_components_to_check)

set(_paraview_missing_components)
foreach (_paraview_component IN LISTS "${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS")
  if (NOT ${CMAKE_FIND_PACKAGE_NAME}_${_paraview_component}_FOUND AND ${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_${_paraview_component})
    list(APPEND _paraview_missing_components
      "${_paraview_component}")
  endif ()
endforeach ()

if (_paraview_missing_components)
  list(REMOVE_DUPLICATES _paraview_missing_components)
  list(SORT _paraview_missing_components)
  string(REPLACE ";" ", " _paraview_missing_components "${_paraview_missing_components}")
  set("${CMAKE_FIND_PACKAGE_NAME}_FOUND" 0)
  set("${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE"
    "Could not find the ${CMAKE_FIND_PACKAGE_NAME} package with the following required components: ${_paraview_missing_components}.")
endif ()
unset(_paraview_missing_components)

set("${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES")
if (NOT DEFINED "${CMAKE_FIND_PACKAGE_NAME}_FOUND")
  # If nothing went wrong, we've successfully found the package.
  set("${CMAKE_FIND_PACKAGE_NAME}_FOUND" 1)
  # Build the `_LIBRARIES` variable.
  foreach (_paraview_component IN LISTS "${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS")
    list(APPEND "${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES"
      "${CMAKE_FIND_PACKAGE_NAME}::${_paraview_component}")
  endforeach ()
  unset(_paraview_component)
endif ()

set(CMAKE_PREFIX_PATH "${${CMAKE_FIND_PACKAGE_NAME}_CMAKE_PREFIX_PATH_save}")
unset("${CMAKE_FIND_PACKAGE_NAME}_CMAKE_PREFIX_PATH_save")

set(CMAKE_MODULE_PATH "${${CMAKE_FIND_PACKAGE_NAME}_CMAKE_MODULE_PATH_save}")
unset(${CMAKE_FIND_PACKAGE_NAME}_CMAKE_MODULE_PATH_save)

# Compatibility with old code.
if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
  set(PARAVIEW_USE_FILE "${CMAKE_CURRENT_LIST_DIR}/paraview-use-file-deprecated.cmake")
elseif (${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION VERSION_LESS 5.7)
  set(PARAVIEW_USE_FILE "${CMAKE_CURRENT_LIST_DIR}/paraview-use-file-compat.cmake")
else ()
  set(PARAVIEW_USE_FILE "${CMAKE_CURRENT_LIST_DIR}/paraview-use-file-error.cmake")
endif ()

# 5.8 renamed these variables, so provide them if 5.8 is not the minimum
# requested.
if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION OR
    ${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION VERSION_LESS "5.8")
  set(PARAVIEW_BUILD_QT_GUI "${PARAVIEW_USE_QT}")
  set(PARAVIEW_ENABLE_PYTHON "${PARAVIEW_USE_PYTHON}")
endif ()
