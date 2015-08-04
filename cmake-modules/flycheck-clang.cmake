cmake_minimum_required(VERSION 2.8)


#get a property list and convert it into a space separated
#list of quoted strings
function(get_property_string_list dir property outvarname)

  get_property(paths DIRECTORY ${dir} PROPERTY ${property})

  foreach(path ${paths})
    if(path_list)
      set(path_list "${path_list} \"${path}\"")
    else()
      set(path_list "\"${path}\"")
    endif()
  endforeach()

  set(${outvarname} ${path_list} PARENT_SCOPE)
endfunction()

# Adds -I and -D to include directories and compiler defines
function(get_raw_flags includes defines outvarname)
  string(REGEX REPLACE "\"([^ ])" "\"-I\\1" raw_includes ${includes})
  string(REGEX REPLACE " " "\n      " raw_includes ${raw_includes})

  string(REGEX REPLACE "\"([^ ])" "\"-D\\1" raw_defines ${defines})
  string(REGEX REPLACE " " "\n      " raw_defines ${raw_defines})

  set(${outvarname} "${raw_includes}\n      ${raw_defines}" PARENT_SCOPE)
endfunction()


# Get include directories and compiler definitions from CMake itself
# and write an Emacs file to pass these to clang when it executes
# in the background
function(write_dir_locals_el dir)
  get_property_string_list(${dir} INCLUDE_DIRECTORIES INCLUDE_PATHS)
  string(REGEX REPLACE " " "\n                                 " PRETTY_PATHS ${INCLUDE_PATHS})

  get_property_string_list(${dir} COMPILE_DEFINITIONS DEFINITIONS)
  string(REGEX REPLACE " " "\n                                 " PRETTY_DEFS ${DEFINITIONS})

  get_raw_flags(${INCLUDE_PATHS} ${DEFINITIONS} AC_CLANG_FLAGS)

  set(EMACS_DIR_LOCALS ${dir}/.dir-locals.el)

  file(WRITE ${EMACS_DIR_LOCALS} ";;; Directory Local Variables
;;; See Info node `(emacs) Directory Variables' for more information.

(
  (nil
    (flycheck-clang-include-path ${PRETTY_PATHS})
    (flycheck-clang-definitions ${PRETTY_DEFS}))
  (c-mode
    (eval setq ac-clang-flags (append ac-clang-flags-c '(
      ${AC_CLANG_FLAGS}
    ))))
  (c++-mode
    (eval setq ac-clang-flags (append ac-clang-flags-c++ '(
      ${AC_CLANG_FLAGS}
    ))))
)
")
endfunction()