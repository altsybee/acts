# Add the examples python binding
pybind11_add_module(AliceActsPythonBindings
  src/AliceActsAlgorithmsModule.cpp
)

target_compile_features(AliceActsPythonBindings PUBLIC cxx_std_20)


target_include_directories(AliceActsPythonBindings
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../ExampleAlgorithms/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../Validation/include
    )

#link_directories(${CMAKE_PREFIX_PATH}/lib64)
#message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

find_library(ACTS_EXAMPLES_FRAMEWORK
    NAMES ActsExamplesFramework
    HINTS ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES lib lib64
    REQUIRED
  )

find_library(ACTS_EXAMPLES_IO_ROOT
    NAMES ActsExamplesIoRoot
    HINTS ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES lib lib64
    REQUIRED
  )

target_link_libraries(AliceActsPythonBindings
  PUBLIC
  ExampleAlgorithms
  #Validation
  ${ACTS_EXAMPLES_FRAMEWORK}
  ${ACTS_EXAMPLES_IO_ROOT}
  Acts::Core
)
