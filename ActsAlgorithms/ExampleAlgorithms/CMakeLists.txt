add_library(ExampleAlgorithms SHARED)

# ----------------------------------------
# Sources
# ----------------------------------------
target_sources(ExampleAlgorithms
    PRIVATE
    src/ExampleAlgorithm.cpp
    src/TrackTruthMatcher.cpp
    src/HitRemoverAlgorithm.cpp
    src/TrackMergerAlgorithm.cpp
        # add more .cpp files here
)

# ----------------------------------------
# Include directories
# ----------------------------------------
target_include_directories(ExampleAlgorithms
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)


find_library(ACTS_EXAMPLES_FRAMEWORK
    NAMES ActsExamplesFramework
    HINTS ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES lib lib64
    REQUIRED
)

message(STATUS "ACTS_EXAMPLES_FRAMEWORK: ${ACTS_EXAMPLES_FRAMEWORK}")


# ----------------------------------------
# Link against Acts
# ----------------------------------------
target_link_libraries(ExampleAlgorithms
    PUBLIC
    Acts::Core
    ${ACTS_EXAMPLES_FRAMEWORK}
  )
