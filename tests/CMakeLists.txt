file(GLOB SOURCE_LIST CONFIGURE_DEPENDS
        "${Saturation_SOURCE_DIR}/tests/*.h"
        "${Saturation_SOURCE_DIR}/tests/*.hpp"
        "${Saturation_SOURCE_DIR}/tests/*.cc"
        "${Saturation_SOURCE_DIR}/tests/*.cpp")

# This target runs all of the tests.
# Make sure CATCH_CONFIG_MAIN is only defined in one test file.
add_executable(test-saturation ${SOURCE_LIST})

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(test-saturation PRIVATE saturation catch2)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME test-saturation COMMAND test-saturation) # Command can be a target

# Add folders
set_target_properties(test-saturation PROPERTIES FOLDER iclue)

add_custom_command(
        TARGET test-saturation
        PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -DSRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data" -DDEST_DIR="${CMAKE_CURRENT_BINARY_DIR}/data" -P "${CMAKE_SOURCE_DIR}/cmake/configure_files.cmake"
)

