# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.
# Optionally glob, but only for CMake 3.12 or later:
file(GLOB HEADER_LIST CONFIGURE_DEPENDS
        "${Saturation_SOURCE_DIR}/include/saturation/*.h"
        "${Saturation_SOURCE_DIR}/include/saturation/*.hpp")

file(GLOB SOURCE_LIST CONFIGURE_DEPENDS
        "${Saturation_SOURCE_DIR}/src/*.h"
        "${Saturation_SOURCE_DIR}/src/*.hpp"
        "${Saturation_SOURCE_DIR}/src/*.cc"
        "${Saturation_SOURCE_DIR}/src/*.cpp")

add_library(saturation ${SOURCE_LIST} ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(saturation PUBLIC "${Saturation_SOURCE_DIR}/include/")

target_link_libraries(saturation PUBLIC nlnum-lib PRIVATE ${Python3_LIBRARIES} cppitertools prettyprint OpenMP::OpenMP_CXX)

# IDEs should put the headers in a nice place
source_group(TREE "${PROJECT_SOURCE_DIR}/include" PREFIX "Header Files" FILES ${HEADER_LIST})
