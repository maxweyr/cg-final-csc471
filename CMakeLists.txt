cmake_minimum_required(VERSION 3.10...3.23)
include("${CMAKE_SOURCE_DIR}/FindGfxLibs.cmake")

# Name of the project
project(10Percent)

# Enforce use of C++17
set(CMAKE_CXX_STANDARD_REQUIRED 17)
set(CMAKE_CXX_STANDARD 17)

# Set OpenGL version requirements
# This will ensure we use OpenGL 4.1 specifically
add_definitions(-DGL_SILENCE_DEPRECATION)
add_definitions(-DGLAD_GLAPI_EXPORT)
add_definitions(-DGLAD_GL_IMPLEMENTATION)
add_definitions(-DGLAD_OPENGL_CORE_VERSION_4_1)

# Use glob to get the list of all source files.
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/ext/*/*.cpp" "${CMAKE_SOURCE_DIR}/ext/glad/src/*.c")

# We don't really need to include header and resource files to build, but it's
# nice to have them show up in IDEs.
file(GLOB_RECURSE HEADERS "src/*.h")
file(GLOB_RECURSE GLSL "resources/*.glsl")

include_directories("ext")
include_directories("ext/glad/include")

# Set up Assimp external dependency
include(ExternalProject)

# Set Assimp version to use
set(ASSIMP_VERSION "v5.4.3")

# Set ext directory for storing dependencies
set(ASSIMP_INSTALL_DIR "${CMAKE_SOURCE_DIR}/ext/assimp")

# Configure Assimp build options
set(ASSIMP_CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${ASSIMP_INSTALL_DIR}
    -DASSIMP_BUILD_TESTS=OFF
    -DASSIMP_BUILD_ASSIMP_TOOLS=OFF
    -DBUILD_SHARED_LIBS=ON
    -DASSIMP_NO_EXPORT=OFF
)

# Check if Assimp is already installed at our desired location
if(EXISTS "${ASSIMP_INSTALL_DIR}/include/assimp/Importer.hpp")
    message(STATUS "Found previously built Assimp at ${ASSIMP_INSTALL_DIR}")
    set(ASSIMP_ALREADY_BUILT TRUE)
else()
    set(ASSIMP_ALREADY_BUILT FALSE)
    
    # Add Assimp as an external project
    ExternalProject_Add(
        assimp_external
        GIT_REPOSITORY https://github.com/assimp/assimp.git
        GIT_TAG ${ASSIMP_VERSION}
        CMAKE_ARGS ${ASSIMP_CMAKE_ARGS}
        PREFIX "${CMAKE_BINARY_DIR}/assimp-build"
        UPDATE_COMMAND ""
    )
endif()

# Set include and library directories for Assimp
set(ASSIMP_INCLUDE_DIRS "${ASSIMP_INSTALL_DIR}/include")

# Set the appropriate library name based on platform
if(WIN32)
    set(ASSIMP_LIBRARIES "${ASSIMP_INSTALL_DIR}/lib/assimp-vc143-mtd.lib")
    set(ASSIMP_DLL "${ASSIMP_INSTALL_DIR}/bin/assimp-vc143-mtd.dll")
elseif(APPLE)
    set(ASSIMP_LIBRARIES "${ASSIMP_INSTALL_DIR}/lib/libassimp.dylib")
else()
    set(ASSIMP_LIBRARIES "${ASSIMP_INSTALL_DIR}/lib/libassimp.so")
endif()

# Create a directory for includes (needed before build)
file(MAKE_DIRECTORY ${ASSIMP_INCLUDE_DIRS})

# Add include directories for Assimp
include_directories(${ASSIMP_INCLUDE_DIRS})

# Set the executable.
add_executable(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS} ${GLSL})

# If we need to build Assimp, add dependency
if(NOT ASSIMP_ALREADY_BUILT)
    add_dependencies(${CMAKE_PROJECT_NAME} assimp_external)
endif()

# Link with Assimp library
target_link_libraries(${CMAKE_PROJECT_NAME} ${ASSIMP_LIBRARIES})

# Helper function included from FindGfxLibs.cmake
findGLFW3(${CMAKE_PROJECT_NAME})
findGLM(${CMAKE_PROJECT_NAME})

# OS specific options and libraries
if(NOT WIN32)
  message(STATUS "Adding GCC style compiler flags")
  add_compile_options("-Wall")
  add_compile_options("-pedantic")
  add_compile_options("-Werror=return-type")
 
  if(APPLE)
    # Add required frameworks for GLFW.
    target_link_libraries(${CMAKE_PROJECT_NAME} "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
  else()
    # Link the Linux OpenGL library
    target_link_libraries(${CMAKE_PROJECT_NAME} "GL" "dl")
  endif()
else()
  # Link OpenGL on Windows
  target_link_libraries(${CMAKE_PROJECT_NAME} opengl32.lib)
  
  # Copy DLL to output directory
  add_custom_command(TARGET ${CMAKE_PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${ASSIMP_DLL}"
      $<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>)
endif()

# Add a target to clean the Assimp installation if needed
add_custom_target(clean_assimp
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${ASSIMP_INSTALL_DIR}
    COMMENT "Cleaning Assimp installation directory"
)