cmake_minimum_required(VERSION 3.4)
include("${CMAKE_SOURCE_DIR}/FindGfxLibs.cmake")

# Name of the project
project(FinalProj)

# Enforce use of C++14
set(CMAKE_CXX_STANDARD_REQUIRED 14)
set(CMAKE_CXX_STANDARD 14)

# Use glob to get the list of all source files.
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/ext/*/*.cpp" "${CMAKE_SOURCE_DIR}/ext/glad/src/*.c")

# We don't really need to include header and resource files to build, but it's
# nice to have them show up in IDEs.
file(GLOB_RECURSE HEADERS "src/*.h" "ext/*/*.h" "ext/glad/*/*.h")
file(GLOB_RECURSE GLSL "resources/*.glsl")

include_directories("ext")
include_directories("ext/glad/include")
include_directories("ext/assimp/include")  # Path to Assimp include directory

# Set the executable.
add_executable(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS} ${GLSL})

# Link with the pre-built Assimp library
if(WIN32)
  target_link_libraries(${CMAKE_PROJECT_NAME} "${CMAKE_SOURCE_DIR}/ext/assimp/lib/assimp-vc143-mt.lib")
  # If using Debug configuration, you might need this instead:
  # target_link_libraries(${CMAKE_PROJECT_NAME} "${CMAKE_SOURCE_DIR}/ext/assimp/lib/assimp-vc143-mtd.lib")
else()
  target_link_libraries(${CMAKE_PROJECT_NAME} "${CMAKE_SOURCE_DIR}/ext/assimp/lib/libassimp.a")
endif()

# Helper function included from FindGfxLibs.cmake
findGLFW3(${CMAKE_PROJECT_NAME})
findGLM(${CMAKE_PROJECT_NAME})

# OS specific options and libraries
if(NOT WIN32)
  message(STATUS "Adding GCC style compiler flags")
  add_compile_options("-Wall")
  add_compile_options("-pedantic")
  add_compile_options("-Werror=return-type")
 
  # TODO: The following links may be uneeded. 
  if(APPLE)
    # Add required frameworks for GLFW.
    target_link_libraries(${CMAKE_PROJECT_NAME} "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
  else()
    #Link the Linux OpenGL library
    target_link_libraries(${CMAKE_PROJECT_NAME} "GL" "dl")
  endif()
else()
  # Link OpenGL on Windows
  target_link_libraries(${CMAKE_PROJECT_NAME} opengl32.lib)
endif()

# Copy required DLLs to the binary directory
if(WIN32)
  add_custom_command(TARGET ${CMAKE_PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${CMAKE_SOURCE_DIR}/ext/assimp/bin/assimp-vc143-mt.dll"
      $<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>)
endif()