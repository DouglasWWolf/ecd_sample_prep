cmake_minimum_required(VERSION 2.8.9)
project(esp_sample_prep)

# This is the name of the final executable
set(exe esp)

# Specify where all of the header files are
include_directories(src)

# Find the names of all the source files
file(GLOB SOURCES "src/*.cpp")

# This is the name of the executable to build
add_executable(${exe} ${SOURCES})

# After the build, strip debug symbols from the target
add_custom_command(
  TARGET ${exe} POST_BUILD
  COMMAND strip ${exe}
  VERBATIM
)
