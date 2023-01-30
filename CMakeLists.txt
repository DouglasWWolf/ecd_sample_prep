cmake_minimum_required(VERSION 2.8.9)
project(esp_sample_prep)

# This is the name of the final executable
set(EXE esp)

# Specify where all of the header files are
include_directories(src)

# Find the names of all the source files
file(GLOB SOURCES "src/*.cpp")

# Specify what source files our executable is built from
add_executable(${EXE} ${SOURCES})

# After the build, strip debug symbols from the target
add_custom_command(
  TARGET ${EXE} POST_BUILD
  COMMAND strip ${EXE}
  VERBATIM
)
