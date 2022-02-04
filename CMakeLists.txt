cmake_minimum_required (VERSION 3.21)
project (PROJECT)
enable_testing()

add_executable(example src/example.c)
add_test(example example)

