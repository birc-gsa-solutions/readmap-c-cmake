cmake_minimum_required(VERSION 3.21)
project(readmap)
enable_testing()
set(CMAKE_C_FLAGS "-O3 -Wall")

add_subdirectory(src)

add_executable(readmap src/readmap.c src/fasta.c src/fastq.c src/sam.c)
target_link_libraries(readmap cstr)
