cmake_minimum_required (VERSION 3.21)
project (PROJECT)
enable_testing()

add_subdirectory(src)

add_executable(readmap src/readmap.c src/fasta.c src/fastq.c src/sam.c)
target_link_libraries(readmap cstr)
