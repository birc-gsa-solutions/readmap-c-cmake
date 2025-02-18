project(cstr
    LANGUAGES   C
    DESCRIPTION "String algorithms in C."
)
add_library(${PROJECT_NAME} STATIC)
add_library(${PROJECT_NAME}::framework ALIAS ${PROJECT_NAME})

configure_file(config.h.in config.h)

file(GLOB SOURCES ./*.h ./*.c)
target_sources(${PROJECT_NAME}
    PRIVATE 
        ${SOURCES}
)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set(${PROJECT_NAME}_PUBLIC_HEADERS
        cstr.h
)

set_target_properties(${PROJECT_NAME}
    PROPERTIES 
        C_STANDARD          23
        C_STANDARD_REQUIRED ON
        C_EXTENSIONS        OFF
        PUBLIC_HEADER       "${${PROJECT_NAME}_PUBLIC_HEADERS}"
        EXPORT_NAME         framework
)

install(TARGETS   ${PROJECT_NAME}
    EXPORT        ${PROJECT_NAME}Targets
    ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT     library
)

target_compile_options(${PROJECT_NAME}
    PRIVATE
        -Wcast-align
        -Wcast-qual
        -Wconversion
        -Wmissing-noreturn
        -Wmissing-prototypes
        -Wno-missing-braces
        -Wshadow
        -Werror
        -Wall
        -O3 # optimisation necessary for CPS trick in search
)

# If we are building test code, then we need to include testlib
# so we can create the unit test functions.
target_link_libraries(${PROJECT_NAME} testlib)
