# set minimum cmake version
cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

# project
project(mixdyn
    VERSION 0.1
    DESCRIPTION "Implicit-Explicit transient dynamic analysis"
    HOMEPAGE_URL https://github.com/PaoloBertin/mixdyn
    LANGUAGES Fortran
    )

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f95 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()

add_subdirectory(src)
               
