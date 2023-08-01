# set minimum cmake version
cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

# project name and language
project(mixdyn
        VERSION 0.1
        DESCRIPTION "Implicit-Explicit transient dynamic analysis"
        LANGUAGES Fortran)

add_subdirectory(src)
               
