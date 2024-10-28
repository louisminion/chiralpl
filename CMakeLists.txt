cmake_minimum_required(VERSION 3.10)
project("chiralpl" LANGUAGES "Fortran")
set (FC ifx)
add_executable("chiralpl" "chiralpl.f90")