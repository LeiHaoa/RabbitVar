cmake_minimum_required(VERSION 3.0)

project(FastVC)
#openmp
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_CXX_COMPILER "icpc")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(EXECUTABLE_OUTPUT_PATH ../)
endif()
#target_link_libraries(FastVC
#  PUBLIC
#    OpenMP::OpenMP_CXX
#  )

#find_package(PkgConfig REQUIRED)
#pkg_check_modules(HTSlib REQUIRED htslib)
#set(HTSlib_LIBRARIES "${HTSlib_LDFLAGS}")
#find_package(htslib REQUIRED)
#if(htslib_find)
#message("=============")
#endif()
message(HTSlib_LIBRARIES)

#link
LINK_DIRECTORIES(/home/old_home/haoz/tools/htslib/lib)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -g -ffast-math ${CMAKE_CXX_FLAGS}")

#add_subdirectory(src)
aux_source_directory(./src SOURCE_LIST)
aux_source_directory(./src/modes SOURCE_LIST)
message(${SOURCE_LIST})
include_directories(./include ./include/modes /home/old_home/haoz/tools/htslib/include)

add_executable (FastVC ${SOURCE_LIST})

target_link_libraries(FastVC hts)

