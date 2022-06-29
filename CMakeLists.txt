cmake_minimum_required(VERSION 3.0)

project(FastVC)
option(STATIC "build dependency-free binary" OFF)

#openmp
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_CXX_COMPILER "g++")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#link
#HTS_PREFIX
if(DEFINED HTS_PREFIX)
else()
	set(HTS_PREFIX "/usr/local/htslib")	
endif()
message("hts dir: " ${HTS_PREFIX})

LINK_DIRECTORIES(${HTS_PREFIX}/lib)
#LINK_DIRECTORIES(/home/user_home/haoz/tools/lib)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -g -ffast-math ${CMAKE_CXX_FLAGS}")

if(STATIC)
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-intel -no-intel-extensions -static-libstdc++ -qopenmp-link=static -lz")
endif ()

#add_subdirectory(src)
aux_source_directory(./src SOURCE_LIST)
aux_source_directory(./src/modes SOURCE_LIST)
message(${SOURCE_LIST})
include_directories(./include ./include/modes ${HTS_PREFIX}/include)

add_executable (FastVC ${SOURCE_LIST})

if(STATIC)
	target_link_libraries(FastVC libhts.a)
else()
	target_link_libraries(FastVC hts)
	#target_link_libraries(FastVC libtcmalloc.a)
endif()

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
install(TARGETS FastVC DESTINATION bin)

