cmake_minimum_required(VERSION 3.8.0) # required because of compiler meta-features such as cxx_std_11
project(Bavarian VERSION 0.01 LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # must come after project()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
	set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)
endif()

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release
		CACHE STRING "Choose the type of build : None Debug Release RelWithDebugLogging."
		FORCE)
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
	if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 4.0)
		set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -Og")
	else()
		set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -O0")
	endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
	set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3")
else()
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL FreeBSD AND ${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
	if (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 9.1 AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 10.0 )
		set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath=/usr/local/lib/gcc9")
	elseif (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 10.0)
		set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=/usr/local/lib/gcc10")
		set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath=/usr/local/lib/gcc10")
	else ()
		message("You may want to set CMAKE_SHARED_LINKER_FLAGS to '-Wl,-rpath=/usr/local/lib/gccX', where 'X' is your g++ major version number.")
	endif()
endif()

set(WARNS_COMPILE_FLAGS -Wextra -Wall -Werror -pedantic-errors -Wformat=2
	-Wpointer-arith -Wcast-qual -Wno-unused-result -Wno-deprecated-declarations) #-Winline)
set(NONWARNS_COMPILE_FLAGS "-march=native")

option(BAVARIAN_WITH_EXECUTION "Build Bavarian with the C++17 execution policies. Relies on the Intel ThreadsBuildingBlocks (tbb) library" OFF)

add_subdirectory(src)
