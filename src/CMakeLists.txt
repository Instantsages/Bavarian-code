cmake_minimum_required(VERSION 3.8.0) # required because of compiler meta-features such as cxx_std_11
project(BavarianDrivers VERSION 0.01 LANGUAGES CXX)

# Set the full lib install path as rpath
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
	set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# Find Intel TBB library, if requested.
if (BAVARIAN_WITH_EXECUTION)
	if(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
		message(FATAL_ERROR
			"BAVARIAN_WITH_EXECUTION is set to On, but the clang compiler does not
			support C++17 execution yet (as of version 10). Aborting.")
	elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 9.1)
		message(FATAL_ERROR
			"BAVARIAN_WITH_EXECUTION is set to On, but g++ < 9.1 does not support
			C++17 execution. Aborting.")
	endif()
		list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/FindTBB)
		find_package(TBB)
		if (NOT TBB_FOUND)
			message(FATAL_ERROR
				"BAVARIAN_WITH_EXECUTION is set to On, but the Intel TBB library
				cannot be found. Aborting.")
		endif()
endif()

set(NETWORKIT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/networkit)
if(EXISTS "${NETWORKIT_SOURCE_DIR}")
	set(NETWORKIT_NATIVE ON CACHE BOOL "Use -march=native in NetworKit")
	set(NETWORKIT_CXX_STANDARD 17 CACHE STRING "C++ standard in NetworKit")
	#set(NETWORKIT_RELEASE_LOGGING ON CACHE BOOL) # To avoid log messages at levels TRACE or DEBUG
	set(NETWORKIT_QUIET_LOGGING ON CACHE BOOL
		"Set default log level to QUIET in NetworKit (can still be changed at runtime)")
	add_subdirectory(${NETWORKIT_SOURCE_DIR})
else()
	message(FATAL_ERROR
		"Missing networkit sources in src/networkit. "
		"Run `git submodule update --init --recursive` to get the submodule and"
		" its submodules.")
endif()

set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/libomp/18.1.3" ${CMAKE_PREFIX_PATH})
find_package(OpenMP REQUIRED)

add_executable(exact exact.cpp)
target_compile_features(exact PRIVATE cxx_std_17)
target_compile_options(exact PRIVATE ${WARNS_COMPILE_FLAGS})
target_compile_options(exact PRIVATE ${NONWARNS_COMPILE_FLAGS})
target_include_directories(exact PRIVATE ${NETWORKIT_SRC_DIR})
# target_compile_options(exact PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(exact PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(exact PRIVATE networkit)
target_link_libraries(exact PRIVATE tlx)
install(TARGETS exact RUNTIME DESTINATION bin)

add_executable(bavarian bavarian.cpp)
target_compile_features(bavarian PRIVATE cxx_std_17)
target_compile_options(bavarian PRIVATE ${WARNS_COMPILE_FLAGS})
target_compile_options(bavarian PRIVATE ${NONWARNS_COMPILE_FLAGS})
target_include_directories(bavarian PRIVATE ${NETWORKIT_SRC_DIR})
# target_compile_options(bavarian PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(bavarian PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(bavarian PRIVATE networkit)
target_link_libraries(bavarian PRIVATE tlx)
if (BAVARIAN_WITH_EXECUTION)
	target_compile_definitions(bavarian PRIVATE BAVARIAN_WITH_EXECUTION)
	target_link_libraries(bavarian PRIVATE ${TBB_LIBRARIES})
endif()
install(TARGETS bavarian RUNTIME DESTINATION bin)

add_executable(progrbavarian progrbavarian.cpp)
target_compile_features(progrbavarian PRIVATE cxx_std_17)
target_compile_options(progrbavarian PRIVATE ${WARNS_COMPILE_FLAGS})
target_compile_options(progrbavarian PRIVATE ${NONWARNS_COMPILE_FLAGS})
target_include_directories(progrbavarian PRIVATE ${NETWORKIT_SRC_DIR})
# target_compile_options(progrbavarian PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(progrbavarian PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(progrbavarian PRIVATE networkit)
target_link_libraries(progrbavarian PRIVATE tlx)
if (BAVARIAN_WITH_EXECUTION)
	target_compile_definitions(progrbavarian PRIVATE BAVARIAN_WITH_EXECUTION)
	target_link_libraries(progrbavarian PRIVATE ${TBB_LIBRARIES})
endif()
install(TARGETS progrbavarian RUNTIME DESTINATION bin)


add_executable(schedule schedule.cpp)
target_compile_features(schedule PRIVATE cxx_std_17)
target_compile_options(schedule PRIVATE ${WARNS_COMPILE_FLAGS})
target_compile_options(schedule PRIVATE ${NONWARNS_COMPILE_FLAGS})
target_include_directories(schedule PRIVATE ${NETWORKIT_SRC_DIR})
# target_compile_options(schedule PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(schedule PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(schedule PRIVATE networkit)
target_link_libraries(schedule PRIVATE tlx)
install(TARGETS schedule RUNTIME DESTINATION bin)
