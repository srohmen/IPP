cmake_minimum_required(VERSION 3.5)

message(STATUS "build type: ${CMAKE_BUILD_TYPE}" )

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
	message(STATUS "found ccache and using it")
	set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
include(subdirlist)

set(PROJECT_LIBRARY_NAME "IPP")
project(${PROJECT_LIBRARY_NAME}Project)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

find_package(VTK REQUIRED COMPONENTS
	IOXML vtkRenderingCore vtkRenderingOpenGL2 vtkRenderingImage
	vtkChartsCore  vtkViewsContext2D vtkRenderingContextOpenGL2
	NO_MODULE
)
include(${VTK_USE_FILE})

find_package(Boost REQUIRED COMPONENTS
	mpi
	iostreams
	filesystem
	program_options
	serialization
	chrono
)
find_package(ZLIB REQUIRED)

#add_definitions(-DDEBUG_OUTPUT)

# add_compile_options appends at the end of command line (which is not desired here all the time)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#add_compile_options(-std=c++1y)
add_compile_options(-Wno-long-long)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

set(IPP_COMPILE_PROFILE OFF CACHE BOOL "enable compile time profiling")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	message(STATUS "using Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}
		-Wno-deprecated-declarations
		-Wno-error-unused-variable
		-Wno-error=unused-variable
		-Wno-error-unused-const-variable
		-Wno-error-unused-private-field
		-Wno-error-unused-function
		-Wno-error=unused-function
		-Wno-error=unused-local-typedefs
		-Wno-inconsistent-missing-override"
		# -Wno-unknown-warning-option
	)
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error"  ) # some no error statements above just do not work...
	add_definitions(-DBOOST_PP_VARIADICS=1)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message(STATUS "using GCC")
	#  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
	#          message(FATAL_ERROR "GCC version must be at least 5.0!")
	#  endif()

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -Wno-deprecated-declarations")
	# -Wno-error=unused-variable
	# -Wno-error=unused-but-set-variable
	# -Wno-error=unused-function
	# -Wno-error=unused-local-typedefs"


	# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-implicit-templates")

	# strange error in palabos code:
	# plb::BlockLattice2D<T, Descriptor>::stream(plb::Box2D)
	# cc1plus: Fehler: assuming signed overflow does not occur when assuming that (X - c) <= X is always true [-Werror=strict-overflow]
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-error=strict-overflow")

	if(IPP_COMPILE_PROFILE)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v -ftime-report") # -Q
	endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	message(STATUS "using Intel C++")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -Wno-deprecated-declarations")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	message(STATUS "using Visual Studio C++")

endif()



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPLB_MPI_PARALLEL -DUSE_MPI")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPLB_MPI_PARALLEL")

if(UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPLB_USE_POSIX")
endif()

# write debug output for all debug builds
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DPLB_DEBUG -DDEBUG_OUTPUT -DIPP_DEBUG" )
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DDEBUG_OUTPUT" )


# enable memcheck for all configurations if flag is enabled
set(IPP_MEMCHECK_SANITIZE FALSE CACHE BOOL "enable memory debugging via sanitize")
if(IPP_MEMCHECK_SANITIZE)
	message(STATUS "enabling memory sanitize checks")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined -fsanitize=leak -fno-omit-frame-pointer" )
endif()

set(IPP_ENABLE_EXPLICIT_TEMPLATES_INSTANTIATION FALSE CACHE BOOL "enable explicit templates in order to speed up compilation time")
if(IPP_ENABLE_EXPLICIT_TEMPLATES_INSTANTIATION)
	add_definitions(-DEXPLICIT_INSTANTS)
endif()

set(IPP_ENABLE_LTO FALSE CACHE BOOL "enable link time optimization (LTO)")
if(IPP_ENABLE_LTO)
	message(STATUS "enabling LTO")
	set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

	# for some reason above does not work with current cmake versions for GCC
	# hacking it into flags for GCC
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		add_compile_options("-flto")
		add_link_options("-flto")
	endif()
endif()


set(IPP_PGO_DIR "" CACHE STRING "profile guided optimization (PGO) directory")


set(IPP_ENABLE_PGO_GEN FALSE CACHE BOOL "enable profile guided optimization (PGO) - generation of profiling data")
if(IPP_ENABLE_PGO_GEN)
	message(STATUS "enabling PGO - generation of profiling")

	# hacking it into flags for GCC
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		set(FLAGS "-fprofile-generate=${IPP_PGO_DIR}")
		message(STATUS "PGO flags: ${FLAGS}")
		add_compile_options(${FLAGS})
		add_link_options(${FLAGS})
	endif()
endif()

set(IPP_ENABLE_PGO_USE FALSE CACHE BOOL "enable profile guided optimization (PGO) - usage of profiling data")
if(IPP_ENABLE_PGO_USE)
	message(STATUS "enabling PGO  - usage of profiling ")

	# hacking it into flags for GCC
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		set(FLAGS "-fprofile-use=${IPP_PGO_DIR}")
		message(STATUS "PGO flags: ${FLAGS}")
		add_compile_options(${FLAGS})
		add_link_options(${FLAGS})
	endif()
endif()


set(IPP_ENABLE_RENDERING FALSE CACHE BOOL "enable rendering engine")

set(IPP_ENABLE_BENCHMARK FALSE CACHE BOOL "enable benchmark timing output")
if(IPP_ENABLE_BENCHMARK)
	message(STATUS "enabling additional benchmark timings")
	add_definitions(-DIPP_ENABLE_BENCHMARK)
endif()


set(IPP_ENABLE_TESTS FALSE CACHE BOOL "enable tests")

set(IPP_ENABLE_FPE FALSE CACHE BOOL "enable floating point exceptions debugging feature")
if(IPP_ENABLE_FPE)
	add_definitions(-DENABLE_FPE)
endif()

set(IPP_ENABLE_CORE_DUMP FALSE CACHE BOOL "enable core dump")
if(IPP_ENABLE_CORE_DUMP)
	message(STATUS "enable core dumps is on")
	add_definitions(-DIPP_ENABLE_CORE_DUMP)
endif()

set(IPP_ENABLE_BENCH_BARRIER FALSE CACHE BOOL "enable mpi barriers for bench timings")
if(IPP_ENABLE_BENCH_BARRIER)
	add_definitions(-DIPP_ENABLE_BENCH_BARRIER)
endif()

set(IPP_MPI_CHECK_SYNC FALSE CACHE BOOL "enable mpi sync checks")
if(IPP_MPI_CHECK_SYNC)
	add_definitions(-DENABLE_MPI_CHECK_SYNC)
endif()

set(IPP_DISTANCE_TRANSFORM_TEST FALSE CACHE BOOL "enable distance transform test")


set(IPP_DEBUG FALSE CACHE BOOL "enable debug outputs")
if(IPP_DEBUG)
	add_definitions(-DIPP_DEBUG)
endif()

set(IPP_DEBUG_ENABLED FALSE CACHE BOOL "enable debug outputs of enabled state")
if(IPP_DEBUG_ENABLED)
	add_definitions(-DIPP_DEBUG_ENABLED)
endif()

set(PALABOS_ROOT_DIR $ENV{PALABOS_ROOT_DIR} CACHE PATH "palabos install root")
set(PHREEQCRM_ROOT_DIR $ENV{PHREEQCRM_ROOT_DIR} CACHE PATH "PhreeqcRM install root")
set(JSONCPP_ROOT_DIR $ENV{JSONCPP_ROOT_DIR} CACHE PATH "jsoncpp install root")

message(STATUS "palabos root: " ${PALABOS_ROOT_DIR})
message(STATUS "phreeqcrm root: " ${PHREEQCRM_ROOT_DIR})
message(STATUS "jsoncpp root: " ${JSONCPP_ROOT_DIR})

find_package(Palabos REQUIRED)
find_package(PhreeqcRM REQUIRED)


set(EXAMPLE_INPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/examples/input" )
set(EXAMPLE_DATABASES_PATH "${CMAKE_CURRENT_SOURCE_DIR}/examples/databases" )
file(GLOB RES_FILES ${EXAMPLE_INPUT_PATH}/* ${EXAMPLE_DATABASES_PATH}/*)

add_subdirectory(lib)
add_subdirectory(examples)
