set(target_name ${PROJECT_LIBRARY_NAME}Transport)

file(GLOB sources
	include/*
	src/*
)
add_library(${target_name} SHARED ${sources})

target_include_directories(${target_name} PUBLIC include)

# TODO: migrate to cmake target approach
target_include_directories(${target_name} PRIVATE
	${PALABOS_INCLUDE_DIR}
	${VTK_INCLUDE_DIRS}
)

target_link_libraries(${target_name}
	PRIVATE
		${PROJECT_LIBRARY_NAME}Common
		${PALABOS_LIBRARIES}
		Boost::filesystem
		Boost::mpi
		${VTK_LIBRARIES}
)

option(IPP_ENABLE_MULTI_PRECISION "enable mpfr multi precision")
if(IPP_ENABLE_MULTI_PRECISION)
	message(STATUS "enabling mpfr multi precision")
	target_compile_definitions(${target_name} PRIVATE IPP_MULTI_PRECISION)

	find_package(MPFR REQUIRED)
	message("MPFR libraries: ${MPFR_LIBRARIES}")
	target_link_libraries(${target_name} PRIVATE ${MPFR_LIBRARIES})
endif()

install(TARGETS ${target_name}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)

