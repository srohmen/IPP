set(target_name ${PROJECT_LIBRARY_NAME}Chemistry)

file(GLOB sources
	include/*
	src/*
)

add_library(${target_name} SHARED ${sources})

target_include_directories(${target_name} PUBLIC include)


# TODO: migrate the target dependencies to consume targets, not variables
target_include_directories(${target_name}
	PRIVATE
		${PHREEQCRM_INCLUDE_DIR}
		${VTK_INCLUDE_DIRS}
)

find_package(GSL REQUIRED)

target_link_libraries(${target_name}
	PUBLIC
		${PROJECT_LIBRARY_NAME}Common
	PRIVATE
		Boost::mpi
		Boost::iostreams
		${PHREEQCRM_LIBRARIES}
		ZLIB::ZLIB
		GSL::gsl
		${VTK_LIBRARIES}
)


install(TARGETS ${target_name}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)

