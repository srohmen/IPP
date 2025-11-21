set(target_name ${PROJECT_LIBRARY_NAME}Utils)

file(GLOB sources
	include/*
	src/*
)
add_library(${target_name} SHARED ${sources})

target_include_directories(${target_name} PUBLIC include)

#required for mpimanager
target_link_libraries(${target_name} PUBLIC Boost::mpi)

install(TARGETS ${target_name}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)

