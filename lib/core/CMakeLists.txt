file(GLOB sources
	include/*
	src/*)
add_library(${PROJECT_LIBRARY_NAME} SHARED ${sources})
target_include_directories(${PROJECT_LIBRARY_NAME} PUBLIC include)

target_link_libraries(${PROJECT_LIBRARY_NAME}
	PUBLIC
		${PROJECT_LIBRARY_NAME}Chemistry
		${PROJECT_LIBRARY_NAME}Transport
		Boost::mpi
		Boost::iostreams
)

install(TARGETS ${PROJECT_LIBRARY_NAME}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)

