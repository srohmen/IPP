set(target_name ${PROJECT_LIBRARY_NAME}Common)

file(GLOB sources
	include/*
	src/*)
add_library(${target_name} SHARED ${sources})

target_include_directories(${target_name} PUBLIC include)

target_link_libraries(${target_name}
	PUBLIC
		${PROJECT_LIBRARY_NAME}Tools
		${PROJECT_LIBRARY_NAME}Utils
)

install(TARGETS ${target_name}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)

