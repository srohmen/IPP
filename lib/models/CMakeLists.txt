set(target_name ${PROJECT_LIBRARY_NAME}Models)

file(GLOB sources
	include/*
	src/*
)
add_library(${target_name} SHARED ${sources})

target_include_directories(${target_name} PUBLIC include)

target_link_libraries(${target_name}
	PUBLIC
		${PROJECT_LIBRARY_NAME}Common
	PRIVATE
		${PROJECT_LIBRARY_NAME}
		${PROJECT_LIBRARY_NAME}Utils
)

install(TARGETS ${target_name}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)


add_library(${target_name}_private INTERFACE)
target_link_libraries(${target_name}_private INTERFACE ${target_name})
target_include_directories(${target_name}_private INTERFACE src)

add_subdirectory(tests)
