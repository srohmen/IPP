set(target_name ${PROJECT_LIBRARY_NAME}Tools)

file(GLOB sources
	include/*
	src/*
)
add_library(${target_name} SHARED ${sources})

target_include_directories(${target_name} PUBLIC include)

target_link_libraries(${target_name} PRIVATE
	${PROJECT_LIBRARY_NAME}Utils
	Boost::filesystem
	Boost::chrono)

install(TARGETS ${target_name}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)



