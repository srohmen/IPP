set(target_name ${PROJECT_LIBRARY_NAME}IO)

file(GLOB sources
	include/*
	src/*
)
add_library(${target_name} SHARED ${sources})

target_include_directories(${target_name} PUBLIC include)

# FIXME: use the cmake target approach
target_include_directories(${target_name} PRIVATE ${VTK_INCLUDE_DIRS} )
target_link_libraries(${target_name}
	PUBLIC
		${PROJECT_LIBRARY_NAME}
	PRIVATE
		${PROJECT_LIBRARY_NAME}Models
		Boost::program_options
		${VTK_LIBRARIES}
)


# link to json lib
find_package(jsoncpp)
if(jsoncpp_DIR)
	include(${jsoncpp_DIR}/jsoncppConfig.cmake)
	target_link_libraries(${target_name} PRIVATE jsoncpp_lib)
else()
	message("WARNING: could not find json cmake files. Trying fallback linking")

	if(JSONCPP_ROOT_DIR)
		target_link_libraries(${target_name} PRIVATE ${JSONCPP_ROOT_DIR}/lib/libjsoncpp.so)
	else()
		target_link_libraries(${target_name} PRIVATE libjsoncpp.so)

		# ugly hack: jsoncpp dropped cmake config file. TODO: write a proper cmake find script
		target_include_directories(${target_name} PRIVATE /usr/include/jsoncpp)
	endif()

endif()

install(TARGETS ${target_name}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)
