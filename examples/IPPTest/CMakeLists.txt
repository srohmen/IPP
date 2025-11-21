set(target_name ${PROJECT_LIBRARY_NAME}Test)

aux_source_directory(. sources)
add_executable(${target_name} ${sources})

target_link_libraries(${target_name}
	PRIVATE
		${PROJECT_LIBRARY_NAME}
		${PROJECT_LIBRARY_NAME}IO
		${PROJECT_LIBRARY_NAME}Visualization
		Boost::program_options
)


add_custom_command(
	TARGET ${target_name} POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E create_symlink ${EXAMPLE_INPUT_PATH} "input"
	COMMAND ${CMAKE_COMMAND} -E create_symlink ${EXAMPLE_DATABASES_PATH} "databases"
)

install(TARGETS ${target_name} RUNTIME DESTINATION bin)
