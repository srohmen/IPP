get_filename_component(target_name ${CMAKE_CURRENT_LIST_DIR} NAME)

aux_source_directory(. SRC_LIST)
add_executable(${target_name} ${SRC_LIST})
target_link_libraries(${target_name}
	PRIVATE
		${PROJECT_LIBRARY_NAME}Models_private
)
