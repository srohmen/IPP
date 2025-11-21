get_subdirs(sub_dirs ${CMAKE_CURRENT_SOURCE_DIR})

list(REMOVE_ITEM sub_dirs
	phreeqc
	input
	databases
)

foreach(sub_dir ${sub_dirs})
	message(STATUS "adding example: ${sub_dir}")
	add_subdirectory(${sub_dir})
endforeach()


