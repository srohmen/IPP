set(target_name ${PROJECT_LIBRARY_NAME}Visualization)

file(GLOB sources
	include/*
	src/*)
add_library(${target_name} SHARED ${sources})

target_include_directories(${target_name} PUBLIC include)

target_link_libraries(${target_name} PRIVATE
	${PROJECT_LIBRARY_NAME}
	${PROJECT_LIBRARY_NAME}Common
)

if(IPP_ENABLE_RENDERING)
	include(${VTK_USE_FILE})
	message("VTK_LIBRARIES: " ${VTK_LIBRARIES})
	target_include_directories(${target_name} PRIVATE ${VTK_INCLUDE_DIRS} )
	target_link_libraries(${target_name} PRIVATE ${VTK_LIBRARIES})
else()
	add_definitions(-DNO_RENDERING)
endif()

install(TARGETS ${target_name}
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION include)
