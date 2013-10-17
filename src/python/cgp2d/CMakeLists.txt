include_directories(
    ${PYTHON_INCLUDE_PATH}
    ${PYTHON_NUMPY_INCLUDE_DIR}
)


add_library(_cgp2d SHARED
    py_cgp2d.cxx
    py_cell0.cxx
    py_cell1.cxx
    py_cell2.cxx
    py_cell0vec.cxx
    py_cell1vec.cxx
    py_cell2vec.cxx
    py_main.cxx
)
target_link_libraries(_cgp2d
    ${Boost_PYTHON_LIBRARIES}
    ${VIGRA_IMPEX_LIBRARY}
    ${VIGRA_NUMPY_CORE_LIBRARY}
    ${VIGRA_NUMPY_IMPEX_LIBRARY}
)
set_target_properties(_cgp2d PROPERTIES
    PREFIX ""
)
get_property(location TARGET _cgp2d PROPERTY LOCATION)
add_custom_command(TARGET _cgp2d
    POST_BUILD        
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        ${location}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

