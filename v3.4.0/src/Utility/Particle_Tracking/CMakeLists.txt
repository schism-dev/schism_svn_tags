add_executable(ptrack_selfe2 ptrack_selfe2.f90)

add_dependencies(utility ptrack_selfe2)
target_link_libraries(ptrack_selfe2 utillib)


