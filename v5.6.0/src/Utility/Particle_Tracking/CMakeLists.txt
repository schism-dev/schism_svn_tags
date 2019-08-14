add_executable(ptrack2 ptrack2.f90)

add_dependencies(utility ptrack2)
target_link_libraries(ptrack2 utillib)


