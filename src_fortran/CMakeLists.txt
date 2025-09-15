set(all_f_files  
              ale_connectivity_mod.F
              grid_struct_multicutcell_mod.F90  
              integer_linked_list.F90  
              multicutcell_solver_mod.F90  
              polygon_mod.F90  
              riemann_solver_mod.F90
              CACHE INTERNAL "")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
set_property(GLOBAL PROPERTY mylib_f_module_dir ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(mylib_f SHARED "${all_f_files}")