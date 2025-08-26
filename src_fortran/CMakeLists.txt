set(all_f_files  
              ale_connectivity_mod.F
              grid_struct_multicutcell_mod.F90  
              integer_linked_list.F90  
              multicutcell_solver_mod.F90  
              polygon_mod.F90  
              riemann_solver_mod.F90
              CACHE INTERNAL "")
add_library(mylib_f STATIC "${all_f_files}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")