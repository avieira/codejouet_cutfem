#include "my_real.inc"

program main
  use multicutcell_solver_mod
  use polygon_mod
  use grid2D_struct_multicutcell_mod
  use ALE_CONNECTIVITY_MOD

  implicit none

  integer :: N2D, NUMELQ, NUMELTG, NUMNOD
  integer, dimension(:,:), allocatable :: IXQ, IXTG
  integer, dimension(:), allocatable :: ITAB
  my_real, dimension(:,:), allocatable :: X
  my_real, dimension(:), allocatable :: xp, yp
  integer(kind=8), dimension(:), allocatable :: limits_polygon
  TYPE(t_ale_connectivity) :: ALE_CONNECT
  type(grid2D_struct_multicutcell), dimension(:, :), allocatable :: grid
  my_real :: threshold, dt, gamma
  integer :: sign
  my_real, dimension(:,:), allocatable :: vely
  my_real, dimension(:,:), allocatable :: velz
  my_real, dimension(:,:), allocatable :: rho
  my_real, dimension(:,:), allocatable :: p
  integer :: i, nb_points_polygon, nb_polygons

  !N2D = 1
  !NUMELQ = 4
  !NUMELTG = 0
  !NUMNOD = 9
  N2D = 1
  NUMELQ = 8
  NUMELTG = 0
  NUMNOD = 15
  !nb_points_polygon = 4
  !nb_polygons = 1
  !nb_points_polygon = 10
  !nb_polygons = 1
  !nb_points_polygon = 8
  !nb_polygons = 2;
  !nb_points_polygon = 9
  !nb_polygons = 2;
  nb_points_polygon = 9
  nb_polygons = 1

  allocate(IXQ(6,NUMELQ))
  allocate(IXTG(5,NUMELTG))
  allocate(ITAB(1))
  allocate(X(3, NUMNOD))
  allocate(ALE_CONNECT%ee_connect%iad_connect(NUMELQ + 1))
  allocate(ALE_CONNECT%ee_connect%connected(4*NUMELQ))
  allocate(vely(NUMELQ, 2))
  allocate(velz(NUMELQ, 2))
  allocate(rho(NUMELQ, 2))
  allocate(p(NUMELQ, 2))
  allocate(xp(nb_points_polygon))
  allocate(yp(nb_points_polygon))
  allocate(limits_polygon(nb_polygons+1))

  threshold = 0.5
  dt = 1e-0
  gamma = 1.4
  sign = 1

  !! Grid definition !!
  !X(:, 1) = (/0., 0., 0./)
  !X(:, 2) = (/0., 1., 0./)
  !X(:, 3) = (/0., 2., 0./)
  !X(:, 4) = (/0., 0., 1./)
  !X(:, 5) = (/0., 1., 1./)
  !X(:, 6) = (/0., 2., 1./)
  !X(:, 7) = (/0., 0., 2./)
  !X(:, 8) = (/0., 1., 2./)
  !X(:, 9) = (/0., 2., 2./)

  !IXQ(2:5, 1) = (/1, 2, 5, 4/)
  !IXQ(2:5, 2) = (/2, 3, 6, 5/)
  !IXQ(2:5, 3) = (/4, 5, 8, 7/)
  !IXQ(2:5, 4) = (/5, 6, 9, 8/)

  !ALE_CONNECT%ee_connect%iad_connect = (/1, 5, 9, 13, 17/)
  !ALE_CONNECT%ee_connect%connected(1:4) = (/-1, 2, 3, -1/)
  !ALE_CONNECT%ee_connect%connected(5:8) = (/-1, -1, 4, 1/)
  !ALE_CONNECT%ee_connect%connected(9:12) = (/1, 4, -1, -1/)
  !ALE_CONNECT%ee_connect%connected(13:16) = (/2, -1, -1, 3/)

  !xp(:) = (/0.5, 1.5, 1.5, 0.5/)
  !yp(:) = (/0.5, 0.5, 1.5, 1.5/)

  X(:, 1) = (/0., 0., 0./)
  X(:, 2) = (/0., 1., 0./)
  X(:, 3) = (/0., 2., 0./)
  X(:, 4) = (/0., 3., 0./)
  X(:, 5) = (/0., 4., 0./)
  X(:, 6) = (/0., 0., 1./)
  X(:, 7) = (/0., 1., 1./)
  X(:, 8) = (/0., 2., 1./)
  X(:, 9) = (/0., 3., 1./)
  X(:,10) = (/0., 4., 1./)
  X(:,11) = (/0., 0., 2./)
  X(:,12) = (/0., 1., 2./)
  X(:,13) = (/0., 2., 2./)
  X(:,14) = (/0., 3., 2./)
  X(:,15) = (/0., 4., 2./)

  IXQ(2:5, 1) = (/1,  2,  7,  6/)
  IXQ(2:5, 2) = (/2,  3,  8,  7/)
  IXQ(2:5, 3) = (/3,  4,  9,  8/)
  IXQ(2:5, 4) = (/4,  5, 10,  9/)
  IXQ(2:5, 5) = (/6,  7, 12, 11/)
  IXQ(2:5, 6) = (/7,  8, 13, 12/)
  IXQ(2:5, 7) = (/8,  9, 14, 13/)
  IXQ(2:5, 8) = (/9, 10, 15, 14/)

  ALE_CONNECT%ee_connect%iad_connect = (/1, 5, 9, 13, 17, 21, 25, 29, 33/)
  ALE_CONNECT%ee_connect%connected(1:4) = (/-1, 2, 5, -1/)
  ALE_CONNECT%ee_connect%connected(5:8) = (/-1, 3, 6, 1/)
  ALE_CONNECT%ee_connect%connected(9:12) = (/-1, 4, 7, 2/)
  ALE_CONNECT%ee_connect%connected(13:16) = (/-1, -1, 8, 3/)
  ALE_CONNECT%ee_connect%connected(17:20) = (/1, 6, -1, -1/)
  ALE_CONNECT%ee_connect%connected(21:24) = (/2, 7, -1, 5/)
  ALE_CONNECT%ee_connect%connected(25:28) = (/3, 8, -1, 6/)
  ALE_CONNECT%ee_connect%connected(29:32) = (/4, -1, -1, 7/)

  !! Polygon definition : list of points, then limits of points in case of several polygons defined !!
  !xp(:) = (/0.1, 0.5, 1.5, 2.5, 3.5, 3.9, 3.5, 2.5, 1.5, 0.5/)
  !yp(:) = (/0.9, 0.5, 0.5, 0.5, 0.5, 1.1, 1.5, 1.5, 1.5, 1.5/)
  !limits_polygon(:) = (/1, 10/)
  !xp(:) = (/1.0, 1.9, 1.9, 1.0, 2.1, 2.5, 3.0, 2.5/)
  !yp(:) = (/0.5, 0.5, 1.5, 1.5, 1.0, 0.5, 1.0, 1.5/)
  !limits_polygon(:) = (/1, 4, 8/)
  !xp(:) = (/1.0, 1.9, 1.9, 1.0,  2.1,  2.5,  3.0,  2.1, 2.4/)
  !yp(:) = (/0.5, 0.5, 1.5, 1.5, 0.75, 0.75, 1.25, 1.25, 1.0/)
  !limits_polygon(:) = (/1, 4, 9/)
  xp(:) = (/1.0, 1.2, 2.8, 3.0, 2.8, 2.5, 2.0, 1.5, 1.2/)
  yp(:) = (/1.0, 0.5, 0.5, 1.0, 1.5, 0.6, 1.0, 0.6, 1.5/)
  limits_polygon(:) = (/1, nb_points_polygon/)

  vely(:, :) = 0.
  velz(:, :) = 0.
  p(:, :) = 1
  rho(:, :) = 1

  call initialize_solver(N2D, NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, ITAB, ALE_CONNECT, &
                              grid, xp, yp, limits_polygon)

  call update_fluid(N2D, NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, ITAB, ALE_CONNECT, &
                          grid, vely, velz, rho, p, gamma, dt, threshold, sign)

  do i=1,4
    write(*,*) "grid(i, :)%lambdan_per_cell = ", grid(i, :)%lambdan_per_cell
    write(*,*) "grid(i, :)%lambdanp1_per_cell = ", grid(i, :)%lambdanp1_per_cell
    write(*,*) "grid(i, 1)%normal_intern_face_space = ", grid(i, 1)%normal_intern_face_space
    write(*,*) "grid(i, 2)%normal_intern_face_space = ", grid(i, 2)%normal_intern_face_space
    write(*,*) "grid(i, 1)%lambda_per_edge = ", grid(i, 1)%lambda_per_edge
    write(*,*) "grid(i, 2)%lambda_per_edge = ", grid(i, 2)%lambda_per_edge
    write(*,*)
  end do


  call deallocate_solver()

  deallocate(IXQ)
  deallocate(IXTG)
  deallocate(ITAB)
  deallocate(X)
  deallocate(ALE_CONNECT%ee_connect%iad_connect)
  deallocate(ALE_CONNECT%ee_connect%connected)
  deallocate(vely)
  deallocate(velz)
  deallocate(rho)
  deallocate(p)
  deallocate(xp)
  deallocate(yp)
end program main


!program main
!  use polygon_mod, only: Point2D, Point3D
!
!  implicit none
!  integer(kind=4), parameter :: n_x = 4
!  integer(kind=4), parameter :: n_y = 4
!  real (kind=8) :: x_v(n_x) = (/0.0, 1.0, 2.0, 3.0/)
!  real (kind=8) :: y_v(n_y) = (/0.0, 1.0, 2.0, 3.0/)
!  integer(kind=4), parameter :: n_x_clipped = 2
!  integer(kind=4), parameter :: n_y_clipped = 2
!  real (kind=8) :: x_v_clipped(n_x_clipped) = (/0.5,2.5/)
!  real (kind=8) :: y_v_clipped(n_y_clipped) = (/0.5,2.5/)
!  real (kind=8), parameter :: dt = 1.0
!  integer(kind=8) :: p_clipper2D, p_clipped2D, polygon2d_from_vertices_fortran !kind=8 mandatory for 64 bits architecture
!  integer(kind=8), parameter :: nb_regions = 2
!  integer(kind=8), parameter :: nb_edges_per_cell = 4
!  integer(kind=8), parameter :: nb_vertices_per_cell = 4
!  integer(kind=8), parameter :: nb_cells = (n_x-1)*(n_y-1)
!  integer(kind=8), parameter :: nb_vertices = n_x*n_y
!  real(kind=8), dimension(nb_edges_per_cell*nb_regions) :: local_small_lambdas
!  integer(kind=8), dimension(nb_regions) :: small_lambdas
!  real(kind=8), dimension(nb_regions, nb_cells) :: big_lambda_n
!  real(kind=8), dimension(nb_regions, nb_cells) :: big_lambda_np1
!  type(Point2D), dimension(nb_vertices_per_cell, nb_cells) :: vec_move_solid
!  type(Point3D), dimension(nb_cells) :: mean_normal
!  logical(kind=1), dimension(nb_cells) :: is_narrowband
!  integer(kind=8) :: i, j, curr_cell
!  integer(kind=8) :: neigh_i, neigh_j, neigh_cell
!  integer(kind=8) :: grbinfo, grb_matrix_setelement, grb_matrix_getelement, grb_matrix_new_fp_fortran
!  real(kind=8), dimension(nb_cells) :: printable_lambda
!
!  do i=1,nb_cells
!    vec_move_solid(:,i)%y = 0.
!    vec_move_solid(:,i)%z = 0.
!  end do
!
!  call launch_grb()
!  do i=1,nb_regions
!    grbinfo = grb_matrix_new_fp_fortran(small_lambdas(i), nb_cells, nb_cells)
!  end do
!  p_clipped2D = polygon2d_from_vertices_fortran(x_v_clipped, n_x_clipped, y_v_clipped, n_y_clipped);
!  do j=1,n_y-1
!    do i=1,n_x-1
!      curr_cell = (n_x-1)*(j-1) + i
!      p_clipper2D = polygon2d_from_vertices_fortran(x_v(i:i+1), n_x_clipped, y_v(j:j+1), n_y_clipped);
!      call compute_lambdas2d_fortran(p_clipper2D, p_clipped2D, vec_move_solid(:, curr_cell), nb_vertices_per_cell, dt, &
!                            local_small_lambdas, &
!                            big_lambda_n(:, curr_cell), &
!                            big_lambda_np1(:, curr_cell), &
!                            mean_normal(curr_cell), is_narrowband(curr_cell))
!    
!      !Western neighbour
!      if (i == 1) then 
!        neigh_i = n_x-1
!      else
!        neigh_i = i-1
!      end if
!      neigh_j = j
!      neigh_cell = (n_x-1)*(neigh_j-1) + neigh_i
!      grbinfo = grb_matrix_setelement(small_lambdas(1), curr_cell, neigh_cell, local_small_lambdas(1))
!      grbinfo = grb_matrix_setelement(small_lambdas(2), curr_cell, neigh_cell, local_small_lambdas(2))
!
!      !Southern neighbour
!      if (j == 1) then 
!        neigh_j = n_y-1
!      else
!        neigh_j = j-1
!      end if
!      neigh_i = i
!      neigh_cell = (n_x-1)*(neigh_j-1) + neigh_i
!      grbinfo = grb_matrix_setelement(small_lambdas(1), curr_cell, neigh_cell, local_small_lambdas(3))
!      grbinfo = grb_matrix_setelement(small_lambdas(2), curr_cell, neigh_cell, local_small_lambdas(4))
!
!      !Eastern neighbour
!      if (i == n_x-1) then 
!        neigh_i = 1
!      else
!        neigh_i = i+1
!      end if
!      neigh_j = j
!      neigh_cell = (n_x-1)*(neigh_j-1) + neigh_i
!      grbinfo = grb_matrix_setelement(small_lambdas(1), curr_cell, neigh_cell, local_small_lambdas(5))
!      grbinfo = grb_matrix_setelement(small_lambdas(2), curr_cell, neigh_cell, local_small_lambdas(6))
!
!      !Northern neighbour
!      if (j == n_y-1) then 
!        neigh_j = 1
!      else
!        neigh_j = j+1
!      end if
!      neigh_i = i
!      neigh_cell = (n_x-1)*(neigh_j-1) + neigh_i
!      grbinfo = grb_matrix_setelement(small_lambdas(1), curr_cell, neigh_cell, local_small_lambdas(7))
!      grbinfo = grb_matrix_setelement(small_lambdas(2), curr_cell, neigh_cell, local_small_lambdas(8))
!
!      call polygon2D_free_fortran(p_clipper2D)
!    end do
!  end do
!
!
!  write(*,*) "small_lambdas(1) = "
!  do i=1,nb_cells
!    do j=1,nb_cells
!      grbinfo = grb_matrix_getelement(small_lambdas(1), i, j, printable_lambda(j))
!    end do
!    print *, printable_lambda
!  end do
!
!  write(*,*) "small_lambdas(2) = "
!  do i=1,nb_cells
!    do j=1,nb_cells
!      grbinfo = grb_matrix_getelement(small_lambdas(2), i, j, printable_lambda(j))
!    end do
!    print *, printable_lambda
!  end do
!
!  write(*,*) "big_lambda_n = "
!  do i=1,nb_cells
!    print *, big_lambda_n(:, i)
!  end do
!  write(*,*) "big_lambda_np1 = "
!  do i=1,nb_cells
!    print *, big_lambda_np1(:, i)
!  end do
!
!  write (*,*) "mean_normal = "
!  do i=1,nb_cells
!    write (*,*) mean_normal(i)%y, mean_normal(i)%z, mean_normal(i)%t
!  end do
!
!  write (*,*) "is_narrowband = ", is_narrowband(:)
!
!  call end_grb()
!end program