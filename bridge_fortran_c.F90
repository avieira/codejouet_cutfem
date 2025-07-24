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
  TYPE(t_ale_connectivity) :: ALE_CONNECT
  type(grid2D_struct_multicutcell), dimension(:, :), allocatable :: grid
  my_real :: threshold, dt, gamma
  integer :: sign
  my_real, dimension(:,:), allocatable :: vely
  my_real, dimension(:,:), allocatable :: velz
  my_real, dimension(:,:), allocatable :: rho
  my_real, dimension(:,:), allocatable :: p
  integer :: i

  N2D = 1
  NUMELQ = 4
  NUMELTG = 0
  NUMNOD = 9

  call launch_grb()

  allocate(IXQ(6,NUMELQ))
  allocate(IXTG(5,NUMELTG))
  allocate(ITAB(1))
  allocate(X(3, NUMNOD))
  allocate(ALE_CONNECT%ee_connect%iad_connect(NUMELQ + 1))
  allocate(ALE_CONNECT%ee_connect%connected(4*NUMELQ))
  allocate(grid(NUMELQ, 2))
  allocate(vely(NUMELQ, 2))
  allocate(velz(NUMELQ, 2))
  allocate(rho(NUMELQ, 2))
  allocate(p(NUMELQ, 2))
  allocate(xp(4))
  allocate(yp(4))

  threshold = 0.5
  dt = 1e-0
  gamma = 1.4
  sign = 1

  X(:, 1) = (/0., 0., 0./)
  X(:, 2) = (/0., 1., 0./)
  X(:, 3) = (/0., 2., 0./)
  X(:, 4) = (/0., 0., 1./)
  X(:, 5) = (/0., 1., 1./)
  X(:, 6) = (/0., 2., 1./)
  X(:, 7) = (/0., 0., 2./)
  X(:, 8) = (/0., 1., 2./)
  X(:, 9) = (/0., 2., 2./)

  IXQ(2:5, 1) = (/1, 2, 5, 4/)
  IXQ(2:5, 2) = (/2, 3, 6, 5/)
  IXQ(2:5, 3) = (/4, 5, 8, 7/)
  IXQ(2:5, 4) = (/5, 6, 9, 8/)

  ALE_CONNECT%ee_connect%iad_connect = (/1, 5, 9, 13, 17/)
  ALE_CONNECT%ee_connect%connected(1:4) = (/-1, 2, 3, -1/)
  ALE_CONNECT%ee_connect%connected(5:8) = (/-1, -1, 4, 1/)
  ALE_CONNECT%ee_connect%connected(9:12) = (/1, 4, -1, -1/)
  ALE_CONNECT%ee_connect%connected(13:16) = (/2, -1, -1, 3/)

  vely(:, :) = 0
  velz(:, :) = 0
  p(:, :) = 1
  rho(:, :) = 1

  xp(:) = (/0.5, 1.5, 1.5, 0.5/)
  yp(:) = (/0.5, 0.5, 1.5, 1.5/)

  call build_clipped_from_pts_fortran(xp, yp, 4)

  do i=1,NUMELQ
    grid(i, :) = makegrid(NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, i)
  end do

  call update_fluid(N2D, NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, ITAB, ALE_CONNECT, &
                          grid, vely, velz, rho, p, gamma, dt, threshold, sign)

  call end_grb()
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