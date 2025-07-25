#include "my_real.inc"
#define NOT_FUSED 0 
#define FUSED 1
#define TARGET_FUSED 2
#define PROBLEMATIC 3
#define PROBLEMATIC_UNSOLVED 4

module multicutcell_solver_mod
  contains

  subroutine primal_to_conservative(gamma, vely, velz, rho, p, W)
    use grid2D_struct_multicutcell_mod
  
    implicit none
    my_real :: gamma
    my_real, intent(in) :: vely, velz, rho, p 
    type(ConservativeState2D), intent(out) :: W
  
    W%rho = rho
    W%rhovy = rho*vely
    W%rhovz = rho*velz
    W%rhoE = rho*(vely*vely + velz*velz) + p/(gamma-1)
  end subroutine primal_to_conservative

  subroutine conservative_to_primal(gamma, vely, velz, rho, p, W)
    use grid2D_struct_multicutcell_mod
  
    implicit none
    my_real :: gamma
    my_real, intent(out) :: vely, velz, rho, p 
    type(ConservativeState2D), intent(in) :: W
  
    rho = W%rho
    vely = W%rhovy/rho
    velz = W%rhovz/rho
    p =  (gamma - 1) * (W%rhoE - 0.5*(W%rhovy*W%rhovy + W%rhovz*W%rhovz)/W%rho)
  end subroutine conservative_to_primal
  
  subroutine adjacency_edge(ALE_CONNECT, i_cell, j_edge, other_cell, other_edge)
    use ALE_CONNECTIVITY_MOD

    implicit none
    TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
    integer(kind=8), intent(in) :: i_cell, j_edge
    integer(kind=8), intent(out) :: other_cell, other_edge

    integer(kind=8) :: IAD2, LGTH, J, IV

    IAD2 = ALE_CONNECT%ee_connect%iad_connect(i_cell)
    other_cell = ALE_CONNECT%ee_connect%connected(IAD2 + j_edge - 1)

    !Look for correct edge on the other side
    other_edge = -1
    if (other_cell>0) then
      IAD2 = ALE_CONNECT%ee_connect%iad_connect(other_cell)
      LGTH = ALE_CONNECT%ee_connect%iad_connect(other_cell+1) - IAD2
      do J=1,LGTH
        IV = ALE_CONNECT%ee_connect%connected(IAD2 + J - 1)
        if (IV == i_cell) then
          other_edge = J 
        end if
      end do
    end if
  end subroutine adjacency_edge

  subroutine hllc_flux(gamma, rhoL, rhoR, velyL, velyR, velzL, velzR, pL, pR, normalVec, rho, rhovy, rhovz, rhoE)
    use grid2D_struct_multicutcell_mod
    use polygon_mod

    implicit none

    my_real :: gamma, rhoL, rhoR, velyL, velyR, velzL, velzR, pL, pR
    type(Point2D) :: normalVec
    my_real :: rho, rhovy, rhovz, rhoE
    my_real :: velpyL, velpzL, velpyR, velpzR !Projected velocities
    my_real :: rhoEL, rhoER, SL, SR, S_star
    my_real :: rho_star, rhou_star, rhov_star, rhoE_star
    my_real :: aL, aR

    velpyL = velyL*normalVec%y + velzL*normalVec%z
    velpyR = velyR*normalVec%y + velzR*normalVec%z
    velpzL = -velyL*normalVec%y + velzL*normalVec%z
    velpzR = -velyR*normalVec%y + velzR*normalVec%z

    !Compute left and right speeds of sound
    aL = sqrt(gamma * pL / rhoL)
    aR = sqrt(gamma * pR / rhoR)
  
    !Compute conservative variables
    rhoEL = 0.5*rhoL*(velpyL*velpyL + velpzL*velpzL) + pL/(gamma-1) !energy
    rhoER = 0.5*rhoR*(velpyR*velpyR + velpzr*velpzr) + pR/(gamma-1)
  
    !Compute wave speeds
    SL = min(velpyL - aL, velpyR - aR)
    SR = max(velpyL + aL, velpyR + aR)
  
    !Compute the contact wave speed
    S_star = (pR - pL + rhoL * velpyL * (SL - velpyL) - rhoR * velpyR * (SR - velpyR)) / &
              (rhoL * (SL - velpyL) - rhoR * (SR - velpyR)) !Equation (37) in Toro slides
  
    !Compute the HLLC flux
    if (0 <= SL) then
      !return flux(UL, normalVec)
      rho   = rhoL * velpyL
      rhovy = rhoL * velyL * velpyL + pL * normalVec%y
      rhovz = rhoL * velzL * velpzL + pL * normalVec%z
      rhoE  = velpzL * (0.5 * rhoL * (velyL*velyL + velzL*velzL) + gamma / (gamma - 1) * pL)
    elseif (SL <= 0 .and. 0 <= S_star) then
      ! Compute the left and right star region states: Equation (39) in Toro slides
      rho_star  = rhoL * (SL - velpyL) / (SL - S_star)
      rhou_star = rho_star * S_star
      rhov_star = rho_star * velpzL
      rhoE_star = (SL - velpyL) / (SL - S_star) * (rhoEL + rhoL * (S_star - velpyL) * (S_star + pL / (rhoL * (SL - velpyL))))
    
      ! Compute the fluxes in the star region
      rho   = rhoL * velpyL + SL * (rho_star - rhoL)
      rhovy = rhoL * velyL * velpyL + pL * normalVec%y + &
              SL*((rhou_star - rhoL * velpyL) * normalVec%y - (rhov_star - rhoL * velpzL) * normalVec%z)
      rhovz = rhoL * velzL * velpzL + pL * normalVec%z + &
              SL*((rhou_star - rhoL * velpyL) * normalVec%z + (rhov_star - rhoL * velpzL) * normalVec%y)
      rhoE  = velpzL * (0.5 * rhoL * (velyL*velyL + velzL*velzL) + gamma / (gamma - 1) * pL) + SL*(rhoE_star - rhoEL)
    elseif (S_star <= 0 .and. 0 <= SR) then
      ! Compute the left and right star region states: Equation (39) in Toro slides
      rho_star  = rhoR * (SR - velpyR) / (SR - S_star)
      rhou_star = rho_star * S_star
      rhov_star = rho_star * velpzR
      rhoE_star = (SR - velpyR) / (SR - S_star) * (rhoER + rhoR * (S_star - velpyR) * (S_star + pR / (rhoR * (SR - velpyR))))
    
      ! Compute the fluxes in the star region
      rho   = rhoR * velpyR + SR * (rho_star - rhoR)
      rhovy = rhoR * velyR * velpyR + pR * normalVec%y + &
              SR*((rhou_star - rhoR * velpyR) * normalVec%y - (rhov_star - rhoR * velpzR) * normalVec%z)
      rhovz = rhoR * velzR * velpzR + pR * normalVec%z + &
              SR*((rhou_star - rhoR * velpyR) * normalVec%z + (rhov_star - rhoR * velpzR) * normalVec%y)
      rhoE  = velpzR * (0.5 * rhoR * (velyR*velyR + velzR*velzR) + gamma / (gamma - 1) * pR) + SR*(rhoE_star - rhoER)
    else
      rho   = rhoR * velpyR
      rhovy = rhoR * velyR * velpyR + pR * normalVec%y
      rhovz = rhoR * velzR * velpzR + pR * normalVec%z
      rhoE  = velpzR * (0.5 * rhoR * (velyR*velyR + velzR*velzR) + gamma / (gamma - 1) * pR)
    end if
  end subroutine hllc_flux
  
  subroutine fuse_cells(NUMELQ, NUMELTG, ALE_CONNECT, grid, threshold, final_target_cells, cell_type)
    use grid2D_struct_multicutcell_mod
    use integer_LL_mod
    use ALE_CONNECTIVITY_MOD

    implicit none
    integer :: NUMELQ, NUMELTG
    TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
    type(grid2D_struct_multicutcell), dimension(:, :) :: grid
    my_real :: threshold
    integer(kind=8), dimension(:, :) :: final_target_cells
    integer(kind=8), dimension(:, :) :: cell_type

    type(ptr_to_integer_LL_), dimension(:, :), allocatable :: target_cells
    logical :: narrowBand
    type(integer_LL_), pointer :: cells_to_be_merged
    type(integer_LL_), pointer :: copy_to_be_merged
    type(integer_LL_), pointer :: targ_c
    type(integer_LL_), pointer :: curr_ptr_int, curr_targ
    integer(kind=8) :: nb_cell, nb_regions, nb_edges
    integer(kind=8) :: i, k, j, c
    integer(kind=8) :: other_face, other_edge
    logical :: global_fusion_happened, fusion_happened
    integer(kind=8) :: size_list
    my_real :: max_lambda
    integer(kind=8) :: i_max_lambd, rand_cell

    nb_cell = size(grid, 1)
    nb_regions = 2

    allocate(target_cells(nb_cell, nb_regions))
    !allocate(targ_c)
    !allocate(cells_to_be_merged)
    nullify(targ_c)
    nullify(cells_to_be_merged)
    nullify(copy_to_be_merged)
    nullify(curr_ptr_int)
    nullify(curr_targ)

    do k = 1,nb_regions
      call integer_LL_destroy(cells_to_be_merged)
      do i = 1,nb_cell !First : detect all regions needing some fusion, and list the possible neighbours
        nullify(target_cells(i, k)%ptr)
        narrowBand = ((((0 < grid(i,k)%lambdan_per_cell / grid(i,k)%area) .and. &
                          (grid(i,k)%lambdan_per_cell / grid(i,k)%area < threshold)) .or. &
                        ((0 < grid(i,k)%lambdanp1_per_cell / grid(i, k)%area) .and. &
                          (grid(i,k)%lambdanp1_per_cell / grid(i, k)%area < threshold))) .and. &
                        grid(i, k)%is_narrowband)
        if (narrowBand) then
            call integer_LL_insert_unique(cells_to_be_merged, i)
            cell_type(i, k) = PROBLEMATIC_UNSOLVED
            nb_edges = max_nb_edges_in_cell(NUMELQ, NUMELTG)
            do j = 1,nb_edges
              if (grid(i, k)%lambda_per_edge(j) > 0) then
                call adjacency_edge(ALE_CONNECT, i, j, other_face, other_edge) 
                if  (other_face > 0) then
                  if ((grid(other_face, k)%lambdan_per_cell > 0.0) .or. (grid(other_face, k)%lambdanp1_per_cell > 0.0)) then
                    call integer_LL_insert_unique(target_cells(i, k)%ptr, other_face)
                  end if
                end if
              end if
            end do
        else
            cell_type(i, k) = NOT_FUSED
            call integer_LL_insert_after(target_cells(i, k)%ptr, i)
        end if
      end do

      !Do fusion: choose one neighbour among acceptable ones (acceptable /= possible), treat pathologic cases, do it again and again until nothing happens.
      global_fusion_happened = .true.
      do while (global_fusion_happened)
        global_fusion_happened = .false.

        fusion_happened = .true.
        do while (fusion_happened)
          fusion_happened = .false.
          call integer_LL_copy(cells_to_be_merged, copy_to_be_merged)

          !Iterate over copy_to_be_merged
          curr_ptr_int => copy_to_be_merged
          do while (associated(curr_ptr_int))
            !Build targ_c, a list of acceptable neighbours among possible ones.
            call integer_LL_destroy(targ_c)
            size_list = 0
            c = curr_ptr_int%val
            if (associated(target_cells(c, k)%ptr)) then
              curr_targ => target_cells(c, k)%ptr
              do while (associated(curr_targ))
                if ((cell_type(curr_targ%val, k) == NOT_FUSED) .or. (cell_type(curr_targ%val, k) == FUSED)) then
                  call integer_LL_insert_after(targ_c, curr_targ%val)
                  size_list = size_list + 1
                end if
                curr_targ => curr_targ%next
              end do
            end if

            !if targ_c is not empty: fuse with the biggest Lambda.
            if (size_list > 0) then
              curr_targ => targ_c
              max_lambda = grid(curr_targ%val, k)%lambdan_per_cell / grid(curr_targ%val, k)%area
              i_max_lambd = curr_targ%val
              curr_targ => curr_targ%next
              do while(associated(curr_targ)) !looking for maximum Lambda
                if (max_lambda < grid(curr_targ%val, k)%lambdan_per_cell / grid(curr_targ%val, k)%area) then
                  max_lambda = grid(curr_targ%val, k)%lambdan_per_cell / grid(curr_targ%val, k)%area
                  i_max_lambd = curr_targ%val
                end if
                curr_targ => curr_targ%next
              end do
              
              call integer_LL_destroy(target_cells(c, k)%ptr)
              call integer_LL_insert_after(target_cells(c, k)%ptr, i_max_lambd)
              
              call integer_LL_delete_value(cells_to_be_merged, c)
              cell_type(c, k) = FUSED
              cell_type(i_max_lambd, k) = FUSED
              fusion_happened = .true.
              global_fusion_happened = .true.
            end if
          
            curr_ptr_int => curr_ptr_int%next
          end do
        end do

        do while (integer_LL_size(cells_to_be_merged) > 0) !Problem : some cells can't be merged with a cell big enough.
          rand_cell = integer_LL_pop(cells_to_be_merged)
          call integer_LL_destroy(target_cells(rand_cell, k)%ptr)
          call integer_LL_insert_after(target_cells(rand_cell, k)%ptr, rand_cell)
          cell_type(rand_cell, k) = PROBLEMATIC

          fusion_happened = .true.
          do while (fusion_happened)
            fusion_happened = .false.
            call integer_LL_destroy(copy_to_be_merged)
            call integer_LL_copy(cells_to_be_merged, copy_to_be_merged)

            curr_ptr_int => copy_to_be_merged
            do while (associated(curr_ptr_int))
              !Build targ_c, a list of acceptable neighbours among possible ones.
              call integer_LL_destroy(targ_c)
              size_list = 0
              c = curr_ptr_int%val
              curr_targ => target_cells(c, k)%ptr
              do while (associated(curr_targ))
                if (cell_type(curr_targ%val, k) == PROBLEMATIC) then
                  call integer_LL_insert_after(targ_c, curr_targ%val)
                  size_list = size_list + 1
                end if
                curr_targ => curr_targ%next
              end do

              !if targ_c is not empty: fuse with the biggest Lambda.
              if (size_list > 0) then
                curr_targ => targ_c
                max_lambda = grid(curr_targ%val, k)%lambdan_per_cell / grid(curr_targ%val, k)%area
                i_max_lambd = curr_targ%val
                curr_targ => curr_targ%next
                do while(associated(curr_targ)) !looking for maximum Lambda
                  if (max_lambda < grid(curr_targ%val, k)%lambdan_per_cell / grid(curr_targ%val, k)%area) then
                    max_lambda = grid(curr_targ%val, k)%lambdan_per_cell / grid(curr_targ%val, k)%area
                    i_max_lambd = curr_targ%val
                  end if
                  curr_targ => curr_targ%next
                end do

                call integer_LL_destroy(target_cells(c, k)%ptr)
                call integer_LL_insert_after(target_cells(c, k)%ptr, i_max_lambd)
                call integer_LL_delete_value(cells_to_be_merged, c)
                cell_type(c, k) = PROBLEMATIC
                fusion_happened = .true.
              end if
              curr_ptr_int => curr_ptr_int%next
            end do
          end do
        end do

        do i = 1,nb_cell
            final_target_cells(i, k) = integer_LL_pop(target_cells(i, k)%ptr)
        end do

        do i = 1,nb_cell
          j = final_target_cells(i, k)
          do while (j /= final_target_cells(j, k))
              j = final_target_cells(j, k)
          end do
          final_target_cells(i, k) = j
          if ((i == j) .and. (cell_type(i, k) == FUSED)) then
              cell_type(i, k) = TARGET_FUSED
          end if
        end do
      end do
    end do

    !Deallocate all temp variables
    do k = 1,nb_regions 
      do i = 1,nb_cell
        call integer_LL_destroy(target_cells(i,k)%ptr)
      enddo
    enddo
    deallocate(target_cells)

    call integer_LL_destroy(targ_c)
    call integer_LL_destroy(cells_to_be_merged)
  end subroutine fuse_cells
  
  subroutine compute_fluxes(NUMELQ, NUMELTG, IXQ, IXTG, X, ALE_CONNECT, gamma, rho, vely, velz, p, fx)
    use grid2D_struct_multicutcell_mod
    use ALE_CONNECTIVITY_MOD

    implicit none

    !Dummy arguments
    integer, intent(in) :: NUMELQ, NUMELTG
    integer, dimension(:,:), intent(in) :: IXQ, IXTG
    my_real, dimension(:,:), intent(in) :: X
    TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
    my_real :: gamma
    my_real, dimension(:,:) :: vely
    my_real, dimension(:,:) :: velz
    my_real, dimension(:,:) :: rho
    my_real, dimension(:,:) :: p
    type(ConservativeFlux2D), dimension(:, :) :: fx

    !Local variables
    integer(kind=8) :: i, j, k
    integer(kind=8) :: nb_cell, nb_regions, nb_edges
    integer(kind=8) :: other_edge, other_face
    type(Point2D), dimension(4) :: normals

    nb_cell = size(vely, 1)
    nb_regions = size(vely, 2)

    do i=1,nb_cell
      call compute_normals(NUMELQ, NUMELTG, IXQ, IXTG, X, i, normals, nb_edges)
      do j=1,nb_edges
        call adjacency_edge(ALE_CONNECT, i, j, other_face, other_edge) 
        do k = 1,nb_regions
          fx(i, k)%rho(j) = 0.
          fx(i, k)%rhovy(j) = 0.
          fx(i, k)%rhovz(j) = 0.
          fx(i, k)%rhoE(j) = 0.
          !normalVector = grid.λ_per_edge[e, r]
          !norm_vec = norm(normalVector)
          if (other_face > -1) then
            call hllc_flux(gamma, rho(i, k), rho(other_face, k), &
                                vely(i, k), vely(other_face, k), velz(i, k), velz(other_face, k), &
                                p(i, k), p(other_face, k), &
                                normals(j), &
                                fx(i, k)%rho(j), fx(i, k)%rhovy(j), fx(i, k)%rhovz(j), fx(i, k)%rhoE(j))
          else !There is no neighbour!
            !Exact flux
            fx(i, k)%rho(j) = 0.
            fx(i, k)%rhovy(j) = p(i, k)*normals(j)%y
            fx(i, k)%rhovz(j) = p(i, k)*normals(j)%z
            fx(i, k)%rhoE(j) = 0.
          end if
        end do
      end do
    end do
  end subroutine compute_fluxes
  
  function max_nb_edges_in_cell(NUMELQ, NUMELTG)
    integer :: NUMELQ, NUMELTG
    integer(kind=8) ::max_nb_edges_in_cell

    max_nb_edges_in_cell = merge(4, 3, NUMELQ>0) ! NUMELQ > 0 ? 4 : 3;
  end function max_nb_edges_in_cell

  function max_nb_points_in_cell(NUMELQ, NUMELTG)
    integer :: NUMELQ, NUMELTG
    integer(kind=8) :: max_nb_points_in_cell

    max_nb_points_in_cell = merge(4, 3, NUMELQ>0) ! NUMELQ > 0 ? 4 : 3;
  end function max_nb_points_in_cell

  subroutine compute_close_cells(NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, grid)
    use polygon_mod
    use grid2D_struct_multicutcell_mod
  
    IMPLICIT NONE
  
    ! INPUT argument
    integer :: NUMELQ, NUMELTG, NUMNOD
    integer, dimension(:,:) :: IXQ, IXTG
    my_real, dimension(:,:) :: X
    ! IN/OUTPUT argument
    type(grid2D_struct_multicutcell), dimension(:, :) :: grid

    !Local variables
    logical, dimension(NUMNOD) :: is_narrowband_pt
    integer :: nb_cell, nb_regions, nb_pts_in_cell, i, k

    nb_cell = size(grid, 1)
    nb_regions = size(grid, 2)
    if (NUMELQ>0) then
      nb_pts_in_cell = 4
    else
      nb_pts_in_cell = 3
    end if
    !First : all points in contact with a cell in a narrowband are marked as true.
    is_narrowband_pt(:) = .false.
    do i = 1,nb_cell
      do k=2,2+nb_pts_in_cell-1
        is_narrowband_pt(IXQ(k, i)) = (is_narrowband_pt(IXQ(k, i)) .or. grid(i, 1)%is_narrowband)
      end do
    end do

    !Second : all cells that have at least one point marked true are close cells
    grid(:,:)%close_cells = .false.
    do i = 1,nb_cell
      do k=2,2+nb_pts_in_cell-1
        grid(i,1)%close_cells = (grid(i,1)%close_cells .or. is_narrowband_pt(IXQ(k, i)))
        grid(i,2:nb_regions)%close_cells = grid(i,1)%close_cells
      end do
    end do
  end subroutine compute_close_cells

  subroutine compute_lambdas(NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, grid, dt)
    use polygon_mod
    use grid2D_struct_multicutcell_mod
  
    IMPLICIT NONE
  
    ! INPUT argument
    integer :: NUMELQ, NUMELTG, NUMNOD
    integer, dimension(:,:) :: IXQ, IXTG
    my_real, dimension(:,:) :: X
    my_real :: dt
    ! IN/OUTPUT argument
    type(grid2D_struct_multicutcell), dimension(:, :) :: grid

    !Local variables
    integer(kind=8) :: nb_cell, nb_regions, nb_edges
    integer(kind=8) :: i, j, k
    my_real, dimension(:), allocatable :: ptr_lambdas_arr
    my_real, dimension(:), allocatable :: ptr_big_lambda_n, ptr_big_lambda_np1
    type(Point3D) :: mean_normal
    logical :: is_narrowband
    my_real, dimension(4) :: ys, zs

    nb_cell = size(grid, 1)
    nb_regions = size(grid, 2)

    allocate(ptr_lambdas_arr(4*nb_regions))
    allocate(ptr_big_lambda_n(nb_regions))
    allocate(ptr_big_lambda_np1(nb_regions))

    do i = 1,nb_cell
      if (grid(i,1)%close_cells) then
        if (NUMELQ>0) then
          nb_edges = 4
        else
          nb_edges = 3
        end if
        ys(1:nb_edges) = X(2, IXQ(2:2+nb_edges-1, i))
        zs(1:nb_edges) = X(3, IXQ(2:2+nb_edges-1, i))
        call build_grid_from_points_fortran(ys, zs, nb_edges) 
        call compute_lambdas2d_fortran(dt, ptr_lambdas_arr, ptr_big_lambda_n, ptr_big_lambda_np1, mean_normal, is_narrowband)

        do j=1,nb_regions
          grid(i,j)%lambdan_per_cell = ptr_big_lambda_n(j)
          grid(i,j)%lambdanp1_per_cell = ptr_big_lambda_np1(j)
          grid(i,j)%normal_intern_face_space%y = mean_normal%y
          grid(i,j)%normal_intern_face_space%z = mean_normal%z
          grid(i,j)%normal_intern_face_time = mean_normal%t
          grid(i,j)%is_narrowband = is_narrowband
          do k=1,nb_edges
            grid(i,j)%lambda_per_edge(k) = ptr_lambdas_arr((k-1)*nb_regions + j)
          end do
        end do
      end if
    end do

    deallocate(ptr_lambdas_arr)
    deallocate(ptr_big_lambda_n)
    deallocate(ptr_big_lambda_np1)

    call compute_close_cells(NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, grid)
  end subroutine compute_lambdas

  subroutine compute_normals(NUMELQ, NUMELTG, IXQ, IXTG, X, i_cell, normals, nb_normals)
    use polygon_mod

    IMPLICIT NONE

    integer(kind=8), intent(in) :: i_cell
    integer, intent(in) :: NUMELQ, NUMELTG
    integer, dimension(:,:), intent(in) :: IXQ, IXTG
    my_real, dimension(:,:), intent(in) :: X
    type(Point2D), dimension(4) :: normals
    integer(kind=8), intent(out) :: nb_normals

    integer(kind=8) :: j 
    type(Point2D) :: P1, P2, P3, P4
    my_real :: norm
    
    if (NUMELQ > 0) then
      nb_normals = 4

      P1%y = X(2, IXQ(2, i_cell))
      P1%z = X(3, IXQ(2, i_cell))
      
      P2%y = X(2, IXQ(3, i_cell))
      P2%z = X(3, IXQ(3, i_cell))

      P3%y = X(2, IXQ(4, i_cell))
      P3%z = X(3, IXQ(4, i_cell))

      P4%y = X(2, IXQ(5, i_cell))
      P4%z = X(3, IXQ(5, i_cell))

      normals(1)%y = (P2%z-P1%z)
      normals(1)%z =-(P2%y-P1%y)
      normals(2)%y = (P3%z-P2%z)
      normals(2)%z =-(P3%y-P2%y)
      normals(3)%y = (P4%z-P3%z)
      normals(3)%z =-(P4%y-P3%y)
      normals(4)%y = (P1%z-P4%z)
      normals(4)%z =-(P1%y-P4%y)
    else
      nb_normals = 3

      P1%y = X(2, IXTG(2, i_cell))
      P1%z = X(3, IXTG(2, i_cell))
      
      P2%y = X(2, IXTG(3, i_cell))
      P2%z = X(3, IXTG(3, i_cell))

      P3%y = X(2, IXTG(4, i_cell))
      P3%z = X(3, IXTG(4, i_cell))

      normals(1)%y = (P2%z-P1%z)
      normals(1)%z =-(P2%y-P1%y)
      normals(2)%y = (P3%z-P2%z)
      normals(2)%z =-(P3%y-P2%y)
      normals(3)%y = (P1%z-P3%z)
      normals(3)%z =-(P1%y-P3%y)
    end if

    do j=1,nb_normals
      !Normalize
      norm = sqrt(normals(j)%y*normals(j)%y + normals(j)%z*normals(j)%z)
      normals(j)%y = normals(j)%y / norm
      normals(j)%z = normals(j)%z / norm
    end do
  end subroutine compute_normals

  function find_ind_encompassing_cell(pt) result(i)
    use polygon_mod
    implicit none
    type(Point2D) :: pt
    integer(kind=8) :: i

    i = 1 !TODO Change this
  end function find_ind_encompassing_cell

  subroutine compute_vec_move_clipped(gamma, rho, vely, velz, p, nb_pts_clipped, &
                                    normalVecy, normalVecz, normalVecEdgey, normalVecEdgez, &
                                    vec_move_clippedy, vec_move_clippedz, pressure_edge)
    use polygon_mod
    use riemann_solver_mod

    implicit none
    my_real :: gamma
    my_real, dimension(:,:) :: rho
    my_real, dimension(:,:) :: vely
    my_real, dimension(:,:) :: velz
    my_real, dimension(:,:) :: p
    integer(kind=8)       :: nb_pts_clipped
    my_real, dimension(:) :: normalVecy
    my_real, dimension(:) :: normalVecz
    my_real, dimension(:) :: normalVecEdgey
    my_real, dimension(:) :: normalVecEdgez
    my_real, dimension(:) :: vec_move_clippedy
    my_real, dimension(:) :: vec_move_clippedz
    my_real, dimension(:) :: pressure_edge

    integer(kind=8) :: eL, eR, k, i
    type(Point2D) :: pt
    my_real :: rhoL, velyL, velzL, pL
    my_real :: rhoR, velyR, velzR, pR
    my_real :: us, vsL, vsR, ps
    my_real :: uLEdge, vLEdgeL, vLEdgeR, pEdgeL
    my_real :: uREdge, vREdgeL, vREdgeR, pEdgeR

    do k = 1,nb_pts_clipped
      call get_clipped_ith_vertex_fortran(k, pt, eR, eL) 
      if ((eR>-1) .and. (eL>-1)) then !point k is part of an edge
        i = find_ind_encompassing_cell(pt) !TODO
        rhoL = rho(i, 1)
        velyL = vely(i, 1)
        velzL = velz(i, 1)
        pL = p(i, 1)
        rhoR = rho(i, 2)
        velyR = vely(i, 2)
        velzR = velz(i, 2)
        pR = p(i, 2)

        call solve_riemann_problem(gamma, rhoL, rhoR, velyL, velyR, velzL, velzR, pL, pR, 2, &
                                    normalVecy(k), normalVecz(k), &
                                    us, vsL, vsR, ps)
      
        call solve_riemann_problem(gamma, rhoL, rhoR, velyL, velyR, velzL, velzR, pL, pR, 2, &
                                  normalVecEdgey(eL), normalVecEdgez(eL), uLEdge, vLEdgeL, vLEdgeR, pEdgeL)
        call solve_riemann_problem(gamma, rhoL, rhoR, velyL, velyR, velzL, velzR, pL, pR, 2, &
                                  normalVecEdgey(eR), normalVecEdgez(eR), uREdge, vREdgeL, vREdgeR, pEdgeR)
        pressure_edge(eL+1) = pEdgeL !+1 because C is 0-indexed and fortran is 1-indexed...
        pressure_edge(eR+1) = pEdgeR !+1 because C is 0-indexed and fortran is 1-indexed...
      
        vec_move_clippedy(k) = us * normalVecy(k) - 0.5 * (vsL + vsR) * normalVecz(k) !Choice of the mean of left and right tangential velocities, another choice could be made!
        vec_move_clippedz(k) = us * normalVecz(k) + 0.5 * (vsL + vsR) * normalVecy(k) !Choice of the mean of left and right tangential velocities, another choice could be made!
      !Transfer phi to fluid flux: TODO !
      else
          vec_move_clippedy(k) = 0.0
          vec_move_clippedz(k) = 0.0
      end if
    end do
  
    !return vec_move_clipped
  end subroutine compute_vec_move_clipped

  subroutine update_fluid(N2D, NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, ITAB, ALE_CONNECT, &
                          grid, vely, velz, rho, p, gamma, dt, threshold, sign)
    use polygon_mod
    use grid2D_struct_multicutcell_mod
    use ALE_CONNECTIVITY_MOD
  
    !Manuel : chercher analy pour les normales
    IMPLICIT NONE
  
    ! INPUT arguments
    integer :: N2D, NUMELQ, NUMELTG, NUMNOD
    integer, dimension(:,:) :: IXQ, IXTG
    integer, dimension(:) :: ITAB
    my_real, dimension(:,:) :: X
    TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
    type(grid2D_struct_multicutcell), dimension(:, :) :: grid
    my_real :: threshold, dt, gamma
    integer :: sign
    ! IN/OUTPUT arguments
    my_real, dimension(:,:) :: vely
    my_real, dimension(:,:) :: velz
    my_real, dimension(:,:) :: rho
    my_real, dimension(:,:) :: p
  
    !Local variables
    integer(kind=8) :: i, j, k
    integer(kind=8) :: nb_cell, nb_regions, nb_edges
    integer(kind=8) :: ind_targ
    type(Point2D) :: pressure_mean_normal
    my_real :: mean_normal_time
    !my_real :: lambdan
    my_real :: lambdanp1
    type(ConservativeState2D), dimension(:, :), allocatable :: W
    type(ConservativeState2D), dimension(:, :), allocatable :: dW
    type(ConservativeFlux2D), dimension(:, :), allocatable :: fx
    my_real, dimension(:,:), allocatable :: lambdan_prev
    !my_real, dimension(:), allocatable :: areas
    integer(kind=8), dimension(:, :), allocatable :: target_cells
    integer(kind=8), dimension(:, :), allocatable :: cell_type
    my_real, dimension(:), allocatable :: normalVecy
    my_real, dimension(:), allocatable :: normalVecz
    my_real, dimension(:), allocatable :: normalVecEdgey
    my_real, dimension(:), allocatable :: normalVecEdgez
    my_real :: min_pos_Se
    my_real, dimension(:), allocatable :: vec_move_clippedy
    my_real, dimension(:), allocatable :: vec_move_clippedz
    my_real, dimension(:), allocatable :: pressure_edge
    integer(kind=8) :: nb_pts_clipped
    integer(kind=8) :: nb_edge_clipped
    my_real :: dx
    my_real :: pressure_mean_normal_time
    type(Point2D) :: mean_normal
    integer(kind=2) :: odd_k
    
    dx = sqrt(minval(grid(:, 1)%area))
    nb_cell = NUMELQ+NUMELTG !size(vely, 1)
    nb_regions = size(vely, 2)
  
    call nb_pts_clipped_fortran(nb_pts_clipped)
    call nb_edge_clipped_fortran(nb_edge_clipped)
    allocate(lambdan_prev(nb_cell, nb_regions))
    allocate(fx(nb_cell, nb_regions))
    !allocate(areas(nb_cell))
    allocate(target_cells(nb_cell, nb_regions))
    allocate(cell_type(nb_cell, nb_regions))
    allocate(W(nb_cell, nb_regions))
    allocate(dW(nb_cell, nb_regions))
    allocate(normalVecy(nb_pts_clipped))
    allocate(normalVecz(nb_pts_clipped))
    allocate(vec_move_clippedy(nb_pts_clipped))
    allocate(vec_move_clippedz(nb_pts_clipped))
    allocate(pressure_edge(nb_edge_clipped))
    allocate(normalVecEdgey(nb_edge_clipped))
    allocate(normalVecEdgez(nb_edge_clipped))
  
    do i = 1,nb_cell
      do k = 1,nb_regions
        call primal_to_conservative(gamma, vely(i, k), velz(i, k), rho(i, k), p(i, k), W(i, k))
        lambdan_prev(i, k) = grid(i, k)%lambdanp1_per_cell
      end do 
    end do
  
    normalVecy(:) = 0.0
    normalVecz(:) = 0.0
    normalVecEdgey(:) = 0.0
    normalVecEdgez(:) = 0.0
    call compute_normals_clipped_fortran(normalVecy, normalVecz, nb_pts_clipped, &
                                          normalVecEdgey, normalVecEdgez, nb_edge_clipped, min_pos_Se)

    pressure_edge(:) = 0.0
    call compute_vec_move_clipped(gamma, rho, vely, velz, p, nb_pts_clipped, &
                                normalVecy, normalVecz, normalVecEdgey, normalVecEdgez, &
                                vec_move_clippedy, vec_move_clippedz, pressure_edge)
  
    call smooth_vel_clipped_fortran(vec_move_clippedy, vec_move_clippedz, min_pos_Se, dt)
    call update_clipped_fortran(vec_move_clippedy, vec_move_clippedz, dt, dx)
  
    call compute_lambdas(NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, grid, dt) 
    !TODO exchange lambdas between procs on neighbouring cells
    call fuse_cells(NUMELQ, NUMELTG, ALE_CONNECT, grid, threshold, target_cells, cell_type) 
    
    call compute_fluxes(NUMELQ, NUMELTG, IXQ, IXTG, X, ALE_CONNECT, gamma, rho, vely, velz, p, fx)
    
    !Compute right hand side for non-fused or target cells
    grid(:, :)%lambdanp1_per_cell_target = 0.
    grid(:, :)%lambdan_per_cell_target = 0.
    
    odd_k = 1
    do k=1,nb_regions
      do i = 1,nb_cell
        ind_targ = target_cells(i, k) !TODO problem here when parallel: what if the target cell is in an other region?
        dW(ind_targ, k)%rho = dW(ind_targ, k)%rho + grid(i, k)%lambdan_per_cell * W(i, k)%rho  !Add explicit term of Euler scheme
        dW(ind_targ, k)%rhovy = dW(ind_targ, k)%rhovy + grid(i, k)%lambdan_per_cell * W(i, k)%rhovy  
        dW(ind_targ, k)%rhovz = dW(ind_targ, k)%rhovz + grid(i, k)%lambdan_per_cell * W(i, k)%rhovz  
        dW(ind_targ, k)%rhoE = dW(ind_targ, k)%rhoE + grid(i, k)%lambdan_per_cell * W(i, k)%rhoE  
    
        !Add numerical fluxes
        nb_edges = max_nb_edges_in_cell(NUMELQ, NUMELTG)
        do j = 1,nb_edges
          dW(ind_targ, k)%rho = dW(ind_targ, k)%rho - dt * fx(i, k)%rho(j) * grid(i, k)%lambda_per_edge(j) 
          dW(ind_targ, k)%rhovy = dW(ind_targ, k)%rhovy - dt * fx(i, k)%rhovy(j) * grid(i, k)%lambda_per_edge(j) 
          dW(ind_targ, k)%rhovz = dW(ind_targ, k)%rhovz - dt * fx(i, k)%rhovz(j) * grid(i, k)%lambda_per_edge(j) 
          dW(ind_targ, k)%rhoE = dW(ind_targ, k)%rhoE - dt * fx(i, k)%rhoE(j) * grid(i, k)%lambda_per_edge(j) 
        end do
    
        !Add swept quantities
        mean_normal = grid(i, 1)%normal_intern_face_space
        mean_normal_time = grid(i, 1)%normal_intern_face_time
        pressure_mean_normal = grid(i, 1)%p_normal_intern_face_space
        pressure_mean_normal_time = grid(i, 1)%p_normal_intern_face_time

        dW(ind_targ, k)%rhovy = dW(ind_targ, 1)%rhovy + odd_k*pressure_mean_normal%y
        dW(ind_targ, k)%rhovz = dW(ind_targ, 1)%rhovz + odd_k*pressure_mean_normal%z
        dW(ind_targ, k)%rhoE  = dW(ind_targ, 1)%rhoE  - odd_k*pressure_mean_normal_time

        !Update size of fluid part
        grid(ind_targ, k)%lambdanp1_per_cell_target = grid(ind_targ, k)%lambdanp1_per_cell_target + grid(i, k)%lambdanp1_per_cell
        grid(ind_targ, k)%lambdan_per_cell_target = grid(ind_targ, k)%lambdan_per_cell_target + grid(i, k)%lambdan_per_cell
        !areas(ind_targ) = areas(ind_targ) + grid(i, 1)%area
      end do
      odd_k = -1
    end do
    
    !Update non-fused or target cells
    do k = 1,nb_regions
        do i = 1,nb_cell
          if (i == target_cells(i, k)) then
            lambdanp1 = grid(i, k)%lambdanp1_per_cell_target
            !lambdan = grid(i)%lambdan_per_cell_target(k)
            if (lambdanp1 > 0.5 * threshold * grid(i, k)%area) then
                W(i, k)%rho   = dW(i, k)%rho / lambdanp1
                W(i, k)%rhovy = dW(i, k)%rhovy / lambdanp1
                W(i, k)%rhovz = dW(i, k)%rhovz / lambdanp1
                W(i, k)%rhoE  = dW(i, k)%rhoE / lambdanp1
            end if
          end if
        end do
    end do
    
    !Broadcast value
    do i = 1,nb_cell
      do k = 1, nb_regions
        if (i /= target_cells(i, k)) then
          W(i, k)%rho   = W(target_cells(i, k), k)%rho
          W(i, k)%rhovy = W(target_cells(i, k), k)%rhovy
          W(i, k)%rhovz = W(target_cells(i, k), k)%rhovz
          W(i, k)%rhoE  = W(target_cells(i, k), k)%rhoE
        end if
        call conservative_to_primal(gamma, vely(i, k), velz(i, k), rho(i, k), p(i, k), W(i, k))
      end do
    end do
  
    !Deallocate local memory
    deallocate(lambdan_prev)
    deallocate(fx)
    !deallocate(areas)
    deallocate(target_cells)
    deallocate(cell_type)
    deallocate(W)
    deallocate(dW)
    deallocate(pressure_edge)
    deallocate(vec_move_clippedy)
    deallocate(vec_move_clippedz)
  
  end subroutine update_fluid
end module multicutcell_solver_mod

#undef NOT_FUSED 
#undef FUSED
#undef TARGET_FUSED
#undef PROBLEMATIC
#undef PROBLEMATIC_UNSOLVED