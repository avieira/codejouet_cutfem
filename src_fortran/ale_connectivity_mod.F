      MODULE ALE_CONNECTIVITY_MOD
      
        IMPLICIT NONE
#include      "my_real.inc"

!     ************    !
!     Connectivity
!     ************    !
      TYPE, PUBLIC :: T_CONNECTIVITY
        INTEGER, DIMENSION(:), ALLOCATABLE :: IAD_CONNECT
        INTEGER, DIMENSION(:), ALLOCATABLE :: CONNECTED
        INTEGER, DIMENSION(:), ALLOCATABLE :: TYPE
      END TYPE T_CONNECTIVITY

!     Extended type : + iface
      TYPE, PUBLIC, EXTENDS(T_CONNECTIVITY) :: T_CONNECTIVITY_EXT1
        INTEGER, DIMENSION(:), ALLOCATABLE :: IFACE2
      END TYPE T_CONNECTIVITY_EXT1


!     ******************    !
!     Extended local ids    !
!     ******************    !
      TYPE :: T_IDGLOB
        INTEGER, DIMENSION(:), ALLOCATABLE :: ID
        INTEGER, DIMENSION(:), ALLOCATABLE :: UID
      END TYPE T_IDGLOB

!     ****************     !
!     ALE connectivity     !
!     ****************     !
      TYPE, PUBLIC :: T_ALE_CONNECTIVITY
!     node-node, node-element, element-element
        TYPE(T_CONNECTIVITY) :: NN_CONNECT, NE_CONNECT
        TYPE(T_CONNECTIVITY_EXT1) :: EE_CONNECT
        LOGICAL :: HAS_NN_CONNECT = .FALSE.
        LOGICAL :: HAS_NE_CONNECT = .FALSE.
        LOGICAL :: HAS_EE_CONNECT = .FALSE.
        LOGICAL :: HAS_IDGLOB = .FALSE.
        LOGICAL :: NALE_ALREADY_COMPUTED = .FALSE.
        INTEGER, DIMENSION(:), ALLOCATABLE :: NALE
        TYPE(T_IDGLOB) :: IDGLOB
      END TYPE T_ALE_CONNECTIVITY
      END MODULE ALE_CONNECTIVITY_MOD