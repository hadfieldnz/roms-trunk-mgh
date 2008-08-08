      SUBROUTINE ana_fsobc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets free-surface open boundary conditions using       !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_fsobc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME( 6)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_fsobc
!
!***********************************************************************
      SUBROUTINE ana_fsobc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Free-surface open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=0.0_r8
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=0.0_r8
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%zeta_south(i)=0.0_r8
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%zeta_north(i)=0.0_r8
        END DO
      END IF
#else
      ana_fsobc.h: No values provided for BOUNDARY(ng)%zeta_xxxx.
#endif

      RETURN
      END SUBROUTINE ana_fsobc_tile
