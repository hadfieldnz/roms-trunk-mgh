      SUBROUTINE ana_respiration (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets respiration rate for hypoxia using analytical  !
!  expressions.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_ocean
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#include "tile.h"
!
      CALL ana_respiration_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           OCEAN(ng) % respiration)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(30)=MyFile
      END IF
!
      RETURN
      END SUBROUTINE ana_respiration
!
!***********************************************************************
      SUBROUTINE ana_respiration_tile (ng, tile, model,                 &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 respiration)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_biology
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: respiration(LBi:,LBj:,:)
#else
      real(r8), intent(out) :: respiration(LBi:UBi,LBj:UBj,N(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set respiration rate (1/day).
!-----------------------------------------------------------------------
!
#if defined CHESAPEAKE_1TERM
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            respiration(i,j,k)=ResRate(ng)
            IF (((i.ge.72).and.(j.le.35)).or.                           &
     &          (((i.ge.61).and.(i.le.71)).and.                         &
     &           ((j.ge. 6).and.(j.le.26)))) THEN
              respiration(i,j,k)=0.0_r8
            END IF
          END DO
        END DO
      END DO
#else
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
	    respiration(i,j,k)=ResRate(ng)
          END DO
        END DO
      END DO
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          respiration)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    respiration)
#endif
!
      RETURN
      END SUBROUTINE ana_respiration_tile
