#include "cppdefs.h"
      MODULE t3dmix4_s_mod
#if defined TS_DIF4 && defined MIX_S_TS && defined SOLVE3D
# ifdef EW_PERIODIC
#  define I_RANGE Istr-1,Iend+1
# else
#  define I_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
# endif
# ifdef NS_PERIODIC
#  define J_RANGE Jstr-1,Jend+1
# else
#  define J_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
# endif
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine computes horizontal biharmonic mixing of tracers    !
!  along S-coordinate levels surfaces.                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC t3dmix4_s

      CONTAINS
!
!***********************************************************************
      SUBROUTINE t3dmix4_s (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 27)
# endif
      CALL t3dmix4_s_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     nrhs(ng), nnew(ng),                          &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % pmon_u,                           &
     &                     GRID(ng) % pnom_v,                           &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     MIXING(ng) % diff4,                          &
     &                     OCEAN(ng) % t)
# ifdef PROFILE
      CALL wclock_off (ng, 27)
# endif
      RETURN
      END SUBROUTINE t3dmix4_s
!
!***********************************************************************
      SUBROUTINE t3dmix4_s_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           nrhs, nnew,                            &
# ifdef MASKING
     &                           umask, vmask,                          &
# endif
     &                           Hz, pmon_u, pnom_v, pm, pn,            &
     &                           diff4,                                 &
     &                           t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs, nnew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: diff4(LBi:,LBj:,:)

      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: diff4(LBi:UBi,LBj:UBj,NT(ng))

      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j, k

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: LapT

# include "set_bounds.h"
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
!
!--------------------------------------------------------------------
!  Compute horizontal biharmonic diffusion along constant S-surfaces.
!  The biharmonic operator is computed by applying the harmonic
!  operator twice.
!--------------------------------------------------------------------
!
!  Compute horizontal tracer flux in the XI- and ETA-directions.
!
          DO j=J_RANGE
            DO i=I_RANGE+1
              FX(i,j)=0.25_r8*(diff4(i,j,itrc)+diff4(i-1,j,itrc))*      &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*pmon_u(i,j)*              &
     &                (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
# ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
# endif
            END DO
          END DO
          DO j=J_RANGE+1
            DO i=I_RANGE
              FE(i,j)=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*      &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*pnom_v(i,j)*              &
     &                (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
# ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
# endif
            END DO
          END DO
!
!  Compute first harmonic operator and multiply by the metrics of the
!  second harmonic operator.
!
          DO j=J_RANGE
            DO i=I_RANGE
              LapT(i,j)=(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))*          &
     &                  pm(i,j)*pn(i,j)/Hz(i,j,k)
            END DO
          END DO
!
!  Apply boundary conditions (except periodic; closed or gradient)
!  to the first harmonic operator.
!
# ifndef EW_PERIODIC
          IF (WESTERN_EDGE) THEN
            DO j=J_RANGE
#  ifdef WESTERN_WALL
              LapT(Istr-1,j)=0.0_r8
#  else
              LapT(Istr-1,j)=LapT(Istr,j)
#  endif
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=J_RANGE
#  ifdef EASTERN_WALL
              LapT(Iend+1,j)=0.0_r8
#  else
              LapT(Iend+1,j)=LapT(Iend,j)
#  endif
            END DO
          END IF
# endif
# ifndef NS_PERIODIC
          IF (SOUTHERN_EDGE) THEN
            DO i=I_RANGE
#  ifdef SOUTHERN_WALL
              LapT(i,Jstr-1)=0.0_r8
#  else
              LapT(i,Jstr-1)=LapT(i,Jstr)
#  endif
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO i=I_RANGE
#  ifdef NORTHERN_WALL
              LapT(i,Jend+1)=0.0_r8
#  else
              LapT(i,Jend+1)=LapT(i,Jend)
#  endif
            END DO
          END IF
# endif
!
!  Compute FX=d(LapT)/d(xi) and FE=d(LapT)/d(eta) terms.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=0.25_r8*(diff4(i,j,itrc)+diff4(i-1,j,itrc))*      &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*pmon_u(i,j)*              &
     &                (LapT(i,j)-LapT(i-1,j))
# ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
# endif
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=0.25_r8*(diff4(i,j,itrc)+diff4(i,j-1,itrc))*      &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*pnom_v(i,j)*              &
     &                (LapT(i,j)-LapT(i,j-1))
# ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
# endif
            END DO
          END DO
!
! Time-step biharmonic, S-surfaces diffusion term.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-                    &
     &                           dt(ng)*pm(i,j)*pn(i,j)*                &
     &                                  (FX(i+1,j)-FX(i,j)+             &
     &                                   FE(i,j+1)-FE(i,j))
            END DO
          END DO
        END DO
      END DO
# undef I_RANGE
# undef J_RANGE
      RETURN
      END SUBROUTINE t3dmix4_s_tile
#endif
      END MODULE t3dmix4_s_mod
