#include "cppdefs.h"
      MODULE hwalk_floats_mod

#if defined NONLINEAR && defined FLOATS && defined FLOAT_HWALK && \
    !defined SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines compute nudging velocities for horizontal random     !
!  walk.                                                               !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: hwalk_floats
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE hwalk_floats (ng, Lstr, Lend, Predictor,               &
     &                         my_thread, nudgX, nudgY)
!***********************************************************************
!
      USE mod_param
      USE mod_floats
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend

      logical, intent(in) :: Predictor
# ifdef ASSUMED_SHAPE
      logical, intent(in) :: my_thread(Lstr:)

      real(r8), intent(inout) :: nudgX(Lstr:)
      real(r8), intent(inout) :: nudgY(Lstr:)
# else
      logical, intent(in) :: my_thread(Lstr:Lend)

      real(r8), intent(inout) :: nudgX(Lstr:Lend)
      real(r8), intent(inout) :: nudgY(Lstr:Lend)
# endif

!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif
      CALL hwalk_floats_tile (ng, Lstr, Lend, knew(ng),                 &
     &                        nfm3(ng), nfm2(ng), nfm1(ng), nf(ng),     &
     &                        nfp1(ng),                                 &
     &                        Predictor, my_thread,                     &
     &                        DRIFTER(ng) % bounded,                    &
     &                        DRIFTER(ng) % Tinfo,                      &
     &                        DRIFTER(ng) % rwalkX,                     &
     &                        DRIFTER(ng) % rwalkY,                     &
     &                        nudgX, nudgY,                             &
     &                        DRIFTER(ng) % track)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif

      RETURN
      END SUBROUTINE hwalk_floats

!
!***********************************************************************
      SUBROUTINE hwalk_floats_tile (ng, Lstr, Lend, knew,               &
     &                              nfm3, nfm2, nfm1, nf, nfp1,         &
     &                              Predictor, my_thread, bounded,      &
     &                              Tinfo, rwalkX, rwalkY,              &
     &                              nudgX, nudgY, track)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
      USE interp_floats_mod
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf
# endif
      USE nrutil, ONLY : gasdev
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend, knew
      integer, intent(in) :: nfm3, nfm2, nfm1, nf, nfp1
      logical, intent(in) :: Predictor
!
# ifdef ASSUMED_SHAPE
      logical, intent(in) :: bounded(:)
      logical, intent(in) :: my_thread(Lstr:)

      real(r8), intent(in) :: Tinfo(0:,:)

      real(r8), intent(inout) :: rwalkX(:)
      real(r8), intent(inout) :: rwalkY(:)
      real(r8), intent(inout) :: nudgX(Lstr:)
      real(r8), intent(inout) :: nudgY(Lstr:)
      real(r8), intent(inout) :: track(:,0:,:)
# else
      logical, intent(in) :: bounded(Nfloats(ng))
      logical, intent(in) :: my_thread(Lstr:Lend)

      real(r8), intent(in) :: Tinfo(0:izrhs,Nfloats(ng))

      real(r8), intent(inout) :: rwalkX(Nfloats(ng))
      real(r8), intent(inout) :: rwalkY(Nfloats(ng))
      real(r8), intent(inout) :: nudgX(Lstr:Lend)
      real(r8), intent(inout) :: nudgY(Lstr:Lend)
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
# endif
!
!  Local variable declarations.
!
# ifdef MASKING
      logical, parameter :: Lmask = .TRUE.
# else
      logical, parameter :: Lmask = .FALSE.
# endif
      integer :: LBi, UBi, LBj, UBj
      integer :: i, l, nfindx
      integer :: ierr

      real(r8) :: HalfDT, xrhs, yrhs
      real(r8) :: cff, cff1, cff2, cff3, cff4
!
! Set tile array bounds.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Compute nudging horizontal velocities for random walk.
!-----------------------------------------------------------------------
!
!  Set float time level index to process.
!
      IF (Predictor) THEN
        nfindx=nf
      ELSE
        nfindx=nfp1
      END IF
!
!  If predictor step, generate random number sequence.
!
      IF (Predictor) THEN
# ifdef DISTRIBUTE
        IF (Master) THEN
          CALL gasdev (rwalkX)
          CALL gasdev (rwalkY)
        END IF
        CALL mp_bcastf (ng, iNLM, rwalkX)
        CALL mp_bcastf (ng, iNLM, rwalkY)
# else
!$OMP MASTER
        CALL gasdev (rwalkX)
        CALL gasdev (rwalkY)
!$OMP END MASTER
!$OMP BARRIER
# endif
      END IF
!
!  Compute nudging velocity coefficients. Use normally distributed
!  random numbers. The diffusivity is currently specified via the
!  preprocessor macro FLOAT_HWALK_DIFF, which should be equal to
!  K/(dx*dy).
!
      cff=2.0_r8/dt(ng)
      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
          nudgX(l)=SQRT(cff*FLOAT_HWALK_DIFF)*rwalkX(l)
          nudgY(l)=SQRT(cff*FLOAT_HWALK_DIFF)*rwalkY(l)
        ELSE
          nudgX(l)=0.0_r8
          nudgY(l)=0.0_r8
        END IF
      END DO
!
!  Interpolate horizontal slopes using nudging velocity coefficients.
!
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfindx, ixrhs, isUbar,            &
     &                    -u2dvar, Lmask, spval, nudgX,                 &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % ubar(:,:,knew),                   &
     &                    my_thread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfindx, iyrhs, isVbar,            &
     &                    -v2dvar, Lmask, spval, nudgY,                 &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % vbar(:,:,knew),                   &
     &                    my_thread, bounded, track)
!
!  If newly relased float, initialize all time levels.
!
      HalfDT=0.5_r8*dt(ng)

      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
          IF (time(ng)-HalfDT.le.Tinfo(itstr,l).and.                    &
     &        time(ng)+HalfDT.gt.Tinfo(itstr,l)) THEN
            xrhs=track(ixrhs,nfindx,l)
            yrhs=track(iyrhs,nfindx,l)
            DO i=0,NFT
              track(ixrhs,i,l)=xrhs
              track(iyrhs,i,l)=yrhs
            END DO
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Time step for horizontal position.
!-----------------------------------------------------------------------
!
!  Assign predictor/corrector weights.
!
      IF (Predictor) THEN
        cff1=8.0_r8/3.0_r8
        cff2=4.0_r8/3.0_r8
      ELSE
        cff1=9.0_r8/8.0_r8
        cff2=1.0_r8/8.0_r8
        cff3=3.0_r8/8.0_r8
        cff4=6.0_r8/8.0_r8
      END IF
!
!  Compute new float horizontal position.
!
# ifdef HWALK_FORWARD
      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
          track(ixgrd,nfp1,l)=track(ixgrd,nf,l)+                        &
     &                        dt(ng)*track(ixrhs,nf,l)
          track(iygrd,nfp1,l)=track(iygrd,nf,l)+                        &
     &                        dt(ng)*track(iyrhs,nf,l)
        END IF
      END DO
# else
      IF (Predictor) THEN
        DO l=Lstr,Lend
          IF (my_thread(l).and.bounded(l)) THEN
            track(ixgrd,nfp1,l)=track(ixgrd,nfm3,l)+                    &
     &                          dt(ng)*(cff1*track(ixrhs,nf  ,l)-       &
     &                                  cff2*track(ixrhs,nfm1,l)+       &
     &                                  cff1*track(ixrhs,nfm2,l))
            track(iygrd,nfp1,l)=track(iygrd,nfm3,l)+                    &
     &                          dt(ng)*(cff1*track(iyrhs,nf  ,l)-       &
     &                                  cff2*track(iyrhs,nfm1,l)+       &
     &                                  cff1*track(iyrhs,nfm2,l))
          END IF
        END DO
      ELSE
        DO l=Lstr,Lend
          IF (my_thread(l).and.bounded(l)) THEN
            track(ixgrd,nfp1,l)=cff1*track(ixgrd,nf  ,l)-               &
     &                          cff2*track(ixgrd,nfm2,l)+               &
     &                          dt(ng)*(cff3*track(ixrhs,nfp1,l)+       &
     &                                  cff4*track(ixrhs,nf  ,l)-       &
     &                                  cff3*track(ixrhs,nfm1,l))
            track(iygrd,nfp1,l)=cff1*track(iygrd,nf  ,l)-               &
     &                          cff2*track(iygrd,nfm2,l)+               &
     &                          dt(ng)*(cff3*track(iyrhs,nfp1,l)+       &
     &                                  cff4*track(iyrhs,nf  ,l)-       &
     &                                  cff3*track(iyrhs,nfm1,l))
          END IF
        END DO
      END IF
# endif
!
!  Zeroth-out nudging velocities coefficients.
!
      DO l=Lstr,Lend
        nudgX(l)=0.0_r8
        nudgY(l)=0.0_r8
      END DO

      RETURN
      END SUBROUTINE hwalk_floats_tile
#endif
      END MODULE hwalk_floats_mod

