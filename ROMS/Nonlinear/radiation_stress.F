#include "cppdefs.h"
#if defined NONLINEAR && defined NEARSHORE
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!================================================== John C. Warner ====!
!                                                                      !
!  This routine computes the radiation stress contributions to the     !
!  momentum equations.                                                 !
!                                                                      !
!=======================================================================
!
# if defined NEARSHORE_MELLOR05
#  include <nearshore_mellor05.h>
# elif defined NEARSHORE_MELLOR08
#  include <nearshore_mellor08.h>
# endif
#else
      MODULE radiation_stress_mod
      END MODULE radiation_stress_mod
#endif
