# svn $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2022 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Programs"

list( APPEND _files
      ROMS/Programs/types.F
)

set ( ROMS_Programs_files
      ${_files}
      PARENT_SCOPE
)
