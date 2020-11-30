! *** Module: gaufre ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Meta-module to make the use of GAUFRE easier
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 04.2018 Created [Yann Pouillon]
! *****************************************************************************
module gaufre

  use gaufre_common
  use gaufre_config
  use gaufre_data
  use gaufre_driver
  use gaufre_io
  use gaufre_orbital

  implicit none

  public

end module gaufre
