! This file holds all the different modules for sparse vectors
! This makes editing easier as we do not need to consider *all* files

! A logical sparse array
!!$module class_lSpData3D
!!$  use class_lData3D
!!$!========================
!!$#define TYPE_NAME  lSpData3D
!!$#define STR_TYPE_NAME "lSpData3D"
!!$#define TYPE_NAME_ lSpData3D_
!!$#define NEW_TYPE newlSpData3D
!!$#define VAR_TYPE lData3D
!!$#define VAR_NEW_TYPE newlData3D
!!$#define VAR_TYPE_TYPE logical
!!$! Do not define precision
!!$#include "class_SpData3D.T90"
!!$!========================
!!$end module class_lSpData3D

module class_iSpData3D
  use class_iData3D
!========================
#define TYPE_NAME  iSpData3D
#define STR_TYPE_NAME "iSpData3D"
#define TYPE_NAME_ iSpData3D_
#define NEW_TYPE newiSpData3D
#define VAR_TYPE iData3D
#define VAR_NEW_TYPE newiData3D
#define VAR_TYPE_TYPE integer
! Do not define precision
#include "class_SpData3D.T90"
!========================
end module class_iSpData3D

!!$module class_sSpData3D
!!$  use class_sData3D
!!$!========================
!!$#define TYPE_NAME  sSpData3D
!!$#define STR_TYPE_NAME  "sSpData3D"
!!$#define TYPE_NAME_ sSpData3D_
!!$#define NEW_TYPE newsSpData3D
!!$#define VAR_TYPE sData3D
!!$#define VAR_NEW_TYPE newsData3D
!!$#define VAR_TYPE_TYPE real
!!$#define PREC sp
!!$#include "class_SpData3D.T90"
!!$!========================
!!$end module class_sSpData3D

module class_dSpData3D
  use class_dData3D
!========================
#define TYPE_NAME  dSpData3D
#define STR_TYPE_NAME  "dSpData3D"
#define TYPE_NAME_ dSpData3D_
#define NEW_TYPE newdSpData3D
#define VAR_TYPE dData3D
#define VAR_NEW_TYPE newdData3D
#define VAR_TYPE_TYPE real
#define PREC dp
#include "class_SpData3D.T90"
!========================
end module class_dSpData3D

!!$module class_cSpData3D
!!$  use class_cData3D
!!$!========================
!!$#define TYPE_NAME  cSpData3D
!!$#define STR_TYPE_NAME  "cSpData3D"
!!$#define TYPE_NAME_ cSpData3D_
!!$#define NEW_TYPE newcSpData3D
!!$#define VAR_TYPE cData3D
!!$#define VAR_NEW_TYPE newcData3D
!!$#define VAR_TYPE_TYPE complex
!!$#define PREC sp
!!$#include "class_SpData3D.T90"
!!$!========================
!!$end module class_cSpData3D

module class_zSpData3D
  use class_zData3D
!========================
#define TYPE_NAME  zSpData3D
#define STR_TYPE_NAME  "zSpData3D"
#define TYPE_NAME_ zSpData3D_
#define NEW_TYPE newzSpData3D
#define VAR_TYPE zData3D
#define VAR_NEW_TYPE newzData3D
#define VAR_TYPE_TYPE complex
#define PREC dp
#include "class_SpData3D.T90"
!========================
end module class_zSpData3D

