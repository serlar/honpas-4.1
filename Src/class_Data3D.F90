! We define all variables in this one
! As the entries are so short this makes more sense
!!$module class_lData3D
!!$!========================
!!$#define TYPE_NAME lData3D
!!$#define STR_TYPE_NAME "lData3D"
!!$#define TYPE_NAME_ lData3D_
!!$#define NEW_TYPE newlData3D
!!$#define VAR_TYPE logical
!!$#define VAR_INIT .false.
!!$#include "class_Data3D.T90"
!!$!========================
!!$end module class_lData3D

module class_iData3D
!========================
#define TYPE_NAME iData3D
#define STR_TYPE_NAME "iData3D"
#define TYPE_NAME_ iData3D_
#define NEW_TYPE newiData3D
#define VAR_TYPE integer
#define VAR_INIT 0
#include "class_Data3D.T90"
!========================
end module class_iData3D


!!$module class_sData3D
!!$!========================
!!$#define TYPE_NAME sData3D
!!$#define STR_TYPE_NAME "sData3D"
!!$#define TYPE_NAME_ sData3D_
!!$#define NEW_TYPE newsData3D
!!$#define VAR_TYPE real
!!$#define PREC sp
!!$#define VAR_INIT 0._sp
!!$#include "class_Data3D.T90"
!!$!========================
!!$end module class_sData3D

module class_dData3D
!========================
#define TYPE_NAME dData3D
#define STR_TYPE_NAME "dData3D"
#define TYPE_NAME_ dData3D_
#define NEW_TYPE newdData3D
#define VAR_TYPE real
#define PREC dp
#define VAR_INIT 0._dp
#include "class_Data3D.T90"
!========================
end module class_dData3D

!!$module class_cData3D
!!$!========================
!!$#define TYPE_NAME cData3D
!!$#define STR_TYPE_NAME "cData3D"
!!$#define TYPE_NAME_ cData3D_
!!$#define NEW_TYPE newcData3D
!!$#define VAR_TYPE complex
!!$#define PREC sp
!!$#define VAR_INIT cmplx(0._sp,0._sp)
!!$#include "class_Data3D.T90"
!!$!========================
!!$end module class_cData3D

module class_zData3D
!========================
#define TYPE_NAME zData3D
#define STR_TYPE_NAME "zData3D"
#define TYPE_NAME_ zData3D_
#define NEW_TYPE newzData3D
#define VAR_TYPE complex
#define PREC dp
#define VAR_INIT cmplx(0._dp,0._dp,dp)
#include "class_Data3D.T90"
!========================
end module class_zData3D


