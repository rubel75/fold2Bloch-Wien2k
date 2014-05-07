! -*- mode: f90 -*-
#ifdef HAVE_VARLEN_STR
#define argstr character(len=:), allocatable
#else
#define fetcharg fetcharg_buf
#define argstr character(len=bufsz)
#endif
