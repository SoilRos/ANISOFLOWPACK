MODULE ANISOFLOW_Interface

    USE ANISOFLOW_Types, ONLY : InputTypeVar
    USE ANISOFLOW_Types, ONLY : RunOptionsVar

    IMPLICIT NONE

CONTAINS

SUBROUTINE GetInputDir(InputDir,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputDir

    PetscBool                       :: InputDirFlg

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_dir",InputDir,     &
        & InputDirFlg,ierr)
    IF (.NOT.InputDirFlg) InputDir="/"
    InputDir=TRIM(InputDir)
END SUBROUTINE GetInputDir

SUBROUTINE GetInputType(InputType,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    Type(InputTypeVar),INTENT(OUT)  :: InputType

    PetscBool                       :: InputTypeFlg
    PetscInt                        :: InputTypeTmp

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type",InputTypeTmp,   &
        & InputTypeFlg,ierr)

    IF (InputTypeFlg) THEN
        IF (InputTypeTmp.EQ.1) THEN
            InputType%Gmtry=1
            InputType%Tplgy=1
            InputType%Cvt=1
            InputType%SteadyBC=1
            InputType%TransientBC=1
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "ERROR: Input_type used is invalid\n",ierr)
            STOP
        END IF
    ELSE
        ! Presseting default InputType
        InputType%Gmtry=1
        InputType%Tplgy=1
        InputType%Cvt=1
        InputType%SteadyBC=1
        InputType%TransientBC=1
        ! Setting InputType from interface
        CALL GetInputTypeGmtry(InputType,ierr)
        CALL GetInputTypeTplgy(InputType,ierr)
        CALL GetInputTypeCvt(InputType,ierr)
        CALL GetInputTypeSteadyBC(InputType,ierr)
        CALL GetInputTypeTransientBC(InputType,ierr)
    END IF

END SUBROUTINE GetInputType

SUBROUTINE GetInputTypeGmtry(InputType,ierr)
    IMPLICIT NONE
#include <petsc/finclude/petscsys.h>
    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeGmtryFlg
    PetscInt                        :: InputTypeGmtryTmp

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type_gmtry",          &
        & InputTypeGmtryTmp,InputTypeGmtryFlg,ierr)

    IF (InputTypeGmtryFlg) THEN
        IF (InputTypeGmtryTmp.EQ.1) THEN
            InputType%Gmtry=InputTypeGmtryTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "ERROR: Input_type_gmtry used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeGmtry

SUBROUTINE GetInputTypeTplgy(InputType,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeTplgyFlg
    PetscInt                        :: InputTypeTplgyTmp

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type_tplgy",          &
        & InputTypeTplgyTmp,InputTypeTplgyFlg,ierr)

    IF (InputTypeTplgyFlg) THEN
        IF (InputTypeTplgyTmp.EQ.1) THEN
            InputType%Tplgy=InputTypeTplgyTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "ERROR: Input_type_tplgy used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeTplgy

SUBROUTINE GetInputTypeCvt(InputType,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeCvtFlg
    PetscInt                        :: InputTypeCvtTmp

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type_cvt",            &
        & InputTypeCvtTmp,InputTypeCvtFlg,ierr)

    IF (InputTypeCvtFlg) THEN
        IF (InputTypeCvtTmp.EQ.1) THEN
            InputType%Cvt=InputTypeCvtTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "ERROR: Input_type_Cvt used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeCvt

SUBROUTINE GetInputTypeSteadyBC(InputType,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeSteadyBCFlg
    PetscInt                        :: InputTypeSteadyBCTmp

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type_bc_steady",      &
        & InputTypeSteadyBCTmp,InputTypeSteadyBCFlg,ierr)

    IF (InputTypeSteadyBCFlg) THEN
        IF (InputTypeSteadyBCTmp.EQ.1) THEN
            InputType%SteadyBC=InputTypeSteadyBCTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "ERROR: Input_type_bc_steady used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeSteadyBC

SUBROUTINE GetInputTypeTransientBC(InputType,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeTransientBCFlg
    PetscInt                        :: InputTypeTransientBCTmp

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type_bc_transient",   &
        & InputTypeTransientBCTmp,InputTypeTransientBCFlg,ierr)

    IF (InputTypeTransientBCFlg) THEN
        IF (InputTypeTransientBCTmp.EQ.1) THEN
            InputType%TransientBC=InputTypeTransientBCTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "ERROR: Input_type_bc_transient used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeTransientBC

SUBROUTINE GetInputFileGmtry(InputFileGmtry,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileGmtry

    PetscBool                       :: InputFileGmtryFlg

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_gmtry",       &
        InputFileGmtry,InputFileGmtryFlg,ierr)

    IF (.NOT.InputFileGmtryFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: Input_file_gmtry command must be used\n",ierr)
        STOP
    END IF

    InputFileGmtry=TRIM(InputFileGmtry)

END SUBROUTINE GetInputFileGmtry

SUBROUTINE GetInputFileTplgy(InputFileTplgy,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileTplgy

    PetscBool                       :: InputFileTplgyFlg

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_tplgy",       &
        InputFileTplgy,InputFileTplgyFlg,ierr)

    IF (.NOT.InputFileTplgyFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: Input_file_tplgy command must be used\n",ierr)
        STOP
    END IF

    InputFileTplgy=TRIM(InputFileTplgy)

END SUBROUTINE GetInputFileTplgy

SUBROUTINE GetInputFileCvt(InputFileCvt,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileCvt

    PetscBool                       :: InputFileCvtFlg

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_cvt",         &
        InputFileCvt,InputFileCvtFlg,ierr)

    IF (.NOT.InputFileCvtFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: Input_file_cvt command must be used\n",ierr)
        STOP
    END IF

    InputFileCvt=TRIM(InputFileCvt)

END SUBROUTINE GetInputFileCvt

SUBROUTINE GetInputFileCvtType(InputFileCvtType,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileCvtType

    PetscBool                       :: InputFileCvtTypeFlg

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_cvt_type",    &
        InputFileCvtType,InputFileCvtTypeFlg,ierr)

    IF (.NOT.InputFileCvtTypeFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: Input_file_cvt_type command must be used\n",ierr)
        STOP
    END IF

    InputFileCvtType=TRIM(InputFileCvtType)

END SUBROUTINE GetInputFileCvtType

SUBROUTINE GetInputFileSteadyBC(InputFileSteadyBC,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileSteadyBC

    PetscBool                       :: InputFileSteadyBCFlg

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_bc_steady",   &
        InputFileSteadyBC,InputFileSteadyBCFlg,ierr)

    IF (.NOT.InputFileSteadyBCFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: Input_file_bc_steady command must be used if you want a steady simulation\n",ierr)
        STOP
    END IF

    InputFileSteadyBC=TRIM(InputFileSteadyBC)

END SUBROUTINE GetInputFileSteadyBC

SUBROUTINE GetRunOptions(RunOptions,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(RunOptionsVar),INTENT(OUT) :: RunOptions

    PetscBool                       :: RunOptionsTimeFlg,RunOptionsSchemeFlg

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Run_options_scheme",        &
        RunOptions%Scheme,RunOptionsSchemeFlg,ierr)
    CALL PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-Run_options_time",         &
        RunOptions%Time,RunOptionsTimeFlg,ierr)

    IF (RunOptionsSchemeFlg) THEN
        IF ((RunOptions%Scheme.GE.3).OR.(RunOptions%Scheme.LE.0)) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "ERROR: Run_options_scheme command must be an integer between 1 and 2\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetRunOptions

END MODULE ANISOFLOW_Interface