MODULE ANISOFLOW_Interface

! ANISOFLOW_Interface it's a module that contains routines that serve as a bridge between user and
! program.

    IMPLICIT NONE

CONTAINS

 !  - GetInputDir: It's a routine that provide an input directory to the program given by the user.
 !    > OUT: InputDir, ierr.
 !      + InputDir: It's a string that specifies the directory of the input directory using Unix file 
 !                  system.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user may provide a direction using "-Input_dir" followed by the path in Unix
 !             file system.      

SUBROUTINE GetInputDir(InputDir,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputDir

    PetscBool                       :: InputDirFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_dir",InputDir,     &
        & InputDirFlg,ierr)
    IF (.NOT.InputDirFlg) InputDir=""
    InputDir=TRIM(InputDir)

END SUBROUTINE GetInputDir

 !  - GetInputType: It's a routine that provides a set of file types to use as a input file.
 !    > OUT: InputType, ierr.
 !      + InputType: It's a collection of integers that define the type of input to be used in the
 !                   program.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user may provide a set of file types to be used using "-Input_type" followed by
 !             an integer that defines the set:
 !                  1: Set of files types where the identifier will be Gmtry=1, Tplgy=1, Cvt=1,
 !                     and BC=1. See InputType definition for more information.

SUBROUTINE GetInputType(InputType,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    Type(InputTypeVar),INTENT(OUT)  :: InputType

    PetscBool                       :: InputTypeFlg
    PetscInt                        :: InputTypeTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_type",InputTypeTmp,   &
        & InputTypeFlg,ierr)

    IF (InputTypeFlg) THEN
        IF (InputTypeTmp.EQ.1) THEN
            InputType%Gmtry=1
            InputType%Tplgy=1
            InputType%Cvt=1
            InputType%Sto=1
            InputType%BC=1
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Input_type used is invalid\n",ierr)
            STOP
        END IF
    ELSE
        ! Presseting default InputType
        InputType%Gmtry=1
        InputType%Tplgy=1
        InputType%Cvt=1
        InputType%Sto=1
        InputType%BC=1
    END IF
    
    ! Setting InputType from interface
    CALL GetInputTypeGmtry(InputType,ierr)
    CALL GetInputTypeTplgy(InputType,ierr)
    CALL GetInputTypeCvt(InputType,ierr)
    CALL GetInputTypeSto(InputType,ierr)
    CALL GetInputTypeBC(InputType,ierr)

END SUBROUTINE GetInputType

 !  - GetInputTypeGmtry: It's a routine that provide a file type to use as input geometry.
 !    > OUT: InputType, ierr.
 !      + InputType: It's a collection of integers that define the type of input to be used in the
 !                   program.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user may provide a set of file types to be used using "-Input_type" followed by
 !             an integer that defines the set:
 !                  1: Set of files types where the identifier will be Gmtry=1, Tplgy=1, Cvt=1,
 !                     and BC=1. See InputType definition for more information.
 !             or use "-Input_type_gmtry" followed by an integer that the file type to use as input 
 !             geometry:
 !                  1 and 2: Defined by Perez. An example is provided in "../ex/Perez/in/sanpck.domnRST"
 !                  1 has a default grid, 2 porvided by the file but both has the same syntax.

SUBROUTINE GetInputTypeGmtry(InputType,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeGmtryFlg
    PetscInt                        :: InputTypeGmtryTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_type_gmtry",          &
        & InputTypeGmtryTmp,InputTypeGmtryFlg,ierr)

    IF (InputTypeGmtryFlg) THEN
        IF (InputTypeGmtryTmp.EQ.1) THEN
            InputType%Gmtry=InputTypeGmtryTmp
        ELSEIF (InputTypeGmtryTmp.EQ.2) THEN
            InputType%Gmtry=InputTypeGmtryTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Input_type_gmtry used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeGmtry

 !  - GetInputTypeTplgy: It's a routine that provides a file type to use as input topology.
 !    > OUT: InputType, ierr.
 !      + InputType: It's a collection of integers that define the type of input to be used in the
 !                   program.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user may provide a set of file types to be used using "-Input_type" followed by
 !             an integer that defines the set:
 !                  1: Set of files types where the identifier will be Gmtry=1, Tplgy=1, Cvt=1,
 !                     and BC=1. See InputType definition for more information.
 !             or use "-Input_type_tplgy" followed by an integer that the file type to use as input 
 !             topology:
 !                  1: Default topology, it doesn't need a file. Every face of the domain is a 
 !                     Neumman condition but on the first layer where the boundary is dirichlet.

SUBROUTINE GetInputTypeTplgy(InputType,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeTplgyFlg
    PetscInt                        :: InputTypeTplgyTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_type_tplgy",          &
        & InputTypeTplgyTmp,InputTypeTplgyFlg,ierr)

    IF (InputTypeTplgyFlg) THEN
        IF (InputTypeTplgyTmp.EQ.1) THEN
            InputType%Tplgy=InputTypeTplgyTmp
        ELSE IF (InputTypeTplgyTmp.EQ.2) THEN
            InputType%Tplgy=InputTypeTplgyTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Input_type_tplgy used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeTplgy

 !  - GetInputTypeCvt: It's a routine that provides a file type to use as input conductivity.
 !    > OUT: InputType, ierr.
 !      + InputType: It's a collection of integers that define the type of input to be used in the
 !                   program.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user may provide a set of file types to be used using "-Input_type" followed by
 !             an integer that defines the set:
 !                  1: Set of files types where the identifier will be Gmtry=1, Tplgy=1, Cvt=1,
 !                     and BC=1. See InputType definition for more information.
 !             or use "-Input_type_cvt" followed by an integer that the file type to use as input 
 !             conductivity:
 !                  1: It's a pair of files that provide block conductivities by zones and another 
 !                     one that provides an zone identifier to each cell. An example is 
 !                     provided in "../ex/Blessent/in/matrix.mprops" and 
 !                     "../ex/Blessent/in/tsim_USMH.asc".

SUBROUTINE GetInputTypeCvt(InputType,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeCvtFlg
    PetscInt                        :: InputTypeCvtTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_type_cvt",            &
        & InputTypeCvtTmp,InputTypeCvtFlg,ierr)

    IF (InputTypeCvtFlg) THEN
        IF (InputTypeCvtTmp.EQ.1) THEN
            InputType%Cvt=InputTypeCvtTmp
        ELSE IF (InputTypeCvtTmp.EQ.2) THEN
            InputType%Cvt=InputTypeCvtTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Input_type_cvt used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeCvt

 !  - GetInputTypeSto: It's a routine that provides a file type to use as input specific storage.
 !    > OUT: InputType, ierr.
 !      + InputType: It's a collection of integers that define the type of input to be used in the
 !                   program.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user may provide a set of file types to be used using "-Input_type" followed by
 !             an integer that defines the set:
 !                  1: Set of files types where the identifier will be Gmtry=1, Tplgy=1, Cvt=1, Sto=1,
 !                     and BC=1. See InputType definition for more information.
 !             or use "-Input_type_Sto" followed by an integer that the file type to use as input 
 !             specific storage:
 !                  1: It's a pair of files that provide block conductivities by zones and another 
 !                     one that provides an zone identifier to each cell. An example is 
 !                     provided in "../ex/Blessent/in/matrix.mprops" and 
 !                     "../ex/Blessent/in/tsim_USMH.asc".

SUBROUTINE GetInputTypeSto(InputType,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeStoFlg
    PetscInt                        :: InputTypeStoTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_type_sto",            &
        & InputTypeStoTmp,InputTypeStoFlg,ierr)

    IF (InputTypeStoFlg) THEN
        IF (InputTypeStoTmp.EQ.1) THEN
            InputType%Sto=InputTypeStoTmp
        ELSE IF (InputTypeStoTmp.EQ.2) THEN
            InputType%Sto=InputTypeStoTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Input_type_sto used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeSto

 !  - GetInputTypeBC: It's a routine that provides a file type to use as input boundary condition.
 !    > OUT: InputType, ierr.
 !      + InputType: It's a collection of integers that define the type of input to be used in the
 !                   program.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user may provide a set of file types to be used using "-Input_type" followed by
 !             an integer that defines the set:
 !                  1: Set of files types where the identifier will be Gmtry=1, Tplgy=1, Cvt=1,
 !                     and BC=1. See InputType definition for more information.
 !             or use "-Input_type_bc" followed by an integer that the file type to use as input 
 !             boundary condition:
 !                  1: It's a file that only provide dirichlet condition and their postion on the 
 !                     grid. An example is provided in 
 !                     "../Blessent/in/grid_400_400.nch_nprop_list.lateral_boundary".

SUBROUTINE GetInputTypeBC(InputType,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(InputTypeVar),INTENT(INOUT)    :: InputType

    PetscBool                       :: InputTypeBCFlg
    PetscInt                        :: InputTypeBCTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_type_bc",      &
        & InputTypeBCTmp,InputTypeBCFlg,ierr)

    IF (InputTypeBCFlg) THEN
        IF (InputTypeBCTmp.EQ.1) THEN
            InputType%BC=InputTypeBCTmp
        ELSEIF (InputTypeBCTmp.EQ.2) THEN
            InputType%BC=InputTypeBCTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Input_type_bc used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeBC

 !  - GetInputFileGmtry: It's a routine that provides a file name to open the geometry file in 
 !                       InputDir.
 !    > OUT: InputFileGmtry, ierr.
 !      + InputFileGmtry: It's a string that specifies file name to open the geometry file in InputDir.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user must provide a file name using "-Input_file_gmtry" followed by geometry file 
 !             name.


SUBROUTINE GetInputFileGmtry(InputFileGmtry,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileGmtry

    PetscBool                       :: InputFileGmtryFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_file_gmtry",       &
        InputFileGmtry,InputFileGmtryFlg,ierr)

    IF (.NOT.InputFileGmtryFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Input_file_gmtry command must be used\n",ierr)
        STOP
    END IF

    InputFileGmtry=TRIM(InputFileGmtry)

END SUBROUTINE GetInputFileGmtry

 !  - GetInputFileTplgy: It's a routine that provides a file name to open the topology file in 
 !                       InputDir.
 !    > OUT: InputFileTplgy, ierr.
 !      + InputFileTplgy: It's a string that specifies file name to open the topology file in InputDir.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user may provide a file name using "-Input_file_tplgy" followed by topology file 
 !             name.


SUBROUTINE GetInputFileTplgy(InputFileTplgy,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileTplgy

    PetscBool                       :: InputFileTplgyFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_file_tplgy",       &
        InputFileTplgy,InputFileTplgyFlg,ierr)

    IF (.NOT.InputFileTplgyFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Input_file_tplgy command must be used\n",ierr)
        STOP
    END IF

    InputFileTplgy=TRIM(InputFileTplgy)

END SUBROUTINE GetInputFileTplgy

 !  - GetInputFilePptByZone: It's a routine that provides a file name to open the property  
 !                           zones file in InputDir. 
 !    > OUT: InputFilePptByZone, InputFilePptByZoneFlg, ierr.
 !      + InputFilePptByZone: It's a string that specifies file name to open the property zones 
 !                             file in InputDir.
 !      + InputFilePptByZoneFlg: ...
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user can provide a file name using "-Input_file_ppt_by_zones" followed by 
 !             property zones file name.

SUBROUTINE GetInputFilePptByZone(InputFilePptByZone,InputFilePptByZoneFlg,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFilePptByZone
    PetscBool,INTENT(OUT)           :: InputFilePptByZoneFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_file_ppt_by_zones",    &
        InputFilePptByZone,InputFilePptByZoneFlg,ierr)

    InputFilePptByZone=TRIM(InputFilePptByZone)

END SUBROUTINE GetInputFilePptByZone


 !  - GetInputFileCvt: It's a routine that provides a file name to open the conductivity file in 
 !                       InputDir.
 !    > OUT: InputFileCvt, ierr.
 !      + InputFileCvt: It's a string that specifies file name to open the conductivity file in 
 !                      InputDir.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user must provide a file name using "-Input_file_cvt" followed by conductiviy
 !             file name.

SUBROUTINE GetInputFileCvt(InputFileCvt,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileCvt

    PetscBool                       :: InputFileCvtFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_file_cvt",         &
        InputFileCvt,InputFileCvtFlg,ierr)

    IF (.NOT.InputFileCvtFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Input_file_cvt command must be used\n",ierr)
        STOP
    END IF

    InputFileCvt=TRIM(InputFileCvt)

END SUBROUTINE GetInputFileCvt

 !  - GetInputFileCvtByZone: It's a routine that provides a file name to open the property  
 !                           zones file in InputDir. 
 !    > OUT: InputFileCvtByZone, InputFileCvtByZoneFlg, ierr.
 !      + InputFileCvtByZone: It's a string that specifies file name to open the property zones 
 !                             file in InputDir.
 !      + InputFileCvtByZoneFlg: ...
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user can provide a file name using "-Input_file_cvt_by_zones" followed by 
 !             property zones file name.

SUBROUTINE GetInputFileCvtByZone(InputFileCvtByZone,InputFileCvtByZoneFlg,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileCvtByZone
    PetscBool,INTENT(OUT)           :: InputFileCvtByZoneFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_file_cvt_by_zones",    &
        InputFileCvtByZone,InputFileCvtByZoneFlg,ierr)

    InputFileCvtByZone=TRIM(InputFileCvtByZone)

END SUBROUTINE GetInputFileCvtByZone

 !  - GetInputFileSto: It's a routine that provides a file name to open the Specific Storage file in 
 !                       InputDir.
 !    > OUT: InputFileSto, ierr.
 !      + InputFileSto: It's a string that specifies file name to open the Specific Storage file in 
 !                      InputDir.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user must provide a file name using "-Input_file_sto" followed by Specific Storage
 !             file name.

SUBROUTINE GetInputFileSto(InputFileSto,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileSto

    PetscBool                       :: InputFileStoFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_file_sto",         &
        InputFileSto,InputFileStoFlg,ierr)

    IF (.NOT.InputFileStoFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Input_file_sto command must be used\n",ierr)
        STOP
    END IF

    InputFileSto=TRIM(InputFileSto)

END SUBROUTINE GetInputFileSto

 !  - GetInputFileStoByZone: It's a routine that provides a file name to open the property  
 !                           zones file in InputDir. 
 !    > OUT: InputFileStoByZone, InputFileStoByZoneFlg, ierr.
 !      + InputFileStoByZone: It's a string that specifies file name to open the property zones 
 !                             file in InputDir.
 !      + InputFileStoByZoneFlg: ...
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user can provide a file name using "-Input_file_sto_by_zones" followed by 
 !             property zones file name.

SUBROUTINE GetInputFileStoByZone(InputFileStoByZone,InputFileStoByZoneFlg,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileStoByZone
    PetscBool,INTENT(OUT)           :: InputFileStoByZoneFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_file_sto_by_zones",    &
        InputFileStoByZone,InputFileStoByZoneFlg,ierr)

    InputFileStoByZone=TRIM(InputFileStoByZone)

END SUBROUTINE GetInputFileStoByZone

 !  - GetInputFileBC: It's a routine that provides a file name to open the boundary condition file
 !                    in InputDir. 
 !    > OUT: InputFileBC, ierr.
 !      + InputFileBC: It's a string that specifies file name to open the boundary condition file in
 !                     InputDir.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user must provide a file name using "-Input_file_bc" followed by boundary
 !             condition file name.

SUBROUTINE GetInputFileBC(InputFileBC,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileBC

    PetscBool                       :: InputFileBCFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Input_file_bc",   &
        InputFileBC,InputFileBCFlg,ierr)

    IF (.NOT.InputFileBCFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Input_file_bc command must be used.\n",ierr)
        STOP
    END IF

    InputFileBC=TRIM(InputFileBC)

END SUBROUTINE GetInputFileBC

 !    

SUBROUTINE GetOutputDir(OutputDir,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: OutputDir

    PetscBool                       :: OutputDirFlg

    CALL PetscOptionsGetString(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Output_dir",OutputDir,     &
        & OutputDirFlg,ierr)
    IF (.NOT.OutputDirFlg) OutputDir=""
    OutputDir=TRIM(OutputDir)

END SUBROUTINE GetOutputDir

 !

SUBROUTINE GetOutputType(OutputType,ierr)

    USE ANISOFLOW_Types, ONLY : OutputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    Type(OutputTypeVar),INTENT(OUT)  :: OutputType

    PetscBool                       :: OutputTypeFlg
    PetscInt                        :: OutputTypeTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Output_type",OutputTypeTmp,   &
        & OutputTypeFlg,ierr)

    IF (OutputTypeFlg) THEN
        IF (OutputTypeTmp.EQ.0) THEN
            OutputType%Tplgy=0
            OutputType%Cvt=0
            OutputType%Sol=0
        ELSE IF (OutputTypeTmp.EQ.1) THEN
            OutputType%Tplgy=1
            OutputType%Cvt=1
            OutputType%Sol=1
        ELSE IF (OutputTypeTmp.EQ.2) THEN
            OutputType%Tplgy=2
            OutputType%Cvt=2
            OutputType%Sol=2
        ELSE IF (OutputTypeTmp.EQ.3) THEN
            OutputType%Tplgy=3
            OutputType%Cvt=3
            OutputType%Sol=3
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Output_type used is invalid\n",ierr)
            STOP
        END IF
    ELSE
        ! Presseting default OutputType
        OutputType%Tplgy=0
        OutputType%Cvt=0
        OutputType%Sol=1
    END IF
    
    ! Setting OutputType from interface
    CALL GetOutputTypeTplgy(OutputType,ierr)
    CALL GetOutputTypeCvt(OutputType,ierr)
    CALL GetOutputTypeSol(OutputType,ierr)

END SUBROUTINE GetOutputType

SUBROUTINE GetOutputTypeTplgy(OutputType,ierr)

    USE ANISOFLOW_Types, ONLY : OutputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(OutputTypeVar),INTENT(INOUT)    :: OutputType

    PetscBool                       :: OutputTypeTplgyFlg
    PetscInt                        :: OutputTypeTplgyTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Output_type_tplgy",          &
        & OutputTypeTplgyTmp,OutputTypeTplgyFlg,ierr)

    IF (OutputTypeTplgyFlg) THEN
        IF (OutputTypeTplgyTmp.EQ.0) THEN
            OutputType%Tplgy=OutputTypeTplgyTmp
        ELSE IF (OutputTypeTplgyTmp.EQ.1) THEN
            OutputType%Tplgy=OutputTypeTplgyTmp
        ELSE IF (OutputTypeTplgyTmp.EQ.2) THEN
            OutputType%Tplgy=OutputTypeTplgyTmp
        ELSE IF (OutputTypeTplgyTmp.EQ.3) THEN
            OutputType%Tplgy=OutputTypeTplgyTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Output_type_tplgy used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetOutputTypeTplgy

SUBROUTINE GetOutputTypeCvt(OutputType,ierr)

    USE ANISOFLOW_Types, ONLY : OutputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(OutputTypeVar),INTENT(INOUT)    :: OutputType

    PetscBool                       :: OutputTypeCvtFlg
    PetscInt                        :: OutputTypeCvtTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Output_type_cvt",          &
        & OutputTypeCvtTmp,OutputTypeCvtFlg,ierr)

    IF (OutputTypeCvtFlg) THEN
        IF (OutputTypeCvtTmp.EQ.0) THEN
            OutputType%Cvt=OutputTypeCvtTmp
        ELSE IF (OutputTypeCvtTmp.EQ.1) THEN
            OutputType%Cvt=OutputTypeCvtTmp
        ELSE IF (OutputTypeCvtTmp.EQ.2) THEN
            OutputType%Cvt=OutputTypeCvtTmp
        ELSE IF (OutputTypeCvtTmp.EQ.3) THEN
            OutputType%Cvt=OutputTypeCvtTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Output_type_cvt used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetOutputTypeCvt

SUBROUTINE GetOutputTypeSol(OutputType,ierr)

    USE ANISOFLOW_Types, ONLY : OutputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(OutputTypeVar),INTENT(INOUT)    :: OutputType

    PetscBool                       :: OutputTypeSolFlg
    PetscInt                        :: OutputTypeSolTmp

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Output_type_sol",          &
        & OutputTypeSolTmp,OutputTypeSolFlg,ierr)

    IF (OutputTypeSolFlg) THEN
        IF (OutputTypeSolTmp.EQ.0) THEN
            OutputType%Sol=OutputTypeSolTmp
        ELSE IF (OutputTypeSolTmp.EQ.1) THEN
            OutputType%Sol=OutputTypeSolTmp
        ELSE IF (OutputTypeSolTmp.EQ.2) THEN
            OutputType%Sol=OutputTypeSolTmp
        ELSE IF (OutputTypeSolTmp.EQ.3) THEN
            OutputType%Sol=OutputTypeSolTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Output_type_sol used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetOutputTypeSol

 !  - GetRunOptions: It's a routine that provides RunOptionsVar structure that contains all options 
 !                   related with the running.
 !    > OUT: RunOptions, ierr.
 !      + RunOptions: It's a collection of integers that that contains all options related with the
 !                    running.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user may provide the running options using "-Run_options_scheme" to set a scheme
 !             stencil and "-Run_options_time" to set the time option followed by an integer to be 
 !             used. See RunOptionsVar for more information.

SUBROUTINE GetRunOptions(RunOptions,ierr)

    USE ANISOFLOW_Types, ONLY : RunOptionsVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(RunOptionsVar),INTENT(OUT) :: RunOptions

    PetscBool                       :: RunOptionsTimeFlg,RunOptionsSchemeFlg

    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Run_options_scheme",        &
        RunOptions%Scheme,RunOptionsSchemeFlg,ierr)
    CALL PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Run_options_time",         &
        RunOptions%Time,RunOptionsTimeFlg,ierr)

    IF (RunOptionsSchemeFlg) THEN
        IF ((RunOptions%Scheme.GE.3).OR.(RunOptions%Scheme.LE.0)) THEN
            CALL PetscSynchronizedPrintf(PETSC_NULL_OBJECT,PETSC_COMM_WORLD,                     &
                & "[ERROR] Run_options_scheme command must be an integer between 1 and 2\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetRunOptions

SUBROUTINE GetVerbose(Verbose,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    PetscBool,INTENT(OUT)           :: Verbose

    PetscBool                       :: VerboseFlg

    CALL PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,"-Verbose",         &
        Verbose,VerboseFlg,ierr)

    IF (.NOT.VerboseFlg) THEN
        Verbose=.TRUE.
    END IF

END SUBROUTINE GetVerbose

END MODULE ANISOFLOW_Interface