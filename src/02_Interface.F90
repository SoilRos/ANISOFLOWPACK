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

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_dir",InputDir,     &
        & InputDirFlg,ierr)
    IF (.NOT.InputDirFlg) InputDir="/"
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

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type",InputTypeTmp,   &
        & InputTypeFlg,ierr)

    IF (InputTypeFlg) THEN
        IF (InputTypeTmp.EQ.1) THEN
            InputType%Gmtry=1
            InputType%Tplgy=1
            InputType%Cvt=1
            InputType%BC=1
        ELSE IF (InputTypeTmp.EQ.2) THEN
            InputType%Gmtry=2
            InputType%Tplgy=2
            InputType%Cvt=2
            InputType%BC=1
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Input_type used is invalid\n",ierr)
            STOP
        END IF
    END IF
    
    ! Setting InputType from interface
    CALL GetInputTypeGmtry(InputType,ierr)
    CALL GetInputTypeTplgy(InputType,ierr)
    CALL GetInputTypeCvt(InputType,ierr)
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
 !                  1: Defined by Blessent. An example is provided in "../ex/Blessent/in/tsim_USMH.asc"
 !                  2: Defined by Perez. An example is provided in "../ex/Perez/in/sanpck.domnRST"

SUBROUTINE GetInputTypeGmtry(InputType,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar

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
        ELSE IF (InputTypeGmtryTmp.EQ.2) THEN
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

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type_tplgy",          &
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

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type_cvt",            &
        & InputTypeCvtTmp,InputTypeCvtFlg,ierr)

    IF (InputTypeCvtFlg) THEN
        IF (InputTypeCvtTmp.EQ.1) THEN
            InputType%Cvt=InputTypeCvtTmp
        ELSE IF (InputTypeCvtTmp.EQ.2) THEN
            InputType%Cvt=InputTypeCvtTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Input_type_Cvt used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetInputTypeCvt

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

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Input_type_bc",      &
        & InputTypeBCTmp,InputTypeBCFlg,ierr)

    IF (InputTypeBCFlg) THEN
        IF (InputTypeBCTmp.EQ.1) THEN
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

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_gmtry",       &
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

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_tplgy",       &
        InputFileTplgy,InputFileTplgyFlg,ierr)

    IF (.NOT.InputFileTplgyFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Input_file_tplgy command must be used\n",ierr)
        STOP
    END IF

    InputFileTplgy=TRIM(InputFileTplgy)

END SUBROUTINE GetInputFileTplgy

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

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_cvt",         &
        InputFileCvt,InputFileCvtFlg,ierr)

    IF (.NOT.InputFileCvtFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Input_file_cvt command must be used\n",ierr)
        STOP
    END IF

    InputFileCvt=TRIM(InputFileCvt)

END SUBROUTINE GetInputFileCvt

 !  - GetInputFileCvtByZones: It's a routine that provides a file name to open the conductivity  
 !                            zones file in InputDir. 
 !    > OUT: InputFileCvtByZones, ierr.
 !      + InputFileCvtByZones: It's a string that specifies file name to open the conductivity zones 
 !                             file in InputDir.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The user must provide a file name using "-Input_file_cvt_by_zones" followed by 
 !             conductiviy zones file name.

SUBROUTINE GetInputFileCvtByZones(InputFileCvtByZones,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: InputFileCvtByZones

    PetscBool                       :: InputFileCvtByZonesFlg

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_cvt_by_zones",    &
        InputFileCvtByZones,InputFileCvtByZonesFlg,ierr)

    IF (.NOT.InputFileCvtByZonesFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Input_file_cvt_by_zones command must be used.\n",ierr)
        STOP
    END IF

    InputFileCvtByZones=TRIM(InputFileCvtByZones)

END SUBROUTINE GetInputFileCvtByZones

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

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Input_file_bc",   &
        InputFileBC,InputFileBCFlg,ierr)

    IF (.NOT.InputFileBCFlg) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Input_file_bc command must be used.\n",ierr)
        STOP
    END IF

    InputFileBC=TRIM(InputFileBC)

END SUBROUTINE GetInputFileBC

 !    

SUBROUTINE GetOuputDir(OuputDir,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    CHARACTER(LEN=200),INTENT(OUT)  :: OuputDir

    PetscBool                       :: OuputDirFlg

    CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,"-Ouput_dir",OuputDir,     &
        & OuputDirFlg,ierr)
    IF (.NOT.OuputDirFlg) OuputDir="/"
    OuputDir=TRIM(OuputDir)

END SUBROUTINE GetOuputDir

 !

SUBROUTINE GetOuputType(OuputType,ierr)

    USE ANISOFLOW_Types, ONLY : OuputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    Type(OuputTypeVar),INTENT(OUT)  :: OuputType

    PetscBool                       :: OuputTypeFlg
    PetscInt                        :: OuputTypeTmp

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Ouput_type",OuputTypeTmp,   &
        & OuputTypeFlg,ierr)

    IF (OuputTypeFlg) THEN
        IF (OuputTypeTmp.EQ.1) THEN
            OuputType%Sol=1
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Ouput_type used is invalid\n",ierr)
            STOP
        END IF
    ELSE
        ! Presseting default OuputType
        OuputType%Sol=1
    END IF
    
    ! Setting OuputType from interface
    CALL GetOuputTypeSol(OuputType,ierr)


END SUBROUTINE GetOuputType

SUBROUTINE GetOuputTypeSol(OuputType,ierr)

    USE ANISOFLOW_Types, ONLY : OuputTypeVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Type(OuputTypeVar),INTENT(INOUT)    :: OuputType

    PetscBool                       :: OuputTypeSolFlg
    PetscInt                        :: OuputTypeSolTmp

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Ouput_type_sol",          &
        & OuputTypeSolTmp,OuputTypeSolFlg,ierr)

    IF (OuputTypeSolFlg) THEN
        IF (OuputTypeSolTmp.EQ.1) THEN
            OuputType%Sol=OuputTypeSolTmp
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
                & "[ERROR] Ouput_type_sol used is invalid\n",ierr)
            STOP
        END IF
    END IF

END SUBROUTINE GetOuputTypeSol

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

    CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-Run_options_scheme",        &
        RunOptions%Scheme,RunOptionsSchemeFlg,ierr)
    CALL PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-Run_options_time",         &
        RunOptions%Time,RunOptionsTimeFlg,ierr)

    IF (RunOptionsSchemeFlg) THEN
        IF ((RunOptions%Scheme.GE.3).OR.(RunOptions%Scheme.LE.0)) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                     &
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

    CALL PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-verbose",         &
        Verbose,VerboseFlg,ierr)

    IF (.NOT.VerboseFlg) THEN
        Verbose=.TRUE.
    END IF

END SUBROUTINE GetVerbose

END MODULE ANISOFLOW_Interface