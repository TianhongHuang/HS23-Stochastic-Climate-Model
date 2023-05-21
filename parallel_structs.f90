!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MOUDLE: parallel_tools
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Author: Tianhong Huang
! Date:   MAY 2023
! Intro:  define and store all parallel compute related paraneters and datatype
!         arrays.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE PARALLEL_STRUCTS
  
    IMPLICIT NONE

    PRIVATE

    TYPE, PUBLIC :: SCHEME

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! large grid (global) variables:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! MPI communicator hosts all processors.
        INTEGER :: comm
        ! # of extra boundary columns (one-side).
        INTEGER :: ovlp
        ! # of processors in the MPI communicator.
        INTEGER :: commSize
        ! MPI datatype of the array.
        INTEGER :: datatype
        ! size of large grid = (xlen, ylen).
        INTEGER :: gridSize(0:1)
        ! space between columns = dx.
        DOUBLE PRECISION :: colSpc
        ! # of rows in each blocks of horizontal slab decomposition.
        INTEGER, ALLOCATABLE :: rowDcmpSizes(:)
        ! # of cols in each blocks of vertical slab decomposition.
        INTEGER, ALLOCATABLE :: colDcmpSizes(:)
        ! # of cols in each blocks of vertical slab decomposition,
        ! contained all the boundary overlaps.
        INTEGER, ALLOCATABLE :: colDcmpSizesOvlp(:)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! sub grids (local) variables:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! processor ID #.
        INTEGER :: procID
        ! column index of reference point.
        INTEGER :: colIdx
        ! column indices of the interior of the vertical slab.
        INTEGER :: vSlabInt(0:1)
        ! size of the sub-grid in the horizontal slab decomposition.
        INTEGER :: hSlabSize(0:1)
        ! size of the sub-grid in the vertical slab decomposition.
        INTEGER :: vSlabSize(0:1)
        ! size of the sub-grid in the vertical slab decomposition,
        ! contained all the boundary overlaps.
        INTEGER :: vSlabSizeOvlp(0:1)
        ! MPI derived datatype for sending boundaries to neightbors,
        ! 0 = left, 1 = right.
        INTEGER :: SEND_BOUNDARIES(0:1)
        ! MPI derived datatype for recieving boundaries from neightbors.
        INTEGER :: RECV_BOUNDARIES(0:1)
        ! physical position of the reference point along the rows.
        DOUBLE PRECISION :: colRef

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! additional fft and spectral derivative variables:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! normalization coefficient for 2D zFFT.
        DOUBLE PRECISION :: norm_2D
        ! holds initialization info for DFFTPACK 1D row FFTs.
        DOUBLE PRECISION, ALLOCATABLE :: WSAVE1(:)
        ! holds initialization info for DFFTPACK 1D col FFTs.
        DOUBLE PRECISION, ALLOCATABLE :: WSAVE2(:)
        ! holds the datatypes needed for ALLTOALLW communication.
        INTEGER, ALLOCATABLE :: datatype_array(:,:)
        ! arrays for using in global re-distribution.
        INTEGER, ALLOCATABLE :: counts(:), displs(:)
        ! holds wavenumber vector in 1st dim, i*k.
        DOUBLE COMPLEX, ALLOCATABLE :: wvNum1(:)
        ! holds wavenumber vector in 2nd dim, i*l.
        DOUBLE COMPLEX, ALLOCATABLE :: wvNum2(:)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! initial flags:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         LOGICAL :: initScheme = .FALSE.

         LOGICAL :: initFFT = .FALSE.

    END TYPE SCHEME
  
END MODULE PARALLEL_STRUCTS
