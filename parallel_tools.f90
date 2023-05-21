!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MOUDLE: parallel_tools
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Author: Tianhong Huang
! Date:   MAY 2023
! Intro:  package contains all subroutines might be used during parallel
!         computing.
! SBRs:   CREATE_SCHEME
!         CREATE_SCHEME_FFT
!         CREATE_SPEC_DRV
!         PARALLEL_ZFFT
!         PARALLEL_IFFT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE PARALLEL_TOOLS
  
    USE MPI
    USE PARALLEL_STRUCTS

    IMPLICIT NONE

    PRIVATE

    INCLUDE 'integer_types.h' !< define standard integer typs.
    INCLUDE 'real_types.h'    !< define standard real-precision types.

    PUBLIC :: CREATE_SCHEME
    PUBLIC :: CREATE_SCHEME_FFT
    PUBLIC :: CREATE_SPEC_DRV
    PUBLIC :: PARALLEL_ZFFT
    PUBLIC :: PARALLEL_IFFT

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: create_scheme
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @param[in] rowCount:    # of total rows.
! @param[in] colSpc:      spacing between columns (dx).
! @param[in] comm:        communicator hosts all processors.
! @param[in] mpiDatatype: datatype of the big grid.
! @param[in] ovlp:        # of extra boundary columns (one-side).
! @param[io] colCount:    # of columns in sub-grid (with ovlps).
! @param[io] colRef:      col-pos of the reference point.
! @param[io] colIdx:      col-idx of the reference point.
! @param[io] sch:         datatype used for grid decmp, parallel fft, etc.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE CREATE_SCHEME(rowCount, colCount, colSpc, colRef, colIdx, &
    comm, mpiDatatype, ovlp, sch)

USE MPI
USE PARALLEL_STRUCTS

IMPLICIT NONE

INTEGER, INTENT(IN) :: rowCount
INTEGER, INTENT(IN) :: comm
INTEGER, INTENT(IN) :: mpiDatatype
INTEGER, INTENT(IN) :: ovlp
INTEGER, INTENT(INOUT) :: colCount
INTEGER, INTENT(INOUT) :: colIdx
TYPE(SCHEME), INTENT(INOUT) :: sch
DOUBLE PRECISION, INTENT(IN) :: colSpc
DOUBLE PRECISION, INTENT(INOUT) :: colRef

INTEGER :: i
INTEGER :: ierror !< argument for MPI subroutine calls.
INTEGER :: stdCount !< standard temp row/col count.

CALL MPI_COMM_RANK(comm, sch%procID, ierror)
CALL MPI_COMM_SIZE(comm, sch%commSize, ierror)

!< set up basic parameters & variables.
sch%ovlp = ovlp
sch%comm = comm
sch%colSpc = colSpc
sch%datatype = mpiDatatype
sch%gridSize(0) = rowCount ! gridSize(0) = xLen of big-grid.
sch%gridSize(1) = colCount ! gridSize(1) = yLen of big-grid.

!< set ovlp = 0 for single processor case.
IF (sch%commSize .EQ. 1) THEN
     sch%ovlp = 0
END IF

!< generate row decomposition array.
!< divide # of rows at most average sense, according to # of threads.
ALLOCATE(sch%rowDcmpSizes(0:sch%commSize-1))
stdCount = rowCount/sch%commSize
DO i = 0, sch%CommSize-1
   IF (i .LT. MOD(rowCount, sch%commSize)) THEN
       sch%rowDcmpSizes(i) = stdCount + 1
   ELSE
       sch%rowDcmpSizes(i) = stdCount
   END IF
END DO

!< generate column decomposition array.
!< divide # of columns at most average sense, according to # of threads.
ALLOCATE(sch%colDcmpSizes(0:sch%commSize-1))
stdCount = colCount/sch%commSize
DO i = 0, sch%CommSize-1
   IF (i .LT. MOD(colCount, sch%commSize)) THEN
       sch%colDcmpSizes(i) = stdCount + 1
   ELSE
       sch%colDcmpSizes(i) = stdCount
   END IF
END DO

!< generate column decomposition array with extra overlap columns.
!< first and last processors' sub-grid only contain one-side bdry condition.
ALLOCATE(sch%colDcmpSizesOvlp(0:sch%commSize-1))
sch%colDcmpSizesOvlp(0) = sch%colDcmpSizes(0) + sch%ovlp
DO i = 1, sch%commSize-2
    sch%colDcmpSizesOvlp(i) = sch%colDcmpSizes(i) + 2*sch%ovlp
END DO
sch%colDcmpSizesOvlp(sch%commSize-1) = sch%colDcmpSizes(sch%commSize-1) &
     + sch%ovlp

!< set up slab sizes, vertically and horizontally.
!< slab size = row # of sub-grid (xLen), col # of sub-grid (yLen).
sch%hSlabSize(0) = sch%rowDcmpSizes(sch%procID)
sch%hSlabSize(1) = sch%gridSize(1)

sch%vSlabSize(0) = sch%gridSize(0)
sch%vSlabSize(1) = sch%colDcmpSizes(sch%procID)

sch%vSlabSizeOvlp(0) = sch%gridSize(0)
sch%vSlabSizeOvlp(1) = sch%colDcmpSizesOvlp(sch%procID)

!< store column indices of interior of vertical slab.
IF (sch%procID .EQ. 0) THEN
    sch%vSlabInt(0) = 0
    sch%vSlabInt(1) = sch%vSlabSizeOvlp(1) - sch%ovlp - 1
ELSE IF (sch%procId .EQ. sch%commSize-1) THEN
    sch%vSlabInt(0) = sch%ovlp
    sch%vSlabInt(1) = sch%vSlabSizeOvlp(1) - 1
ELSE
    sch%vSlabInt(0) = sch%ovlp
    sch%vSlabInt(1) = sch%vSlabSizeOvlp(1) - sch%ovlp - 1
END IF

!< return colCount as # of col in sub-grid (with ovlps).
colCount = sch%vSlabSizeOvlp(1)

!< calculate x-pos of sub-grid reference point.
IF ((sch%commSize .EQ. 1) .OR. (sch%procID .EQ. 0)) THEN
    sch%colRef = colRef
    sch%colIdx = colIdx
ELSE
    sch%colRef = colRef + colSpc * &
        (SUM(sch%colDcmpSizes(0:sch%procID-1)) - sch%ovlp)
    sch%colIdx = colIdx + SUM(sch%colDcmpSizes(0:sch%procID-1)) - sch%ovlp
    colRef = sch%colRef
    colIdx = sch%colIdx
END IF

! create SEND and RECV boundary datatypes.
IF ((sch%commSize .GT. 1) .AND. (sch%ovlp .GT. 0)) THEN
    ! in our model simulations we always have ovlp = 0, just skip this part.
END IF

! set scheme as initialized.
sch%initScheme = .TRUE.

END SUBROUTINE CREATE_SCHEME


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: create_scheme_fft
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @param[io] sch: 1. execute part of CFF2DF subroutine.
!                 2. add datatype subarrays needed for parallel fft.
!                 3. allocate arrays needed for ALLTOALLW comm.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE CREATE_SCHEME_FFT(sch)

USE MPI
USE PARALLEL_STRUCTS

IMPLICIT NONE

TYPE(SCHEME), INTENT(INOUT) :: sch

! create arrays needed for 1D zFFT (part of CFFT2DF).
ALLOCATE(sch%WSAVE1(4*sch%gridSize(0)+15))
ALLOCATE(sch%WSAVE2(4*sch%gridSize(1)+15))
CALL ZFFTI(sch%gridSize(0), sch%WSAVE1)
CALL ZFFTI(sch%gridSize(1), sch%WSAVE2)

! calculate normalization coeffs.
sch%norm_2D = DBLE(1/DBLE(PRODUCT(sch%gridSize)))

! determine size of wave number vectors.
ALLOCATE(sch%wvNum1(0:sch%gridSize(0)-1))
ALLOCATE(sch%wvNum2(0:sch%vSlabSizeOvlp(1)-1))

! create datatype subarrays for ALLTOALLW comm during parallel fft.
ALLOCATE(sch%datatype_array(sch%commSize,0:1))
CALL SUBARRAY_DATATYPE(sch%vSlabSize, 0, sch%commSize, &
    sch%datatype, sch%datatype_array(:,0))
CALL SUBARRAY_DATATYPE(sch%hSlabSize, 1, sch%commSize, &
    sch%datatype, sch%datatype_array(:,1))

! allocate specific arrays needed for ALLTOALLW comm.
ALLOCATE(sch%counts(sch%commSize))
sch%counts = 1
ALLOCATE(sch%displs(sch%commSize))
sch%displs = 0

! set scheme as initialized for zFFT.
sch%initFFT = .TRUE.

CONTAINS
    ! [SBR] AVG_DECMP: generate the most average sense decomposition
    !                  along a single dimension of specific 2D array.
    ! @param[in] nelems: # of elements along target dimension.
    ! @param[in] nparts: # of parts to divide dimension into.
    ! @param[in] pidx:   part idx of output block.
    ! @param[ot] sidx:   start index of output block.
    ! @param[ot] nelemspart: # of elements in output block.
    SUBROUTINE AVG_DECMP(nelems, nparts, pidx, nelemspart, sidx)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: pidx
        INTEGER, INTENT(IN) :: nelems
        INTEGER, INTENT(IN) :: nparts
        INTEGER, INTENT(OUT) :: sidx
        INTEGER, INTENT(OUT) :: nelemspart

        INTEGER :: q, r

        q = nelems / nparts
        r = MOD(nelems, nparts)
        IF (pidx < r) THEN
            nelemspart = q + 1
            sidx = nelemspart * pidx
        ELSE
            nelemspart = q
            sidx = nelemspart * pidx + r
        END IF

    END SUBROUTINE AVG_DECMP

    ! [SBR] SUBARRAY_DATATYPE: create subarray datatypes for local sub-grids.
    ! @param[in] axis:     axis to partition, row = 0, col = 1.
    ! @param[in] sizes:    size of original local sub-grids.
    ! @param[in] nparts:   number of parts divided.
    ! @param[in] datatype: MPI datatype descriptor.
    ! @param[ot] datatype_array: output subarray datatypes.
    SUBROUTINE SUBARRAY_DATATYPE(sizes, axis, nparts, datatype, datatype_array)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: axis
        INTEGER, INTENT(IN) :: nparts
        INTEGER, INTENT(IN) :: datatype
        INTEGER, INTENT(IN) :: sizes(0:1)
        INTEGER, INTENT(OUT) :: datatype_array(0:nparts-1)

        INTEGER :: n, s, p
        INTEGER :: ierror
        INTEGER :: subsizes(0:1)
        INTEGER :: substarts(0:1)

        subsizes = sizes !< set initial sizes of subarrays.
        substarts = 0

        DO p = 0, nparts-1
            CALL AVG_DECMP(sizes(axis), nparts, p, n, s)
            subsizes(axis) = n
            substarts(axis) = s
            CALL MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, substarts, &
                 MPI_ORDER_FORTRAN, datatype, datatype_array(p), ierror)
            CALL MPI_TYPE_COMMIT(datatype_array(p), ierror)
        END DO

    END SUBROUTINE SUBARRAY_DATATYPE
 
END SUBROUTINE CREATE_SCHEME_FFT


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: create_spec_drv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @param[in] order1:   order of derivative in the 1st dimension.
! @param[in] order2:   order of derivative in the 2nd dimension.
! @param[io] spec_drv: spectral derivatives matrix.
! @param[io] sch:      add wavenumber vectors to datatype sch.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE CREATE_SPEC_DRV(order1, order2, spec_drv, sch)

USE MPI
USE PARALLEL_STRUCTS

IMPLICIT NONE

INTEGER, INTENT(IN) :: order1
INTEGER, INTENT(IN) :: order2
DOUBLE COMPLEX, ALLOCATABLE, INTENT(INOUT) :: spec_drv(:,:)
TYPE(SCHEME), INTENT(INOUT) :: sch

INTEGER :: i, j
INTEGER :: stIdx
DOUBLE COMPLEX, ALLOCATABLE :: wvNum1Sub(:)
DOUBLE COMPLEX, ALLOCATABLE :: wvNum2Sub(:)
DOUBLE COMPLEX, ALLOCATABLE :: wvNum2All(:)

! create array for wavenumbers along the 1st dimension.
DO i = 0, sch%gridSize(0)-1
    sch%wvNum1(i) = DCMPLX(0.0, -sch%gridSize(0)/2 + i)
END DO

! create array for wavenumbers along the 2nd dimension.
ALLOCATE(wvNum2All(0:sch%gridSize(1)-1))
DO i = 0, sch%gridSize(1)-1
    wvNum2All(i) = DCMPLX(0.0, -sch%gridSize(1)/2 + i)
END DO

! store the values we need for the right sub-grid.
IF (sch%procID .EQ. 0) THEN
   stIdx = 0
ELSE
   stIdx = sch%colIdx
END IF
DO i = 0, sch%vSlabSizeOvlp(1)-1
   sch%wvNum2(i) = wvNum2All(stIdx+i)
END DO

DEALLOCATE(wvNum2All)

! get the spectral derivatives along the 1st dimension.
ALLOCATE(wvNum1Sub(0:sch%vSlabSizeOvlp(0)-1))
wvNum1Sub = (0.0, 0.0)
IF (order1 .NE. 0) THEN
    wvNum1Sub = sch%wvNum1
    ! if gridsize(0) even and derivative order odd, zero out highest wavenumber
    IF ((MOD(sch%gridSize(0),2) .EQ. 0) .AND. (MOD(order1,2) .EQ. 1)) THEN
        wvNum1Sub(0) = 0.0
    END IF
    ! power up wavenumber vector and adjust floating issue
    wvNum1Sub = wvNum1Sub ** (DBLE(order1))
    IF (MOD(order1, 2) .EQ. 0) THEN
        wvNum1Sub = DBLE(wvNum1Sub)
    ELSE
        wvNum1Sub = DCMPLX(0.0, DBLE((0.0, -1.0)*wvNum1Sub))
    END IF
END IF

! get the spectral derivatives along the 2nd dimension.
ALLOCATE(wvNum2Sub(0:sch%vSlabSizeOvlp(1)-1))
wvNum2Sub = (0.0, 0.0)
IF (order2 .NE. 0) THEN
    wvNum2Sub = sch%wvNum2
    ! if gridsize(1) even and derivative order odd, zero out highest wavenumber.
    IF (sch%procID .EQ. 0) THEN
        IF ((MOD(sch%gridSize(1), 2) .EQ. 0) .AND. (MOD(order2, 2) .EQ. 1)) THEN
            wvNum2Sub(0) = 0.0
        END IF
    END IF
    ! power up wavenumber vector and adjust floating issue
    wvNum2Sub = wvNum2Sub ** (DBLE(order2))
    IF (MOD(order2, 2) .EQ. 0) THEN
        wvNum2Sub = DBLE(wvNum2Sub)
    ELSE
        wvNum2Sub = DCMPLX(0.0, DBLE((0.0, -1.0)*wvNum2Sub))
    END IF
END IF

! get the final spec_drv array.
ALLOCATE(spec_drv(0:sch%vSlabSizeOvlp(0)-1,0:sch%vSlabSizeOvlp(1)-1))
DO i = 0, sch%vSlabSizeOvlp(0)-1
DO j = 0, sch%vSlabSizeOvlp(1)-1
    spec_drv(i,j) = wvNum1Sub(i) + wvNum2Sub(j)
END DO
END DO

! deallocate local temp arrays.
DEALLOCATE(wvNum1Sub)
DEALLOCATE(wvNum2Sub)

END SUBROUTINE CREATE_SPEC_DRV


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: parallel_zfft
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @param[io] subGrid: the local sub-grid data which taking zfft.
! @param[io] sch:     using sub-grid info from sch%.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL_ZFFT(subGrid, sch)

USE MPI
USE PARALLEL_STRUCTS

IMPLICIT NONE

TYPE(SCHEME), INTENT(IN) :: sch
DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(0:,0:)

INTEGER :: i, j
INTEGER :: ierror
DOUBLE COMPLEX, ALLOCATABLE :: arrSign(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: arrTemp(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: factor_x(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: factor_y(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: subGridInt(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: subGridInt_T(:,:)

DOUBLE PRECISION :: pi_dp = 4.0_dp * ATAN(1.0_dp)

! modify input arrays
ALLOCATE(arrTemp(0:sch%vSlabSizeOvlp(0)-1, 0:sch%vSlabSizeOvlp(1)-1))
ALLOCATE(arrSign(0:sch%vSlabSizeOvlp(0)-1, 0:sch%vSlabSizeOvlp(1)-1))

DO i = 0, sch%vSlabSizeOvlp(0)-1
DO j = 0, sch%vSlabSizeOvlp(1)-1
    arrSign(i,j) = DBLE((-1.0,0.0)**(i+j+sch%colIdx))
END DO
END DO
arrTemp = arrSign * subGrid

! zFFT the columns (interior only, part of CFFT2DF).
DO j = sch%vSlabInt(0), sch%vSlabInt(1)
    CALL ZFFTF(sch%vSlabSizeOvlp(0), arrTemp(:,j), sch%WSAVE1)
END DO

! use ALLTOALLW communication to obtain transpose (la Dalcin et. al. 2019).
ALLOCATE(subGridInt(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
ALLOCATE(subGridInt_T(0:sch%hSlabSize(0)-1, 0:sch%hSlabSize(1)-1))
subGridInt = arrTemp(:,sch%vSlabInt(0):sch%vSlabInt(1))

CALL MPI_ALLTOALLW(subGridInt, sch%counts, sch%displs, sch%datatype_array(:,0), &
    subGridInt_T, sch%counts, sch%displs, sch%datatype_array(:,1), sch%comm, ierror)

! zFFT the rows.
DO i = 0, sch%hslabSize(0)-1
   CALL ZFFTF(sch%hSlabSize(1), subGridInt_T(i,:), sch%WSAVE2)
END DO

! transpose and copy back.
CALL MPI_ALLTOALLW(subGridInt_T, sch%counts, sch%displs, sch%datatype_array(:,1), &
    subGridInt, sch%counts, sch%displs, sch%datatype_array(:,0), sch%comm, ierror)

! modify output arrays.
ALLOCATE(factor_x(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
ALLOCATE(factor_y(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
DO i = 0, sch%vSlabSize(0)-1
DO j = 0, sch%vSlabSize(1)-1
    factor_x(i,j) = ZEXP((0.0,-1.0)*pi_dp*(0.5*sch%gridSize(0)-i)) &
                    / sch%gridSize(0)
    factor_y(i,j) = ZEXP((0.0,-1.0)*pi_dp*(0.5*sch%gridSize(1)-j-sch%colIdx)) &
                    / sch%gridSize(1)
    subGridInt(i,j) = factor_x(i,j) * factor_y(i,j) * subGridInt(i,j)
END DO
END DO

subGrid(:,sch%vSlabInt(0): sch%vSlabInt(1)) = subGridInt

! deallocate local temp array.
DEALLOCATE(arrSign)
DEALLOCATE(arrTemp)
DEALLOCATE(factor_x)
DEALLOCATE(factor_y)
DEALLOCATE(subGridInt)
DEALLOCATE(subGridInt_T)

END SUBROUTINE PARALLEL_ZFFT


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: parallel_ifft
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @param[io] subGrid: the local sub-grid data which taking inverse zfft.
! @param[io] sch:     using sub-grid info from sch%.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL_IFFT(subGrid, sch)

USE MPI
USE PARALLEL_STRUCTS

IMPLICIT NONE

TYPE(SCHEME), INTENT(IN) :: sch
DOUBLE COMPLEX, INTENT(INOUT) :: subGrid(0:,0:)

INTEGER :: i, j
INTEGER :: ierror
DOUBLE COMPLEX, ALLOCATABLE :: arrSign(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: arrTemp(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: factor_x(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: factor_y(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: subGridInt(:,:)
DOUBLE COMPLEX, ALLOCATABLE :: subGridInt_T(:,:)

DOUBLE PRECISION :: pi_dp = 4.0_dp * ATAN(1.0_dp)

! modify input arrays
ALLOCATE(arrTemp(0:sch%vSlabSizeOvlp(0)-1, 0:sch%vSlabSizeOvlp(1)-1))
ALLOCATE(arrSign(0:sch%vSlabSizeOvlp(0)-1, 0:sch%vSlabSizeOvlp(1)-1))

DO i = 0, sch%vSlabSizeOvlp(0)-1
DO j = 0, sch%vSlabSizeOvlp(1)-1
    arrSign(i,j) = DBLE((-1.0,0.0)**(i+j+sch%colIdx))
END DO
END DO
arrTemp = arrSign * subGrid

! inverse FFT the columns.
DO j = sch%vSlabInt(0), sch%vSlabInt(1)
    CALL ZFFTB(sch%vSlabSizeOvlp(0), arrTemp(:,j), sch%WSAVE1)
END DO

! use ALLTOALLW communication to obtain transpose (la Dalcin et. al. 2019).
ALLOCATE(subGridInt(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
ALLOCATE(subGridInt_T(0:sch%hSlabSize(0)-1, 0:sch%hSlabSize(1)-1))
subGridInt = arrTemp(:,sch%vSlabInt(0):sch%vSlabInt(1))

CALL MPI_ALLTOALLW(subGridInt, sch%counts, sch%displs, sch%datatype_array(:,0), &
    subGridInt_T, sch%counts, sch%displs, sch%datatype_array(:,1), sch%comm, ierror)

! inverse FFT the rows.
DO i = 0, sch%hSlabSize(0)-1
   CALL ZFFTB(sch%hSlabSize(1), subGridInt_T(i,:), sch%WSAVE2)
END DO

! transpose and copy back
CALL MPI_ALLTOALLW(subGridInt_T, sch%counts, sch%displs, sch%datatype_array(:,1), &
    subGridInt, sch%counts, sch%displs, sch%datatype_array(:,0), sch%comm, ierror)

! multiple norm constant
subGridInt = sch%norm_2D * subGridInt

! modify output arrays.
ALLOCATE(factor_x(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
ALLOCATE(factor_y(0:sch%vSlabSize(0)-1, 0:sch%vSlabSize(1)-1))
DO i = 0, sch%vSlabSize(0)-1
DO j = 0, sch%vSlabSize(1)-1
    factor_x(i,j) = ZEXP((0.0, 1.0)*pi_dp*(0.5*sch%gridSize(0)-i)) &
                    * sch%gridSize(0)
    factor_y(i,j) = ZEXP((0.0, 1.0)*pi_dp*(0.5*sch%gridSize(1)-j-sch%colIdx)) &
                    * sch%gridSize(1)
    subGridInt(i,j) = factor_x(i,j) * factor_y(i,j) * subGridInt(i,j)
END DO
END DO

subGrid(:,sch%vSlabInt(0): sch%vSlabInt(1)) = subGridInt

! deallocate local temp array.
DEALLOCATE(arrSign)
DEALLOCATE(arrTemp)
DEALLOCATE(factor_x)
DEALLOCATE(factor_y)
DEALLOCATE(subGridInt)
DEALLOCATE(subGridInt_T)

END SUBROUTINE PARALLEL_IFFT


END MODULE PARALLEL_TOOLS
