program part1
    use arpack_eig

    implicit none

    integer, dimension(45,30) :: matrix
    integer, dimension(30,45) :: T_matrix
    real :: eigval(45, 1), eigvec(1, 45)

    matrix = generatedExampleMatrix()
    T_matrix = TRANSPOSE(matrix)

    A = matmul(matrix, T_matrix)

    call find_eigens(eigval, eigvec, 45, 1, 'LM')

    print *, eigval

    print *, "---------------------------------------------------"
    print *,size(multiplyMatrixByScalar((1/(eigval(1,1))), T_matrix))

    contains
        function generatedExampleMatrix() result (generatedMatrix)
            implicit none
            integer :: i,j
            integer, dimension(45,30) :: generatedMatrix

            do i = 1,45
                do j = 1,30
                    generatedMatrix(i,j) = i**2 + j**2
                end do
            end do
        
            
        end function generatedExampleMatrix

        function multiplyMatrixByScalar(scalar, Tmatrix) result (scaledMatrix)
            real, intent(in) :: scalar
            real, dimension(30,45) :: scaledMatrix
            integer, dimension(30,45), intent(in) :: Tmatrix
            integer :: i, j

            do i = 1,30
                do j = 1,45
                    scaledMatrix(i, j) = 1*Tmatrix(i,j)
                end do
            end do

        end function multiplyMatrixByScalar
    

end program part1