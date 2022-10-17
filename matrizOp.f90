
module matrixUtilities
    use arpack_eig
    implicit none
    contains 
        
        
        
        function getIdentity(n) result(matrix)
        
            implicit none
        
            !Declaraciones
        
            integer::n
            real*16,dimension(n,n)::matrix
            integer :: i,k
        
            do i=1,n
                do k=1,n
        
                    if(i==k) then
                        matrix(i,k)=1
                    else  
                        matrix(i,k)=0
                    end if
                    
                end do
        
            end do
        
        
        end function

  

        
       

        function norm(matrix) 

            real*16,dimension(:,:)::matrix
            double complex::norm
            integer,dimension(2)::dims
            integer::i,k
            norm=0.0


            do i=1,dims(1)
                do k=1,dims(2)
                    norm =norm + abs(matrix(i,k))**2

                end do

            end do
            norm=sqrt(norm)


        end function

        function initialValue(matrix) result(X0)
            
            real*16,dimension(:,:)::matrix
            real*16,dimension(:,:),allocatable::X0
            integer,dimension(2)::dims
            real*16::maxProper

            dims=shape(matrix)

            maxProper=getMaxProperValue(matrix,dims(1),dims(2))
            X0=((1/(maxProper**2))*transpose(matrix))

        end function

        function getMaxProperValue(matrix,m,n) result(value)

            use arpack_eig
            implicit none
            integer::m,n
            real*16::value
            real*16, dimension(m,n) :: matrix
            real*16, dimension(n,m) :: T_matrix
            real :: eigval(m, 1), eigvec(1, m)

            T_matrix=transpose(matrix)

            

            A = matmul(T_matrix, matrix)

            call find_eigens(eigval, eigvec, n, 1, 'LM')

            value=maxval(eigval)

        end function

        subroutine printMatrix(matrix)
            real*16 ,dimension (:,:),intent(in) :: matrix  
            integer, dimension(2) :: matrixSize
            integer :: i
    
            matrixSize =  shape(matrix)
    
            do i = 1, matrixSize(1)
    
                write(*,*)  (matrix(i,:))
            end do
    
    
        end subroutine
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

       

end module matrixUtilities



