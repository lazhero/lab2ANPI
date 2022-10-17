include 'arpack_wrapper.f90'
module matrixUtilities
    use arpack_eig
    implicit none
    contains 
        function  getMatrixDimensions(matrix)
            implicit none
            real,dimension(:,:)::matrix
            integer, dimension(2) :: getMatrixDimensions
        
        
            getMatrixDimensions=shape(matrix)
        
        end function
        
        
        function getIdentity(n) result(matrix)
        
            implicit none
        
            !Declaraciones
        
            integer::n
            real,dimension(n,n)::matrix
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

        function mulByNumber(matrix,number) result(resultMatrix)

            
            real,dimension(:,:)::matrix
            real::number
            real,dimension(:,:),allocatable::resultMatrix
            integer,dimension(2)::matrixSize
            integer::n,m,i,k


            matrixSize=getMatrixDimensions(matrix)
            n=matrixSize(1)
            m=matrixSize(2)
            allocate(resultMatrix(n,m))

            

           

            do i=1,n
                do k=1,m
                    resultMatrix(i,k)=number* matrix(i,k)
                end do
            end do

          


        end function

        function addMatrix(matrix1,matrix2) result(resultMatrix)

            real,dimension(:,:)::matrix1,matrix2
            integer,dimension(2)::dims
            real,dimension(:,:),allocatable::resultMatrix
            integer::i,k

            dims=getMatrixDimensions(matrix1)

            allocate(resultMAtrix(dims(1),dims(2)))
            do i=1,dims(1)
                do k=1,dims(2)
                    resultMatrix(i,k)=matrix1(i,k)+matrix2(i,k)
                end do
            end do
            
        end function

        function substractMatrix(matrix1,matrix2) result(resultMatrix)
            real,dimension(:,:)::matrix1,matrix2
            real,dimension(:,:),allocatable::resultMatrix,invMatrix2

            invMAtrix2=mulByNumber(matrix2,-1.0000000000)
            resultMatrix=addMAtrix(matrix1,invMatrix2)
            
        end function

        function mulMatrix(matrix1,matrix2) result(resultMatrix)

            real,dimension(:,:)::matrix1,matrix2
            real,dimension(:,:),allocatable::resultMatrix

            resultMatrix=matmul(matrix1,matrix2)



        end function


        function transpuesta(matrix) result(resultMatrix)

            real,dimension(:,:)::matrix
            real,dimension(:,:),allocatable::resultMatrix
            integer,dimension(2)::dims
            integer::i,k

            dims=getMatrixDimensions(matrix)
            allocate(resultMatrix(dims(1),dims(2)))
            
            do i=1,dims(1)

                do k=1,dims(2)
                    resultMatrix(k,i)=matrix(i,k)

                end do
            end do
        end function

        function norm(matrix) 

            real,dimension(:,:)::matrix
            real::norm
            integer,dimension(2)::dims
            integer::i,k
            norm=0.0


            do i=1,dims(1)
                do k=1,dims(1)
                    norm =norm + abs(matrix(i,k))**2

                end do

            end do
            norm=sqrt(norm)


        end function

        function initialValue(matrix) result(X0)

            real,dimension(:,:)::matrix
            real,dimension(:,:),allocatable::matrixT
            real,dimension(:,:),allocatable::X0
            integer,dimension(2)::dims

            dims=getMatrixDimensions(matrix)


            matrixT=transpuesta(matrix)
            X0=matrix


        end function

        subroutine printMatrix(matrix)
            real, dimension (:,:) :: matrix  
            integer, dimension(2) :: matrixSize
            integer :: i
    
            matrixSize =  getMatrixDimensions(matrix)
    
            do i = 1, matrixSize(1)
    
                write(*,*)  (matrix(i,:))
            end do
    
    
        end subroutine

       

end module matrixUtilities





