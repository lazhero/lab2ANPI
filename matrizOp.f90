

!
!
!Laboratorio 2 ANPI
!
!
! Luis Andrey Zúñiga
! Brian Wagemans Alvarado
! Adrian Gonzalez  Jimenez
!
module matrixUtilities
    use arpack_eig
    implicit none
    contains 
        
        !Descripcion: Obtiene una matriz cuadrada con n filas y n columnas
        !Entradas: 
        !             *n un número entero positivo indicando el numero de files y columnas
        ! Salida: 
        !             *Matriz es la matriz identidad solicitada
        !
        
        function getIdentity(n) result(matrix)
        
            implicit none
        
            !Declaraciones
        
            integer::n
            real*16,dimension(n,n)::matrix
            integer :: i,k
        
            do i=1,n
                do k=1,n
        
                    if(i==k) then
                        matrix(i,k)=1 ! si se encuentra en la diagona
                    else  
                        matrix(i,k)=0 ! si no esta en la diagonal
                    end if
                    
                end do
        
            end do
        
        
        end function

  

        
       
        ! Descripcion: Obtiene la norma de frobenius para una matriz indicada
        ! Entradas: 
        !           *matriz: la matriz a la cual calcular la norma
        ! Salidas:
        !           *norm: la norma de frobenius
        !
        function norm(matrix) 

            real*16,dimension(:,:)::matrix
            double precision::norm
            integer,dimension(2)::dims
            integer::i,k
            norm=0.0d0

            dims=shape(matrix)
            do i=1,dims(1)
                do k=1,dims(2)
                    norm =norm + abs(matrix(i,k))**2 !calculo de la suma
                    
                end do

            end do
            norm=sqrt(norm)
        


        end function


        ! Descripcion: Calcula el valor X0 para los métodos iterativos para calcular la pseudo inversa
        ! Entradas:
        !           * matriz: La matriz en base a la cual calcular el valor inicial de la aproximacion
        ! Salidas:
        !           * X0: La aproximacion
        !
        function initialValue(matrix) result(X0)
            
            real*16,dimension(:,:)::matrix
            real*16,dimension(:,:),allocatable::X0
            integer,dimension(2)::dims
            real*16::maxProper

            dims=shape(matrix) !obtiene las dimensiones de la matriz

            maxProper=getMaxProperValue(matrix,dims(1),dims(2))
            
            X0=((1/(maxProper**2))*transpose(matrix))
    
           

        end function

        ! Descripcion: Obtiene el maximo valor propio de una matriz indicada de dimensiones mxn
        ! Entradas:
        !           * matrix: La matrix en la cual se debe determinar el maximo valor propio
        !           * m: el numero de filas de la matriz
        !           * n: el numero de columnas de la matriz
        !
        !
        !

        function getMaxProperValue(matrix,m,n) result(value)

            use arpack_eig
            implicit none
            integer::m,n
            real*16::value
            real*16, dimension(m,n) :: matrix
            real*16, dimension(n,m) :: T_matrix
            real :: eigval(m, 1), eigvec(1, m)

            T_matrix=transpose(matrix)

            

            A = matmul(T_matrix, matrix) !matriz A para el calculo de los Valores propios de AT*A
            
            !Este es el metodo creado por Melkiades en  https://gist.github.com/Melkiades/485d9dc3d0fc4c6a9d8ad47466f70913
            call find_eigens(eigval, eigvec, n, 1, 'LM')

            value=maxval(eigval) !obtiene el valor maximo de los posibles

        end function
        !
        ! Descripcion: Es utilizado para imprimir una matriz dada
        ! Entradas: 
        !           *matriz: la matriz a imprimir
        !
        !
        !
        !
        subroutine printMatrix(matrix)
            real*16 ,dimension (:,:),intent(in) :: matrix  
            integer, dimension(2) :: matrixSize
            integer :: i
    
            matrixSize =  shape(matrix)
    
            do i = 1, matrixSize(1)
    
                write(*,*)  (matrix(i,:))
            end do
    
    
        end subroutine

        ! Descripcion: Funcion utilzado para generar la matriz utilizada para el lab
        ! Entradas: Ninguna
        ! Salidas: 
        !          *generatedMatrix: La matriz para el ejemplo
        !
        !
        !
        !
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



