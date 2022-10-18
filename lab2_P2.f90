module part2_example
    use matrixUtilities
    implicit none
    

    contains 

    function generateMatrixAplication() result(matrix)
        implicit none
        real*16,dimension(10,10)::matrix  !matriz de mxn
        integer::m,n
        integer::i,k

        m=10
        n=10

        do i=1,m
            do k=1,n

                matrix(i,k)=5  !algun valor que quiera poner
            end do

        end do
        

    end function


    function generateMatrixSolution() result(matrix)
        implicit none
        real*16,dimension(10,1)::matrix  !matriz de =mx1
        integer::m
        integer::i

        m=10

     
        do i=1,m

            matrix(i,1)=5  !algun valor que quiera poner
        end do

       

    end function


end module part2_example


program part2
    use matrixUtilities
    use part2_example
    use pseudo
    implicit none

    real*16,dimension(:,:),allocatable::mat,pse,x,b ! mat *x=b
    real*16::times
    double precision::error
    integer::iterations
    integer::iterMax

    iterMax=250


    print *,"Se corre la parte 2"
    mat=generateMatrixAplication()
    b=generateMatrixSolution()

    pse=newton(mat,iterMax,times,iterations,error)
    print *,"_______________________________________________________________________________________"
    !Para mostrar la matriz quitar el comentario en la siguiente linea
    !call printMatrix(x)
    print *, "Metodo de Newton"
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error
    print *,"La solucion corresponde a:"
    x=matmul(pse,b)
    call printMatrix(x)


    pse=chebyshev(mat,iterMax,times,iterations,error)
    print *,"_______________________________________________________________________________________"
    !Para mostrar la matriz quitar el comentario en la siguiente linea
    !call printMatrix(x)
    print *, "Metodo de ChebyShev"
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error
    print *,"La solucion corresponde a:"
    x=matmul(pse,b)
    call printMatrix(x)


    pse=homier(mat,iterMax,times,iterations,error)
    print *,"_______________________________________________________________________________________"
    !Para mostrar la matriz quitar el comentario en la siguiente linea
    !call printMatrix(x)
    print *, "Metodo de Homier"
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error
    print *,"La solucion corresponde a:"
    x=matmul(pse,b)
    call printMatrix(x)


    pse =toutonian(mat,iterMax,times,iterations,error)
    print *,"_______________________________________________________________________________________"
    !Para mostrar la matriz quitar el comentario en la siguiente linea
    !call printMatrix(x)
    print *, "Metodo de Toutonia"
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error
    print *,"La solucion corresponde a:"
    x=matmul(pse,b)
    call printMatrix(x)
   
    
    
   


    


    
     
    

end program part2