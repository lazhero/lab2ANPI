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
        matrix(1,1) = 1

        matrix(2,2) = 1

        matrix(3,1) = -1
        matrix(3,2) = -1
        matrix(3,3) = 2

        matrix(4,2) = -2
        matrix(4,4) = 3
        matrix(4,5) = -1

        matrix(5,4) = -2
        matrix(5,5) = 7
        matrix(5,6) = -1
        matrix(5,10) = -4

        matrix(6,5) = -1
        matrix(6,6) = 2
        matrix(6,7) = -1

        matrix(7,1) = -12
        matrix(7,6) = -3
        matrix(7,7) = 19
        matrix(7,8) = -4

        matrix(8,7) = -1
        matrix(8,8) = 2
        matrix(8,9) = -1

        matrix(9,8) = -1
        matrix(9,9) = 4
        matrix(9,10) = -3

        matrix(10,5) = -1
        matrix(10,9) = -1
        matrix(10,10) = 2

    end function


    function generateMatrixSolution() result(matrix)
        implicit none
        real*16,dimension(10,1)::matrix  !matriz de =mx1
        integer::m
        integer::i

        m=10

     
        matrix(2,1)=12

       

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