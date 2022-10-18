program part1
    use matrixUtilities
    use pseudo
    implicit none

    real*16,dimension(:,:),allocatable::mat,x
    real*16::times
    double precision::error
    integer::iterations
    integer::iterMax

    iterMax=250


    print *,"Se corre la parte 1"
    mat=generatedExampleMatrix()

    x=newton(mat,iterMax,times,iterations,error)
    print *,"_______________________________________________________________________________________"
    !Para mostrar la matriz quitar el comentario en la siguiente linea
    !call printMatrix(x)
    print *, "Metodo de Newton"
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error


    x=chebyshev(mat,iterMax,times,iterations,error)
    print *,"_______________________________________________________________________________________"
    !Para mostrar la matriz quitar el comentario en la siguiente linea
    !call printMatrix(x)
    print *, "Metodo de ChebyShev"
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error

    x=homier(mat,iterMax,times,iterations,error)
    print *,"_______________________________________________________________________________________"
    !Para mostrar la matriz quitar el comentario en la siguiente linea
    !call printMatrix(x)
    print *, "Metodo de Homier"
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error


    x =toutonian(mat,iterMax,times,iterations,error)
    print *,"_______________________________________________________________________________________"
    !Para mostrar la matriz quitar el comentario en la siguiente linea
    !call printMatrix(x)
    print *, "Metodo de Toutonia"
    print *, "Iteraciones: ",iterations
    print *, "Tiempo: ",times
    print *, "Error: ",error
   
    
    
   


    


    
     
    

end program part1