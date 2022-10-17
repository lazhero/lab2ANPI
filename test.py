import numpy as np

def getMatrizExample():
    matrix=[]
    vector=[]

    for i in range(45):
        vector=[]
        for j in range(30):
            vector.append(pow(i+1,2)+pow(j+1,2))
        matrix.append(vector)
    return matrix





def initValue(A):
    At=np.transpose(A)
    M=np.dot(At,A)
    eigenValue=np.linalg.eigvalsh(M)
    myMax=max(eigenValue)
    return np.multiply(At,(1/(pow(myMax,2))))

A=getMatrizExample()
X=initValue(A)




def homeir(A,iterMax):

    X=initValue(A)
    I=np.identity(len(A))
    error=0
    for i in range(iterMax):
        Y=np.dot(A,X)
        s=np.subtract(I,Y)
        t=np.add(I,np.subtract(np.multiply(2,I),Y))
        t=np.dot(t,t)
        X=np.dot(X,np.add(I,np.multiply(0.5,np.multiply(s,t))))
        
        Y=np.dot(A,X)
        error=np.linalg.norm(np.subtract(np.dot(Y,A),A),ord='fro')
        print(error)
        if error<0.00001:
            break
    return X,i,error


def newton(A,iterMax):

    X=initValue(A)
    I=np.identity(len(A))
    error=0
    for i in range(iterMax):
        Y=np.dot(A,X)
       
        X=np.subtract(np.multiply(X,2),np.dot(X,Y))
        
        Y=np.dot(A,X)
        error=np.linalg.norm(np.subtract(np.dot(Y,A),A),ord='fro')
        print(error)
        if error<0.00001:
            break
    return X,i,error




X,iteraciones,error=homeir(getMatrizExample(),250)
print(X)
print(iteraciones)
print(error)