#####################################################   ATENÇÃO   ####################################################
#### Antes de executar, verificar qual é o f(x) usado no código! Ele deve ser compatível com a função em questão. ####
######################################################################################################################

import math
import numpy

## Método da bisseção
def met_bissecao (a, b):
    tol = 10**(-3)
    contador = 1
    
    while abs (a-b) > tol:
        xi = (a+b)/2
        #y = funcao_y_n1 (xi)
        #y = funcao_y_n2 (xi)
        y = funcao (xi)
        
        if y > 0:
            b = xi
            
        else:
            a = xi
            
        contador += 1
            
    print ("Raiz: " + str (round (xi, 5)) + "\n")
    print ("Tolerância atingida: " + str (abs (a-b)) + "\n")
    print ("Número de iterações: " + str (contador))
        
    return ("\n* Fim *")


## Método de Newton
def met_Newton_raiz (x_inicial):
    tol = 10**(-3)
    n_iter = 1000
    
    xi = x_inicial
    
    for k in range (n_iter):
        #razao = funcao_y_n1 (xi)/derivada_y_n1 (xi)
        #razao = funcao_y_n2 (xi)/derivada_y_n2 (xi)
        razao = funcao (xi)/derivada (xi)
        
        xj = xi - razao
        
        tolk = abs (xi - xj)
        
        if tolk < tol:
            print ("Raiz: " + str (round (xj, 5)) + "\n")
            print ("Tolerância atingida: " + str (tolk) + "\n")
            print ("Número de iterações: " + str (k + 1))
            return ("\n* Fim *")
        
        else:
            xi = xj
          
    return ("Convergência não alcançada.")


## Método da secante
def met_secante (x_inicial):
    tol = 10**(-3)
    n_iter = 1000
    delta_x = 0.001
    
    x0 = x_inicial
    x1 = x0 + delta_x
    
    #f0 = funcao_y_n1 (x0)
    #f0 = funcao_y_n2 (x0)
    f0 = funcao (x0)
    
    for k in range (n_iter):
        #f1 = funcao_y_n1 (x1)
        #f1 = funcao_y_n2 (x1)
        f1 = funcao (x1)
        
        razao = (x1 - x0)/(f1 - f0)
        
        x2 = x1 - f1*razao
        
        tolk = abs (x2 - x1)
        
        if tolk < tol:
            print ("Raiz: " + str (round (x1, 5)) + "\n")
            print ("Tolerância atingida: " + str (tolk) + "\n")
            print ("Número de iterações: " + str (k + 1))
            return ("\n* Fim *")
        
        else:
            f0 = f1
            x0 = x1
            x1 = x2
          
    return ("Convergência não alcançada.")


## Método da interpolação inversa
def met_interp_inv (pontos):
    x1, x2, x3 = pontos
    
    x_inicial = math.inf
    tol = 10**(-4)
    n_iter = 1000
    
    for k in range (n_iter):
        #y1 = funcao_y_n1 (x1)
        #y2 = funcao_y_n1 (x2)
        #y3 = funcao_y_n1 (x3)
        y1 = funcao_y_n2 (x1)
        y2 = funcao_y_n2 (x2)
        y3 = funcao_y_n2 (x3)
        
        xk = (y2*y3*x1)/((y1 - y2)*(y1 - y3)) + (y1*y3*x2)/((y2 - y1)*(y2 - y3)) + (y1*y2*x3)/((y3 - y1)*(y3 - y2))
        tolk = abs (xk - x_inicial)
        
        if tolk < tol:
            print ("Raiz: " + str (round (xk, 5)) + "\n")
            print ("Tolerância atingida: " + str (tolk) + "\n")
            print ("Número de iterações: " + str (k + 1))
            return ("\n* Fim *")
        
        else:
            x_inicial = xk
            
            if abs (y3) > abs (y1) and abs (y3) > abs (y2):
                x3 = xk
                
            elif abs (y2) > abs (y1) and abs (y2) > abs (y3):
                x2 = xk
                
            else:
                x1 = xk
                
            l = [x1, x2, x3]
            l.sort ()
            x1, x2, x3 = l
            
    return ("Convergência não alcançada.")


## Método de Newton para sistemas de equações
def met_Newton_sistemas (pontos):
    xi, yi = pontos
    x = [[xi], [yi]]
    
    tol = 10**(-3)
    n_iter = 1000
    
    for k in range (n_iter):
        j = [[0, 0], [0, 0]]
        f = [[0], [0]]
    
        # J
        #derivada_1 = derivada_y1_n3 (x [0][0], x [1][0], x [2][0])
        #derivada_2 = derivada_y2_n3 (x [0][0], x [1][0], x [2][0])
        #derivada_3 = derivada_y3_n3 (x [0][0], x [1][0], x [2][0])
        '''derivada_1 = derivada_y1_n4 (x [0][0], x [1][0], x [2][0])
        derivada_2 = derivada_y2_n4 (x [0][0], x [1][0], x [2][0])
        derivada_3 = derivada_y3_n4 (x [0][0], x [1][0], x [2][0])'''
        
        derivada_1 = derivada_y1 (x [0][0], x [1][0])
        derivada_2 = derivada_y2 (x [0][0], x [1][0])
        
        j [0][0] = derivada_1 [0]
        j [0][1] = derivada_1 [1]
        #j [0][2] = derivada_1 [2]
        j [1][0] = derivada_2 [0]
        j [1][1] = derivada_2 [1]
        '''j [1][2] = derivada_2 [2]
        j [2][0] = derivada_3 [0]
        j [2][1] = derivada_3 [1]
        j [2][2] = derivada_3 [2]'''
        
        j_matriz = numpy.matrix (j)
        j_inversa = numpy.linalg.inv (j_matriz)
        
        # f (Xk)
        #funcao_1 = funcao_y1_n3 (x [0][0], x [1][0], x [2][0])
        #funcao_2 = funcao_y2_n3 (x [0][0], x [1][0], x [2][0])
        #funcao_3 = funcao_y3_n3 (x [0][0], x [1][0], x [2][0])
        
        #th1, th2 = 0, 3
        #th1, th2 = 0.75, 6.5
        '''th1, th2 = 0, 11.667
        funcao_1 = funcao_y1_n4 (x [0][0], x [1][0], x [2][0])
        funcao_2 = funcao_y2_n4 (x [0][0], x [1][0], x [2][0], th1)
        funcao_3 = funcao_y3_n4 (x [0][0], x [1][0], x [2][0], th2)'''
        
        funcao_1 = funcao_y1 (x [0][0], x [1][0])
        funcao_2 = funcao_y2 (x [0][0], x [1][0])
        
        f [0][0] = funcao_1
        f [1][0] = funcao_2
        
        f_matriz = numpy.matrix (f)
        
        # J**(-1) * f(Xk)        
        delta_x = numpy.matmul (j_inversa, f_matriz)
        
        # Xk+1        
        x_novo = numpy.subtract (x, delta_x)
        
        tolk = numpy.linalg.norm (delta_x)/numpy.linalg.norm (x_novo)
        
        if tolk < tol:
            print ("Raízes: " + str (x_novo) + "\n")
            print ("Tolerância atingida: " + str (tolk) + "\n")
            print ("Número de iterações: " + str (k + 1))
            return ("\n* Fim *")
        
        else:
            x_matriz = x_novo
            x = numpy.ndarray.tolist (x_matriz)
            
    return ("Convergência não alcançada.")
    

## Método de Broyden
def met_Broyden (pontos):
    xi, yi = pontos
    x = [[xi], [yi]]
    
    tol = 10**(-3)
    n_iter = 1000
    
    # B0
    j = [[0, 0], [0, 0]]
    
    #derivada_1 = derivada_y1_n3 (x [0][0], x [1][0], x [2][0])
    #derivada_2 = derivada_y2_n3 (x [0][0], x [1][0], x [2][0])
    #derivada_3 = derivada_y3_n3 (x [0][0], x [1][0], x [2][0])
    '''derivada_1 = derivada_y1_n4 (x [0][0], x [1][0], x [2][0])
    derivada_2 = derivada_y2_n4 (x [0][0], x [1][0], x [2][0])
    derivada_3 = derivada_y3_n4 (x [0][0], x [1][0], x [2][0])'''
    
    derivada_1 = derivada_y1 (x [0][0], x [1][0])
    derivada_2 = derivada_y2 (x [0][0], x [1][0])
    
    j [0][0] = derivada_1 [0]
    j [0][1] = derivada_1 [1]
    #j [0][2] = derivada_1 [2]
    j [1][0] = derivada_2 [0]
    j [1][1] = derivada_2 [1]
    #j [1][2] = derivada_2 [2]
    '''j [2][0] = derivada_3 [0]
    j [2][1] = derivada_3 [1]
    j [2][2] = derivada_3 [2]'''
    
    j = [[1, 0], [0, 1]]
    
    j_matriz = numpy.matrix (j)
    b0 = j_matriz
    
    for k in range (n_iter):
        fa = [[0], [0]]

        # fa (Xk)
        #funcao_1a = funcao_y1_n3 (x [0][0], x [1][0], x [2][0])
        #funcao_2a = funcao_y2_n3 (x [0][0], x [1][0], x [2][0])
        #funcao_3a = funcao_y3_n3 (x [0][0], x [1][0], x [2][0])
        
        #th1, th2 = 0, 3
        #th1, th2 = 0.75, 6.5
        '''th1, th2 = 0, 11.667
        funcao_1a = funcao_y1_n4 (x [0][0], x [1][0], x [2][0])
        funcao_2a = funcao_y2_n4 (x [0][0], x [1][0], x [2][0], th1)
        funcao_3a = funcao_y3_n4 (x [0][0], x [1][0], x [2][0], th2)'''
        
        funcao_1a = funcao_y1 (x [0][0], x [1][0])
        funcao_2a = funcao_y2 (x [0][0], x [1][0])
        
        fa [0][0] = funcao_1a
        fa [1][0] = funcao_2a
        #fa [2][0] = funcao_3a
        
        fa_matriz = numpy.matrix (fa)
        
        # B0**(-1) * f(Xk)     
        delta_x = numpy.matmul (numpy.linalg.inv (b0), fa_matriz)
        
        # Xk+1        
        x_novo = numpy.subtract (x, delta_x)
        x_novo_lista = numpy.ndarray.tolist (x_novo)
        
        tolk = numpy.linalg.norm (delta_x)/numpy.linalg.norm (x_novo)
        
        if tolk < tol:
            print ("Raízes: " + str (x_novo) + "\n")
            print ("Tolerância atingida: " + str (tolk) + "\n")
            print ("Número de iterações: " + str (k + 1))
            return ("\n* Fim *")
        
        else:
            x_matriz = x_novo
            x = numpy.ndarray.tolist (x_matriz)
            
            # fb (Xk+1)
            fb = [[0], [0]]
    
            # f (Xk)
            #funcao_1b = funcao_y1_n3 (x_novo_lista [0][0], x_novo_lista [1][0], x_novo_lista [2][0])
            #funcao_2b = funcao_y2_n3 (x_novo_lista [0][0], x_novo_lista [1][0], x_novo_lista [2][0])
            #funcao_3b = funcao_y3_n3 (x_novo_lista [0][0], x_novo_lista [1][0], x_novo_lista [2][0])
            
            '''funcao_1b = funcao_y1_n4 (x_novo_lista [0][0], x_novo_lista [1][0], x_novo_lista [2][0])
            funcao_2b = funcao_y2_n4 (x_novo_lista [0][0], x_novo_lista [1][0], x_novo_lista [2][0], th1)
            funcao_3b = funcao_y3_n4 (x_novo_lista [0][0], x_novo_lista [1][0], x_novo_lista [2][0], th2)'''
            
            funcao_1b = funcao_y1 (x_novo_lista [0][0], x_novo_lista [1][0])
            funcao_2b = funcao_y2 (x_novo_lista [0][0], x_novo_lista [1][0])
            
            fb [0][0] = funcao_1b
            fb [1][0] = funcao_2b
            #fb [2][0] = funcao_3b
            
            fb_matriz = numpy.matrix (fb)
            
            y = numpy.subtract (fb_matriz, fa_matriz)
            
            delta_x_neg = numpy.multiply (delta_x, [[-1]])
            
            pt1 = numpy.subtract (y, numpy.matmul (b0, delta_x_neg))
            pt2 = numpy.transpose (delta_x_neg)
            pt3 = numpy.matmul (pt2, delta_x_neg)
            razao = numpy.divide (numpy.matmul (pt1, pt2), pt3)
            
            # B1
            b1 = numpy.add (b0, razao)
            b0 = b1
            
    return ("Convergência não alcançada.")


## Ajuste de curvas não-lineares
def ajuste (pontos):
    b0, b1, b2 = pontos
    b = [[b0], [b1], [b2]]
    
    tol = 10**(-4)
    n_iter = 1000
    
    for k in range (n_iter):
        j = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        f = [[0], [0], [0]]
    
        # J
        derivada_1 = derivada_y1_n5 (b [0][0], b [1][0], b [2][0])
        derivada_2 = derivada_y2_n5 (b [0][0], b [1][0], b [2][0])
        derivada_3 = derivada_y3_n5 (b [0][0], b [1][0], b [2][0])
        
        j [0][0] = derivada_1 [0]
        j [0][1] = derivada_1 [1]
        j [0][2] = derivada_1 [2]
        j [1][0] = derivada_2 [0]
        j [1][1] = derivada_2 [1]
        j [1][2] = derivada_2 [2]
        j [2][0] = derivada_3 [0]
        j [2][1] = derivada_3 [1]
        j [2][2] = derivada_3 [2]
        
        j_matriz = numpy.matrix (j)
        j_matriz_transposta = numpy.transpose (j_matriz)
        
        # f (Xk)
        funcao_1 = funcao_y1_n5 (b [0][0], b [1][0], b [2][0])
        funcao_2 = funcao_y2_n5 (b [0][0], b [1][0], b [2][0])
        funcao_3 = funcao_y3_n5 (b [0][0], b [1][0], b [2][0])
        
        f [0][0] = funcao_1
        f [1][0] = funcao_2
        f [2][0] = funcao_3
        
        f_matriz = numpy.matrix (f)
        
        # (Jt*J)**(-1) * Jt * f(Xk)
        pt1 = numpy.linalg.inv (numpy.matmul (j_matriz_transposta, j_matriz))
        pt2 = numpy.matmul (j_matriz_transposta, f_matriz)
        delta_b = numpy.matmul (pt1, pt2)
        
        # Xk+1        
        b_novo = numpy.subtract (b, delta_b)
        
        tolk = numpy.linalg.norm (delta_b)/numpy.linalg.norm (b_novo)
        
        if tolk < tol:
            print ("Parâmetros: " + str (b_novo) + "\n")
            print ("Tolerância atingida: " + str (tolk) + "\n")
            print ("Número de iterações: " + str (k + 1))
            return ("\n* Fim *")
        
        else:
            b_matriz = b_novo
            b = numpy.ndarray.tolist (b_matriz)
            
    return ("Convergência não alcançada.")


## Função - Questão 1
pontos_n1 = [-500, -50, -10]

def funcao_y_n1 (xi):
    g = 9.806
    k = 0.00341
    c = math.sqrt(g*k)
    
    return ((math.log (math.cosh (xi*c), 10)) - 50)

def derivada_y_n1 (xi):
    g = 9.806
    k = 0.00341
    c = math.sqrt(g*k)
    
    return c*math.tanh (c*xi)/math.log (10)


## Função - Questão 2
pontos_n2 = [-10.6, -10.5, -10]

def funcao_y_n2 (xi):
    return (4*math.cos (xi) - math.e**(2*xi))

def derivada_y_n2 (xi):
    return (-4*math.sin(xi) - 2*math.e**(2*xi))


## Sistema - Questão 3
pontos_n3 = [1, 2, 3]

def funcao_y1_n3 (xi, yi, zi):
    return (16*xi**4 + 16*yi**4 + zi**4 - 16)

def funcao_y2_n3 (xi, yi, zi):
    return (xi**2 + yi**2 + zi**2 - 3)

def funcao_y3_n3 (xi, yi, zi):
    return (xi**3 - yi + zi - 1)

def derivada_y1_n3 (xi, yi, zi):
    return (64*xi**3, 64*yi**3, 4*zi**3)

def derivada_y2_n3 (xi, yi, zi):
    return (2*xi, 2*yi, 2*zi)

def derivada_y3_n3 (xi, yi, zi):
    return (3*xi**2, -1, 1)


## Sistema - Questão 4
pontos_n4_Newton = [1, 2, 3]
pontos_n4_Broyden = [0.5, 0.5, 0.5]

def funcao_y1_n4 (xi, yi, zi):
    return (2*yi**2 + xi**2 + 6*zi**2 - 1)

def funcao_y2_n4 (xi, yi, zi, th1):
    return (8*yi**3 + 6*yi*xi**2 + 36*yi*xi*zi + 108*yi*zi**2 - th1)

def funcao_y3_n4 (xi, yi, zi, th2):
    return (60*yi**4 + 60*(yi**2)*(xi**2) + 576*(yi**2)*xi*zi + 2232*(yi**2)*(zi**2) + 252*(zi**2)*(xi**2) + 1296*(zi**3)*xi + 3348*(zi)**4 + 24*(xi**3)*zi + 3*xi - th2)

def derivada_y1_n4 (xi, yi, zi):
    return (2*xi, 4*yi, 12*zi)

def derivada_y2_n4 (xi, yi, zi):
    return (12*xi*yi + 36*yi*zi, 24*yi**2 + 6*xi**2 + 36*xi*zi + 108*zi**2, 36*xi*yi + 216*yi*zi)

def derivada_y3_n4 (xi, yi, zi):
    return (120*xi*yi**2 + 576*(yi**2)*zi + 504*xi*zi**2 + 1296*zi**3 + 72*(xi**2)*zi + 3, 240*yi**3 + 120*(xi**2)*yi + 1152*xi*yi*zi + 4464*yi*zi**2, 576*xi*yi**2 + 4464*(yi**2)*zi + 504*(xi**2)*zi + 3888*xi*zi**2 + 13392*zi**3 + 24*xi**3)


## Função - Questão 5
pontos_n5 = [1, 0.25, 3.5]

def funcao_y1_n5 (b0, b1, b2):
    return (b0 + b1*1**b2 - 1)

def funcao_y2_n5 (b0, b1, b2):
    return (b0 + b1*2**b2 - 2)

def funcao_y3_n5 (b0, b1, b2):
    return (b0 + b1*3**b2 - 9)

def derivada_y1_n5 (b0, b1, b2):
    return (1, 1, 0)

def derivada_y2_n5 (b0, b1, b2):
    return (1, 2**b2, b1*(2**b2)*math.log(2))

def derivada_y3_n5 (b0, b1, b2):
    return (1, 3**b2, b1*(3**b2)*math.log(3))


## Sistema - QUESTÃO 7
pontos_Newton = [1, 1]
pontos_Broyden = [1, 1]

def funcao_y1 (xi, yi):
    return (xi - yi + 2)

def funcao_y2 (xi, yi):
    return (math.e**xi + yi - 5)

def derivada_y1 (xi, yi):
    return (1, -1)

def derivada_y2 (xi, yi):
    return (math.e**xi, 1)


## QUESTÃO 4
def funcao (xi):
    return (2**xi - 1/xi - 2*xi)

def derivada (xi):
    return (-2 + 1/xi**2 + 2**xi * math.log(2))


## QUESTÃO 4
#print (met_bissecao (1, 3))
#print (met_Newton_raiz (1))
#print (met_secante (1))
    

## QUESTÃO 7
#print (met_Newton_sistemas (pontos_Newton))
print (met_Broyden (pontos_Broyden))