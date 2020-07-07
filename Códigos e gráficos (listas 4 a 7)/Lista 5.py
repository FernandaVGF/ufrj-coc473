#####################################################   ATENÇÃO   ####################################################
#### Antes de executar, verificar qual é o f(x) usado no código! Ele deve ser compatível com a função em questão. ####
######################################################################################################################

import math
import numpy

## Integração polinomial
def int_poli (n, ponto_a, ponto_b):
    # X (nós de integração)
    x = [[0] for r in range (n)]
    
    if n == 1:
        x [0][0] = (ponto_a + ponto_b)/2
        
    else:
        delta = (ponto_b - ponto_a)/(n-1)
        
        x [0][0] = ponto_a
        
        for i in range (1, n):
            x [i][0] = ponto_a + i*delta
   
    # V (matriz de Vandermonde)
    v = [[0 for r in range (n)] for t in range (n)]
    
    for i in range (n):
        for j in range (n):
            v [i][j] = (x [j][0])**i
            
    # B (lado direito)
    b = [[0] for r in range (n)]
    
    for i in range (n):
        b [i][0] = (ponto_b**(i+1) - ponto_a**(i+1))/(i+1)
        
    # W (vetor de pesos)
    v_matriz = numpy.matrix (v)
    b_matriz = numpy.matrix (b)
    
    w = numpy.matmul (numpy.linalg.inv (v_matriz), b_matriz)
    
    w_lista = numpy.ndarray.tolist (w)
    
    # Integral
    soma = 0
    
    for i in range (n):
        #funcao = funcao_y_n2 (x [i][0])
        #funcao = funcao_y1_n3 (x [i][0])
        #funcao = funcao_y2_n3 (x [i][0])
        #funcao = funcao_y1_n4 (x [i][0])
        #funcao = funcao_y2_n4 (x [i][0])
        #funcao = funcao_y_n5 (x [i][0])
        funcao = funcao_y_n6 (x [i][0])
        
        fi = funcao
        soma += w_lista [i][0]*fi
       
    return soma


## Quadratura de Gauss  
# Pesos
pesos = [
        [[2], [0]],
        [[1, 1], [-0.5773502691896257, 0.5773502691896257]],
        [[0.8888888888888888, 0.5555555555555556, 0.5555555555555556], [0, -0.7745966692414834, 0.7745966692414834]],
        [[0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538], [0.3399810435848563, -0.3399810435848563, 0.8611363115940526, -0.8611363115940526]],
        [[0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891], [0, -0.5384693101056831, 0.5384693101056831, 0.9061798459386640, -0.9061798459386640]],
        [[0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704], [0.6612093864662645, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 0.9324695142031521]],
        [[0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766, 0.1294849661688697, 0.1294849661688697], [0, 0.4058451513773972, -0.4058451513773972, -0.7415311855993945, 0.7415311855993945, -0.9491079123427585, 0.9491079123427585]],
        [[0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.3137066458778873, 0.2223810344533745, 0.2223810344533745, 0.1012285362903763, 0.1012285362903763], [-0.1834346424956498, 0.1834346424956498, -0.5255324099163290, 0.5255324099163290, 0.7966664774136267, -0.7966664774136267, -0.9602898564975363, 0.9602898564975363]],
        [[0.3302393550012598, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744, 0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354], [0, -0.8360311073266358,0.8360311073266358, -0.9681602395076261, 0.9681602395076261, -0.3242534234038089, 0.3242534234038089, -0.6133714327005904, 0.6133714327005904]],
        [[0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2692667193099963, 0.2190863625159820, 0.2190863625159820, 0.1494513491505806, 0.1494513491505806, 0.0666713443086881, 0.0666713443086881], [-0.1488743389816312, 0.1488743389816312, -0.4333953941292472, 0.4333953941292472, -0.6794095682990244, 0.6794095682990244, -0.8650633666889845, 0.8650633666889845, -0.9739065285171717, 0.9739065285171717]]
        ]

# Função
def quad_Gauss (n, a, b):
    soma = 0
    l = b - a

    for i in range (n):
        zi = pesos [n-1][1][i]
        wi = pesos [n-1][0][i]
        xi = 0.5*(a + b + zi*l)
        
        #funcao = funcao_y_n2 (xi)
        #funcao = funcao_y1_n3 (xi)
        #funcao = funcao_y2_n3 (xi)
        #funcao = funcao_y1_n4 (xi)
        #funcao = funcao_y2_n4 (xi)
        #funcao = funcao_y_n5 (xi)
        funcao = funcao_y_n6 (xi)
        
        soma += funcao*wi

    soma = 0.5*l*soma
    
    return soma


## Função (Polinomial + Gauss)
def escolha_usuario (opcao, n, a, b):
    if opcao == 0:
        return int_poli (n, a, b)
    
    elif opcao == 1:
        return quad_Gauss (n, a, b)


## Função - Questão 2
def funcao_y_n2 (x):
    return (1/math.sqrt (2*math.pi)) * math.e**(-(x**2)/2)


## Funções - Questão 3
def funcao_y1_n3 (w):
    sn = 2
    c1 = 1
    c2 = 0.05
    return (sn/(math.sqrt((1 - (w/c1)**2)**2 + (2*c2*(w/c1))**2))**2)

def funcao_y2_n3 (w):
    sn = 2
    c1 = 1
    c2 = 0.05
    return (w**2)*(sn/(math.sqrt((1 - (w/c1)**2)**2 + (2*c2*(w/c1))**2))**2)


## Função - Questão 4
def funcao_y1_n4 (w):
    hs = 3
    tz = 5
    sn = 4 * math.pi**3 * hs**2 * math.e**((-16 * math.pi**3)/(w**4 * tz**4)) / (w**5 * tz **4)
    c1 = 1
    c2 = 0.05
    return (sn/(math.sqrt((1 - (w/c1)**2)**2 + (2*c2*(w/c1))**2))**2)

def funcao_y2_n4 (w):
    hs = 3
    tz = 5
    sn = 4 * math.pi**3 * hs**2 * math.e**((-16 * math.pi**3)/(w**4 * tz**4)) / (w**5 * tz **4)
    c1 = 1
    c2 = 0.05
    return (w**2)*(sn/(math.sqrt((1 - (w/c1)**2)**2 + (2*c2*(w/c1))**2))**2)
    

## Função - Questão 5
def funcao_y_n5 (x):
    return (2 + 2*x - x**2 + 3*x**3)


## Função - Questão 6
def funcao_y_n6 (x):
    return (1/(1 + x**2))


## Função - QUESTÃO 5
def funcao_y (x):
    return (math.e**((-x**2)/2))


n = int (input ("Número de pontos de integração (1 a 10): "))
opcao = 0

## Questão 2
#a, b = 0, 1
#a, b = 0, 5
#print (escolha_usuario (opcao, n, a, b))


## Questão 3
#a, b = 0, 10
#print (escolha_usuario (opcao, n, a, b))


## Questão 4
#a, b = 0, 10
#print (escolha_usuario (opcao, n, a, b))


## QUESTÃO 5
a, b = 0, 3
print (escolha_usuario (opcao, n, a, b))


## Questão 5
referencia = 194.67
'''for n in range (1, 100):
    a, b = 0, 4
    
    if round (escolha_usuario (opcao, n, a, b), 2) == referencia:
        print ("Valor da integral: " + str (escolha_usuario (opcao, n, a, b)))
        print ("\nMenor número de pontos de integração: " + str (n))
        break'''
        

## Questão 6
#a, b = 0, 3
#print ("Polinomial/Gauss: " + str (escolha_usuario (opcao, n, a, b)))

'''m = (b + a)/2

# Regra do Ponto Médio
mf = funcao_y_n6 (m)*(b - a)
        
# Regra do Trapézio
tf = (funcao_y_n6 (a) + funcao_y_n6 (b))*(b - a)/2
    
# Regra de Simpson
sf = (funcao_y_n6 (a) + 4*funcao_y_n6 (m) + funcao_y_n6 (b))*(b - a)/6'''

#print ("Regra do Ponto Médio: " + str (mf))
#print ("Regra do Trapézio: " + str (tf))
#print ("Regra de Simpson: " + str (sf))