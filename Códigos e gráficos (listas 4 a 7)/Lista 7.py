#####################################################   ATENÇÃO   ####################################################
#### Antes de executar, verificar qual é o f(x) usado no código! Ele deve ser compatível com a função em questão. ####
######################################################################################################################

import math

## Diferença central
def dif_central (x, v_analitico):
    delta_x = 0.01
    
    numerador = funcao_y1 (x + delta_x) - funcao_y1 (x - delta_x)
    #numerador = funcao_y2 (x + delta_x) - funcao_y2 (x - delta_x)
    #numerador = funcao_y3 (x + delta_x) - funcao_y3 (x - delta_x)
    
    derivada = numerador/(2*delta_x)
            
    print ("Valor encontrado: " + str (derivada) + "\n")
    print ("Valor analítico: " + str (v_analitico) + "\n")
        
    return ("* Fim *")


## Passo para frente
def passo_frente (x, v_analitico):
    delta_x = 0.01
    
    numerador = funcao_y1 (x + delta_x) - funcao_y1 (x)
    #numerador = funcao_y2 (x + delta_x) - funcao_y2 (x)
    #numerador = funcao_y3 (x + delta_x) - funcao_y3 (x)
    
    derivada = numerador/delta_x
            
    print ("Valor encontrado: " + str (derivada) + "\n")
    print ("Valor analítico: " + str (v_analitico) + "\n")
        
    return ("* Fim *")


## Passo para frente
def passo_atras (x, v_analitico):
    delta_x = 0.01
    
    numerador = funcao_y1 (x) - funcao_y1 (x - delta_x)
    #numerador = funcao_y2 (x) - funcao_y2 (x - delta_x)
    #numerador = funcao_y3 (x) - funcao_y3 (x - delta_x)
    
    derivada = numerador/delta_x
            
    print ("Valor encontrado: " + str (derivada) + "\n")
    print ("Valor analítico: " + str (v_analitico) + "\n")
        
    return ("* Fim *")


## Interpolação de Richard
def interp_Richard (x, p, v_analitico):
    delta_x1 = 0.5
    delta_x2 = 0.25
    
    #numerador_d1 = funcao_y1 (x + delta_x1) - funcao_y1 (x)
    #numerador_d1 = funcao_y2 (x + delta_x1) - funcao_y2 (x)
    numerador_d1 = funcao_y3 (x + delta_x1) - funcao_y3 (x)
    
    d1 = numerador_d1/delta_x1
    
    #numerador_d2 = funcao_y1 (x + delta_x2) - funcao_y1 (x)
    #numerador_d2 = funcao_y2 (x + delta_x2) - funcao_y2 (x)
    numerador_d2 = funcao_y3 (x + delta_x2) - funcao_y3 (x)
    
    d2 = numerador_d2/delta_x2
    
    q = delta_x1/delta_x2
    
    rch = d1 + (d1 - d2)/(1/q**p - 1)
    
    print ("Valor encontrado: " + str (rch) + "\n")
    print ("Valor analítico: " + str (v_analitico) + "\n")
        
    return ("* Fim *")


## Funções
#v_analitico = 26.95021
def funcao_y1 (x):
    return (x**3 + 1/(math.e**x))

#v_analitico = 0.70999 
def funcao_y2 (x):
    return (x**(1/3) + math.log (x))

v_analitico = 0.11373
def funcao_y3 (x):
    return (1 - math.e**(-(x/5)**2))


## Questão 1
#x = 3   
#x = 2   
x = 6

p = 1
#p = 2

#print (dif_central (x, v_analitico))
#print (passo_frente (x, v_analitico))
#print (passo_atras (x, v_analitico))
print (interp_Richard (x, p, v_analitico))