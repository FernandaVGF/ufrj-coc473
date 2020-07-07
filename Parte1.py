from math import floor, sin, cos, pi, sqrt, atan

class Matriz:

    def __init__(self, m=0, n=0, l=[], todos_objetos_float=False):
        if (m!=0 and n!=0 and l!=[]):
            self.definir(m,n,l, todos_objetos_float)

    def definir(self, m, n, l, todos_objetos_float=False):
        self.n = n
        self.m = m
        if (type(l) == type([])):
            if (self.n * self.m == len(l)):
                self.l = l
                if todos_objetos_float == False:
                    for contador, valor in enumerate(self.l):
                        self.l[contador] = float(valor)
            else:
                self.n = 0
                self.m = 0
                return -1
        elif (l=="I"):
            self.l = []
            for i in range(self.numeroElementos()):
                if (self.ij_posicao(i)[0]==self.ij_posicao(i)[1]):
                    self.l.append(float(1))
                else:
                    self.l.append(float(0))
        else:
            self.l = [float(l) for i in range(self.numeroElementos())]

    def ij_posicao(self, pos):
        '''
        Retorna a linha e coluna do elemento na posição pos da lista
        :param pos: posição
        :return: [linha, coluna]
        '''
        div = floor((pos*1.0)/self.m)

        return [div+1,-1*div*self.m + pos + 1]

    def mult(self, m1, m2):
        m1_M = m1.m
        m1_N = m1.n
        m2_M = m2.m
        m2_N = m2.n

        if (m1_N != m2_M):
            return -1
        else:
            l_gerada = []

            # Lista com elementos das colunas de M2
            colunas = []
            for c2 in range(1, m2_N + 1):
                coluna_atual = []
                for l2 in range (1, m2_M+1):
                    coluna_atual.append(m2.obterElemento(l2, c2))
                colunas.append(coluna_atual)


            for linha_m1 in range (1, m1_M + 1):

                #Lista com elementos da primeira linha de M1
                linha_atual = []
                for n in range (1, m1_N+1):
                    linha_atual.append(m1.obterElemento(linha_m1, n))

                #Percorre cada coluna de M2 já cacheada
                for coluna in colunas:
                    elem = 0
                    for counter, value in enumerate(coluna):
                        elem = elem + (value*linha_atual[counter])
                    l_gerada.append(elem)

            self.definir(m1_M, m2_N, l_gerada)

    def arquivo(self, path):
        '''
        Preenchido linha por linha, da esquerda para direita
        '''
        self.m = 0
        self.n = 0
        self.l = []
        a = 0
        f = open(path, "r")
        for l in f:
            if a == 2:
                self.l.append(float(l))
            elif a == 0:
                self.m = int(l)
                a = a + 1
            elif a == 1:
                self.n = int(l)
                a = a + 1

        return self.mostrarOrdem()

    def mostrarOrdem(self):
        if (self.m == self.n & self.n == 0): return -1
        return [self.m, self.n]

    def numeroElementos(self):
        return self.m * self.n

    def obterElemento(self, a, b):
        return self.l[((a - 1) * self.n + (b - 1))]

    def modificarElemento(self, a, b, c):
        self.l[((a - 1) * self.n + (b - 1))] = float(c)
        return True

    def gerarY(self):
        self.criarY()
        for i in range(1, self.m + 1):
            res = self.b.obterElemento(i, 1)

            for j in range(1, self.n + 1):
                if (i >= j):
                    if (i == j):
                        pass
                    else:
                        res = res - self.y.obterElemento(j, 1) * self.lu.obterElemento(i, j)

                else:
                    pass

            self.y.modificarElemento(i, 1, res)

    def gerarX(self):
        self.criarX()
        for i in range(self.m, 0, -1):
            res = self.y.obterElemento(i, 1)
            div = 1

            for j in range(1, self.n + 1):
                if (i <= j):
                    if (i==j):
                        div = self.lu.obterElemento(i,j)
                    else:
                        res = res - self.x.obterElemento(j, 1) * self.lu.obterElemento(i, j)

            self.x.modificarElemento(i, 1, res/div)

    def mostrarMatriz(self):
        coluna = []
        for i in range (1, self.m+1):
            linha = []
            for j in range (1, self.n+1):
                linha.append(round(self.obterElemento(i,j),2))
            coluna.append(linha)

        for elem in coluna:
            print (elem)
        #print (coluna)

    def mostrarMatriz_noRound(self):
        coluna = []
        for i in range (1, self.m+1):
            linha = []
            for j in range (1, self.n+1):
                linha.append(self.obterElemento(i,j))
            coluna.append(linha)
        for elem in coluna:
            print (elem)
        #print (coluna)

    def determinante(self):
        ordem = self.mostrarOrdem()
        if (ordem == [1, 1]):
            return self.obterElemento(1, 1)
        elif (ordem == [2, 2]):
            return (self.obterElemento(1, 1) * self.obterElemento(2, 2)) - (
                        self.obterElemento(1, 2) * self.obterElemento(2, 1))
        else:
            det = 0
            for j in range(1, self.n + 1):
                det = det + self.obterElemento(1, j) * self.cofator(1, j)
            return det

    def cofator(self, i, j):
        m = Matriz()
        m_l = []
        for contador, valor in enumerate(self.l):
            if (self.pertence(contador, i, j) == False):
                m_l.append(valor)
        m.definir(self.m - 1, self.n - 1, m_l)
        return pow(-1, i + j) * m.determinante()

    def pertence(self, pos, a, b):
        at = (pos - b + 1 + self.n) / self.n
        if (at == int(at)) and (at > 0) and (at <= self.m):
            return True

        bt = (pos - (a - 1) * self.n + 1)
        if (bt == int(bt)) and (bt > 0) and (bt <= self.n):
            return True

        return False

    def matrizLU(self, lu):
        self.lu = lu

    def mostrarLU(self):
        self.lu.mostrarMatriz()

    def matrizB(self, b):
        self.b = b

    def mostrarB(self):
        self.b.mostrarMatriz()

    def criarY(self):
        self.y = Matriz()
        self.y.definir(self.m, 1, 0)

    def criarX(self):
        self.x = Matriz()
        self.x.definir(self.m, 1, 0)

    def mostrarY(self):
        self.y.mostrarMatriz()

    def mostrarX(self):
        self.x.mostrarMatriz()

    def gerarLU(self):
        for k in range(1, self.n - 1 + 1):                
            for i in range(k + 1, self.n + 1):
                self.modificarElemento(i, k, self.obterElemento(i, k) / self.obterElemento(k, k))
                
            self.mostrarMatriz ()

            for j in range(k + 1, self.n + 1):
                for i in range(k + 1, self.n + 1):
                    self.modificarElemento(i, j, self.obterElemento(i, j) - (
                                self.obterElemento(i, k) * self.obterElemento(k, j)))

    def powerMethod(self, l, av, tol):

        if (type(1)==type(l) or len(l)==self.n):
            x = Matriz()
            y = Matriz()

            av = float(av)
            tol = float(tol)

            x.definir(self.n, 1, l)
            r = tol + 1.0

            while (r>tol):
                ava = av*1.0
                y.mult(self, x)
                av = y.obterElemento(1,1)

                for counter, value in enumerate(y.l, 1):
                    x.modificarElemento(counter, 1, value/av)

                r = (abs(av-ava)*1.0)/av

        else:
            return -1

        x.mostrarMatriz()
        print ("Autovalor: " + str(av))
        print ("Tolerância Atingida: " + str(r))

    def add_sub (self, m1, m2, add=True):
        if (m1.n == m2.n and m1.m == m2.m):
            self.m, self.n, self.l = m1.m, m1.n, []

            for a,b in zip(m1.l, m2.l):
                if (add):
                    self.l.append(a+b)
                else:
                    self.l.append(a-b)

        else:
            return -1

    def norma (self, p=2):
        if (self.m>=1 and self.n==1):
            soma = 0.0
            for i in range (1, self.m+1):
                soma = soma + pow(abs(self.obterElemento(i,1))*1.0, p)

            return pow(soma, 1/(p*1.0))

        else:
            return -1

    def simetrica_AntSim(self, sim=True):

        for i in range(1, self.m+1):
            for j in range (1, self.n+1):
                if (i!=j):
                    if (sim):
                        if self.obterElemento(i,j) != self.obterElemento(j,i): return False
                    else:
                        if self.obterElemento(i, j) != -self.obterElemento(j, i): return False
        return True

    def jacobi(self, l, tol):
        if (type(1)==type(l) or len(l)==self.n):
            x0 = Matriz()
            x1 = Matriz()

            #x0.definir(self.n, 1, l)
            x1.definir(self.n, 1, l)

            r = tol+1
            sub = Matriz()
            while(r>tol):
                x0.definir(x1.m, x1.n, x1.l[:])
                for i in range (1, self.n+1):
                    soma_AijXj0 = 0
                    for j in range (1, self.n +1):
                        if (j!=i):
                            soma_AijXj0 = soma_AijXj0 + (self.obterElemento(i,j)*x0.obterElemento(j,1))
                    soma_AijXj0 = self.b.obterElemento(i,1) - soma_AijXj0
                    soma_AijXj0 = soma_AijXj0/self.obterElemento(i,i)
                    x1.modificarElemento(i,1,soma_AijXj0)

                sub.add_sub(x1, x0, False)
                r = sub.norma()*1.0/x1.norma()

            x1.mostrarMatriz()
            print("Tolerância Atingida: " + str(r))
        else:
            return -1

    def jacobi_convergencia_linha(self):

        for i in range (1, self.m+1):
            aii = abs(self.obterElemento(i,i))*1.0

            soma_linha = 0.0
            soma_coluna = 0.0

            for j in range (1, self.n +1):
                if (j!=i):
                    soma_linha = soma_linha + (abs(self.obterElemento(i,j))*1.0)
                    soma_coluna = soma_coluna + (abs(self.obterElemento(j, i)) * 1.0)

            if aii<soma_linha or aii<soma_coluna: return False

        return True
    
    def gaussSeidel(self, l, tol):
        if (type(1)==type(l) or len(l)==self.n):
            x0 = Matriz()
            x1 = Matriz()

            #x0.definir(self.n, 1, l)
            x1.definir(self.n, 1, l)

            r = tol+1
            sub = Matriz()
            while(r>tol):
                x0.definir(x1.m, x1.n, x1.l[:])
                for i in range (1, self.n+1):
                    soma_AijXj0 = 0
                    for j in range (1, self.n +1):
                        if (j!=i):
                            if j == 1:
                                soma_AijXj0 = soma_AijXj0 + (self.obterElemento(i,j)*x0.obterElemento(j,1))
                            else:
                                soma_AijXj0 = soma_AijXj0 + (self.obterElemento(i,j)*x1.obterElemento(j,1))
                    soma_AijXj0 = self.b.obterElemento(i,1) - soma_AijXj0
                    soma_AijXj0 = soma_AijXj0/self.obterElemento(i,i)
                    x1.modificarElemento(i,1,soma_AijXj0)

                sub.add_sub(x1, x0, False)
                r = sub.norma()*1.0/x1.norma()

            x1.mostrarMatriz()
            print("Tolerância Atingida: " + str(r))
        else:
            return -1

    def transposta(self):
        a = Matriz()
        a.definir(self.n, self.m, 0)

        for i in range (1, self.m+1):
            for j in range (1, self.n+1):
                a.modificarElemento(j,i, self.obterElemento(i,j))

        return a

    def pertenceSubMatriz(self, pos, k):
        if (k==self.n): return True
        for m in range(k+1, self.n+1):
            if (self.pertence(pos, m, m)): return False
        return True

    def positivaDefinida(self):

        if (self.m!=self.n or self.n==0): return False

        for k in range (1, self.n-1 +1):
            num_elem = k*k
            num_elem_inseridos = 0
            lm = []
            for counter, value in enumerate(self.l):
                if self.pertenceSubMatriz(counter, k):
                    lm.append(value)
                    num_elem_inseridos = num_elem_inseridos + 1

                if num_elem==num_elem_inseridos: break;

            m = Matriz()
            m.definir(k,k,lm)
            if m.determinante()<=0: return False

        return self.determinante()>0

    def trocaLinha(self, l1, l2):
        if (l1>self.m or l2>self.m): return False

        for n in range(self.n):
            self.l[self.n*(l1-1) + n], self.l[self.n*(l2-1) + n] = self.l[self.n*(l2-1) + n], self.l[self.n*(l1-1) + n]

        return True

    def pivoteamento(self, i, j):

        if (self.obterElemento(i,j)!=0.0): return False
        for m in range (i+1, self.m+1):
            if (self.obterElemento(m,j)!=0.0):
                self.trocaLinha(i,m)
                return [i,m]

        return False

    def eliminacaoGauss(self, fim=True, lu = 0):

        combLinhas = Matriz()
        
        l = Matriz()
        l.definir(self.m, self.m, "I")

        for m in range(1, self.m):
            p = self.pivoteamento(m,m)
            if (type(p)==type([])):
                self.b.trocaLinha(p[0], p[1])

            elem = self.obterElemento(m,m)

            combLinhas.definir(self.m, self.m, "I")

            for m_rest in range(m+1, self.m+1):
                combLinhas.modificarElemento(m_rest, m, -1*self.obterElemento(m_rest,m)/elem)
                l.modificarElemento(m_rest, m, self.obterElemento(m_rest,m)/elem)

            self.mult(combLinhas, self)
            
            if lu == 0:
                self.b.mult(combLinhas, self.b)

        if (fim):
            self.criarX()
            self.substituicaoReversa(self.x, self.b)
            
        if lu != 0:
            return l

    def substituicaoReversa(self, x, b):
        l0 = self.l[:]
        for m in range (self.m, 0, -1):
            i = self.m
            a_analisar = self.m + 1 - m
            c = b.obterElemento(m,1)
            while (a_analisar>0):
                el = l0.pop()
                if a_analisar==1:
                    x.modificarElemento(m,1,c/el)
                    a_analisar = 0
                    for p in range (m-1):
                        l0.pop()
                else:
                    c = c - x.obterElemento(i,1)*el
                    i = i - 1
                    a_analisar = a_analisar - 1

    def eliminacaoGaussJordan(self, fim=True):
        self.eliminacaoGauss(False)

        combLinhas = Matriz()

        for m in range(self.m, 1, -1):
            p = self.pivoteamento(m,m) ###
            if (type(p)==type([])): ###
                self.b.trocaLinha(p[0], p[1]) ###

            elem = self.obterElemento(m, m)

            combLinhas.definir(self.m, self.m, "I")

            for m_rest in range(m-1, 0, -1):
                combLinhas.modificarElemento(m_rest, m, -1 * self.obterElemento(m_rest, m) / elem)

            self.mult(combLinhas, self)
            self.b.mult(combLinhas, self.b)

        if (fim):
            self.criarX()
            self.substituicaoDiagonal(self.x, self.b)

    def substituicaoDiagonal(self, x, b):
        l = []
        for i in range(1,self.m+1):
            l.append(b.obterElemento(i,1)/self.obterElemento(i,i))

        x.definir(self.m, 1, l)

    def copiaMatriz(self, a):
        self.n, self.m, self.l = a.n, a.m, a.l[:]

    def inversa(self, mant=0, quad=1, det=1):
        if (mant!=0): bkp = self.l[:]
        if (quad!=1 and self.m!=self.n): return False
        if (det!=1 and self.determinante()==0): return False

        b = Matriz()
        b.definir(self.m, self.m, "I")
        self.matrizB(b)

        self.eliminacaoGaussJordan(False)

        m_final = Matriz()
        m_final.definir(self.m, self.m, 0)

        for i in range (1, self.m+1):
            m_final.modificarElemento(i,i, 1/self.obterElemento(i,i))

        self.b.mult(m_final, self.b)

        if (mant!=0):
            self.l = bkp
        else:
            self.mult(m_final, self)

        return self.b

    def adicionarLinha(self, elementos):

        try:
            n = self.n
        except:
            n = 0


        if (len(elementos)==n):
            for el in elementos:
                self.l.append(float(el))
            self.m = self.m + 1
            return True
        elif (n==0):
            self.l = []
            for el in elementos:
                self.l.append(float(el))
            self.m = 1
            self.n = len(elementos)
        else:
            return False

    def const_e(self):
        return 2.7183

    def regressaoLinear_coordXY(self):
        som_x, som_y, som_xy, som_x_quad, x_ant = 0,0,0,0,0

        for counter, value in enumerate(self.l):
            if int(counter) % 2 == 0:
                #Par - X
                som_x = som_x + value
                som_x_quad = som_x_quad + pow(value,2)
                x_ant = value
            else:
                #Ímpar - Y
                som_y = som_y + value
                som_xy = som_xy + (x_ant*value)

        a_r = Matriz(2,2,[self.m, som_x, som_x, som_x_quad])
        c_r = Matriz(2,1,[som_y, som_xy])

        b_r = Matriz()
        b_r.mult(a_r.inversa(), c_r)
        self.matrizB(b_r)
        return self.b

    def substituicaoFrente(self, y, b):

        #SELF*Y = B

        div = 0

        for i in range(1, self.m + 1):
            res = b.obterElemento(i, 1)

            for j in range(1, self.n + 1):
                if (i >= j):
                    if (i==j):
                        div = self.obterElemento(i,j)
                    else:
                        res = res - y.obterElemento(j, 1) * self.obterElemento(i, j)

                else:
                    pass

            y.modificarElemento(i, 1, res/div)

    def decomposicao_lu (self, b):
        self.matrizB(b)
        t = self.eliminacaoGauss(lu = 1)
        self.mostrarMatriz()
        
        y = Matriz(self.m, 1, 0)
        t.substituicaoFrente (y, b)
        y.mostrarMatriz ()
        
        x = Matriz(self.m, 1, 0)
        self.substituicaoReversa (x, y)
        
        x.mostrarMatriz ()
    
    def decomposicaoCholesky2 (self, b):
        l_cho = Matriz(self.m, self.m, 0)

        for i in range(1, self.m +1):
            for k in range(1, i + 1):
                soma_tmp = sum(l_cho.obterElemento(i,j) * l_cho.obterElemento(k,j) for j in range(k))

                if (i == k):
                    l_cho.modificarElemento(i,k, sqrt(self.obterElemento(i,i) - soma_tmp))

                else:

                    l_cho.modificarElemento(i,k, (1.0 / l_cho.obterElemento(k,k)* (self.obterElemento(i,k) - soma_tmp)))

        u_cho = Matriz()
        u_cho.copiaMatriz(l_cho.transposta())

        #print("Matriz L")
        #l_cho.mostrarMatriz()
        #print("\n")
        #print("Matriz U")
        #u_cho.mostrarMatriz()
        #print("\n")

        y = Matriz(self.m, 1, 0)
        l_cho.substituicaoFrente(y, b)
        #y.mostrarMatriz()

        #print("\n")

        x = Matriz(self.m, 1, 0)
        u_cho.substituicaoReversa(x, y)
        #x.mostrarMatriz()

        return x
        
        #Decomposição de Cholesky 1
        #a = np.matrix ("5, -4, 1, 0; -4, 6, -4, 1; 1, -4, 6, -4; 0, 1, -4, 5")
        #a = np.matrix ("16 9 8 7 6 5 4 3 2 1; 9 17 9 8 7 6 5 4 3 2; 8 9 18 9 8 7 6 5 4 3; 7 8 9 19 9 8 7 6 5 4; 6 7 8 9 18 9 8 7 6 5; 5 6 7 8 9 17 9 8 7 6; 4 5 6 7 8 9 16 9 8 7; 3 4 5 6 7 8 9 15 9 8; 2 3 4 5 6 7 8 9 14 9; 1 2 3 4 5 6 7 8 9 13")
        #b = np.matrix ("-1; 2; 1; 3")
        #b = np.matrix ("4; 0; 8; 0; 12; 0; 8; 0; 4; 0")
        '''
        dim_a = len (a)
        
        a_cho_tmp = np.array (a.copy ())
        a_cho = []
        
        for i in a_cho_tmp:
            a_cho += [list (i)]
        
        l_cho_tmp = np.zeros ((dim_a, dim_a))
        l_cho = []
        
        for i in l_cho_tmp:
            l_cho += [list (i)]
        
        for i in range (dim_a):
            for k in range (i + 1):
                soma_tmp = sum (l_cho [i][j] * l_cho [k][j] for j in range (k))
                
                if (i == k):
                    l_cho [i][k] = sqrt (a_cho [i][i] - soma_tmp)
                    
                else:
                    l_cho [i][k] = (1.0 / l_cho [k][k] * (a_cho [i][k] - soma_tmp))
        
        u_cho = np.matrix.transpose (np.array (l_cho)) # l_cho é uma lista de listas!
            
        print ("Matriz L") 
        print (np.array (l_cho))
        print ("\n")
        print ("Matriz U")
        print (u_cho)
        print ("\n")
        
        y = np.linalg.solve (l_cho, b)
        x = np.linalg.solve (u_cho, y)
        print (x)
        '''

    def idMaiorElem_teste(self, modulo=True):
        '''
        Retorna o [elemento, posicao_l, linha, coluna] do maior elemento em módulo fora da diagonal
        :return:
        '''
        res = -1
        for counter, value in enumerate(self.l):
            ij = self.ij_posicao(counter)
            if (ij[0]!=ij[1]):
                if (type(res)==type([])):
                    if (modulo):
                        if abs(value)>=abs(res[0]):
                            res = [value, counter, ij[0], ij[1]]

                    else:
                        if value >= res[0]:
                            res = [value, counter, ij[0], ij[1]]

                else:
                    res = [value, counter, ij[0], ij[1]]

        return res

    def idMaiorElem(self, modulo=True):
        '''
        PERCORRE POR COLUNA DE CIMA PRA BAIXO, DA ESQUERDA PARA DIREITA
        Retorna o [elemento, posicao_l, linha, coluna] do maior elemento em módulo fora da diagonal
        :return:
        '''
        res = -1
        for j in range(1,self.n+1):
            for i in range (1, self.m+1):
                value = self.obterElemento(i,j)
                counter = -1
                ij = [i,j]
                if (ij[0]!=ij[1]):
                    if (type(res)==type([])):
                        if (modulo):
                            #print(value, res[0], value==res[0], abs(value)==abs(res[0]))
                            if abs(value)>=abs(res[0]):
                                res = [value, counter, ij[0], ij[1]]
                            #elif abs(value)==abs(res[0]) and ij[0]==ij[1]:
                            #    res = [value, counter, ij[0], ij[1]]

                        else:
                            if value >= res[0]:
                                res = [value, counter, ij[0], ij[1]]
                            #elif value==res[0] and ij[0]==ij[1]:
                            #    res = [value, counter, ij[0], ij[1]]
                    else:
                        res = [value, counter, ij[0], ij[1]]

        return res

    def jacobi_auto(self, tol):
        a = Matriz()
        a.copiaMatriz(self)

        x = Matriz(self.m, self.m, "I")

        k = 1
        while (k!=0):
            elemMaior = a.idMaiorElem()

            aii, ajj, aij = a.obterElemento(elemMaior[2], elemMaior[2]), a.obterElemento(elemMaior[3], elemMaior[3]), a.obterElemento(elemMaior[2], elemMaior[3])

            if (aii==ajj):
                fi = pi/4
                cos_fi = round(cos(fi),5)
                sen_fi = round(sin(fi),5)
                n_sen_fi = (-1)*sen_fi
            else:
                fi = 0.5 * atan(2.0 * aij/(aii-ajj))
                cos_fi = round(cos(fi),5)
                sen_fi = round(sin(fi),5)
                n_sen_fi = (-1) * sen_fi



            p = Matriz(self.m, self.m, "I")
            p.modificarElemento(elemMaior[2], elemMaior[2], cos_fi)
            p.modificarElemento(elemMaior[3], elemMaior[3], cos_fi)
            p.modificarElemento(elemMaior[2], elemMaior[3], n_sen_fi)
            p.modificarElemento(elemMaior[3], elemMaior[2], sen_fi)

            a.mult(a, p)
            x.mult(x, p)
            a.mult(p.transposta(), a)

            #print("K = " + str(k))
            #print ("Termo a ser zerado: (" + str(elemMaior[2]) + ", " + str(elemMaior[3]) + ")")
            #print("FI: " + str(fi))
            #print ("Matriz P"); p.mostrarMatriz_noRound(); print ("---------");
            #print ("Matriz A"); a.mostrarMatriz_noRound(); print ("---------");
            #print ("Matriz X"); x.mostrarMatriz_noRound(); print ("---------");
            #print ("FINALIZADO K\n---------")
            #Checar Convergência
            if (abs(a.idMaiorElem()[0])<=tol):
                k=0
            else:
                k = k + 1

        #print ("FIM PROCESSO\n---------------")

        #print("MATRIZ AUTOVALORES")
        #a.mostrarMatriz()
        #print("MATRIZ AUTOVETORES")
        #x.mostrarMatriz()


        for counter, value in enumerate(a.l):
            ij = a.ij_posicao(counter)
            if (ij[0]!=ij[1]):
                a.l[counter] = 0


        return [a,x]
        

###Lista 1
        
##1
#a = Matriz(4,4,[5,-4,1,0,-4,6,-4,1,1,-4,6,-4,0,1,-4,5])
#b = Matriz(4,1,[-1,2,1,3])
#a = Matriz(10,10,[16,9,8,7,6,5,4,3,2,1,9,17,9,8,7,6,5,4,3,2,8,9,18,9,8,7,6,5,4,3,7,8,9,19,9,8,7,6,5,4,6,7,8,9,18,9,8,7,6,5,5,6,7,8,9,17,9,8,7,6,4,5,6,7,8,9,16,9,8,7,3,4,5,6,7,8,9,15,9,8,2,3,4,5,6,7,8,9,14,9,1,2,3,4,5,6,7,8,9,13])
#b = Matriz(10,1,[-1,2,1,3])

#VERIFICAÇÃO (positiva definida)
'''
print(a.positivaDefinida())
'''

#VERIFICAÇÃO - Cholesky (simétrica)
'''
print(a.simetrica_AntSim())
'''

#Decomposição LU + solução AX=B
'''
lu = Matriz()
lu.definir(a.m, a.n, a.l[:])
lu.gerarLU()
a.matrizLU(lu)
a.mostrarLU()
a.matrizB(b)
a.gerarY()
a.gerarX()
a.mostrarX()
'''

#Decomposição de Cholesky
'''
x = a.decomposicaoCholesky2(b)
x.mostrarMatriz()
'''


###Lista 2
#a = Matriz(3,3,[3,2,0,2,3,-1,0,-1,3])
#b = Matriz(3,1,[1,-1,1])

##1
#VERIFICAÇÃO (diagonal dominante)
'''
print(a.jacobi_convergencia_linha())
'''

#Power Method
'''
a.powerMethod([1,1,1], 1, pow(10, -3))
'''

##2
#VERIFICAÇÃO (simétrica)
'''
print(a.simetrica_AntSim())
'''

#Jacobi para autovalores e autovetores

'''a = Matriz(3,3,[1,0.2, 0, 0.2, 1, 0.5, 0, 0.5, 1])
res = a.jacobi_auto(0.00001)
res[0].mostrarMatriz_noRound()
print()
res[1].mostrarMatriz_noRound()'''


##3
#VERIFICAÇÃO - Jacobi (diagonal dominante)
'''
print(a.jacobi_convergencia_linha())
'''

#VERIFICAÇÃO - Gauss-Seidel (positiva definida)
'''
print(a.positivaDefinida())
'''

#VERIFICAÇÃO - Gauss-Seidel (simétrica)
'''
a = Matriz(3,3,[1,2,4,2,3,1,4,1,2])
print(a.simetrica_AntSim())
'''

#Jacobi e Gauss-Seidel
'''
a.matrizB(b)
usuario = input ("Jacobi (0) ou Gauss-Seidel (1): ")
if int (usuario) == 0:
    a.jacobi(1, 0.00001)
else:
    a.gaussSeidel(1, 0.00001)
'''

##4 (ver Lista 1 - 1)


###Lista 3

##9
'''
l = []
iteracao, parar = 1, 1

while int (parar) == 1:
    xi = input ("Coordenada x" + str (iteracao) + ": ")
    yi = input ("Coordenada y" + str (iteracao) + ": ")
    
    l.append (int (xi))
    l.append (int (yi))
    
    iteracao += 1
    
    parar = input ("Inserir mais coordenadas? (1 = SIM): ")

a = Matriz(3,2,l)
a.regressaoLinear_coordXY().mostrarMatriz()
'''


###Caso com pivoteamento
#a = Matriz (3, 3, [0, 1, 1, 1, 2, 1, 1, 1, -1])
#b = Matriz (3, 1, [4, 7, 3])

##Eliminação de Gauss
'''
a.matrizB(b)
a.eliminacaoGauss()
a.mostrarMatriz()
print('\n')
a.b.mostrarMatriz()
print('\n')
a.mostrarX()
'''

##Decomposição LU
'''
a.decomposicao_lu (b)
'''

##Eliminação de Gauss-Jordan
'''
a.matrizB(b)
a.eliminacaoGaussJordan()
a.mostrarX()
'''