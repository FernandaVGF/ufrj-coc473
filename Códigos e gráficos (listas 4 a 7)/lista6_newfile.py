import math

def generate_csv(filename, l, ti=0):
	a = ""
	c = ""
	cont = 0
	t = ti
	for i in l:
		if cont!=0:
			a = a + ";" + str(i)
			c = c + ";" + str(t)
			t = t + 1
		else:
			a = str(i)
			c = str(t)
			t = t + 1
			cont=1

	file = open(filename, 'w')
	file.write(c)
	file.write(a)
	file.close()

def rungekutta_1_euler(ki, kf, h, xi, func):
	x=[xi]
	for i in range(ki+1,kf):
		tk = i*h
		xk = x[-1]
		k1 = func(tk, xk)
		x.append(xk + k1*h)
	return x

def rungekutta_2(ki, kf, h, xi, func):
	x=[xi]
	for i in range(ki,kf):
		tk = i*h
		xk = x[-1]
		k1 = func(tk, xk)
		k2 = func(tk+h, xk+h*k1)
		x.append(xk+(k1+k2)*1.0*(h/2))
	return x

def rungekuta_4(ki, kf, h, xi, func):
	x=[xi]
	for i in range(ki,kf):
		tk = i*h
		xk = x[-1]
		k1 = func(tk, xk)
		k2 = func(tk+1.0*(h/2), xk+1.0*(h/2)*k1)
		k3 = func(tk+1.0*(h/2), xk+1.0*(h/2)*k2)
		k4 = func(tk+h, xk+h*k3)
		x.append(xk+(k1+2*k2+2*k3+k4)*1.0*(h/6))
	return x

def edo_serietaylor(t0, x0, x0_l, ki, kf, h, func):
	x=[x0]
	x_l = [x0_l]
	x_ll = []
	for i in range(ki+1,kf+1):
		tk = t0 + 1.0*h*(i-1)
		x_ll.append(func(tk, x[-1], x_l[-1]))
		x.append(x[-1] + x_l[-1]*h + (x_ll[-1]/2)*h*h)
		x_l.append(x_l[-1] + x_ll[-1]*h)
	return x

def rungekuta_nystrom(t0, x0, x0_l, ki, kf, h, func):
	x = [x0]
	x_l = [x0_l]
	x_ll = []
	meio_h = 0.5*h
	for i in range(ki,kf):
		tk = i*h + t0
		k1 = meio_h*func(tk, x[-1], x_l[-1])
		q = meio_h * (x_l[-1] + 0.5*k1)
		k2 = meio_h*func(tk + meio_h, x[-1] + q, x_l[-1] + k1)
		k3 = meio_h*func(tk + meio_h, x[-1] + q, x_l[-1] + k2)
		l = h * (x_l[-1] + k3)
		k4 = meio_h * func(tk + h, x[-1] + l, x_l[-1] + 2*k3)
		x.append(x[-1] + 1.0*h*(x_l[-1] + (1.0/3)*(k1+k2+k3)))
		x_l.append(x_l[-1] + (1.0/3)*(k1+2*k2+2*k3+k4))
	return x

def exemplo1(t,x):
	#x' = t + x
	return 1.0*t + x

def exemplo2(t, z, z_l):
	return (-1.0)*9.80665-1.0*z_l*abs(z_l)
	
def ex1(t,y):
	return (-2.0)*t*y*y

def questao6(t,y):
	return (-10.0)*t*y*y

def questao9(t, y, y_l):
	m=1.0
	c=0.2
	k=1.0
	w=0.5
	f = 2 + math.cos(3.0*w*t)
	return (f - k*y - c*y_l)/m
	

#print(rungekutta_1_euler(0,3,0.1,0.0,exemplo1))
#print(rungekutta_2(0,2,0.1,0.0,exemplo1))
#print(rungekutta_4(0,2,0.1,0.0,exemplo1))
#print(edo_serietaylor(0, 0.0, 0.0, 0, 10, 0.1, exemplo2))
#rungekuta_nystrom(t0, x0, x0_l, ki, kf, h, func):
#print(rungekuta_nystrom(0, 0.0, 0.0, 0, 10, 0.1, exemplo2))

#l = rungekuta_nystrom(0, 0.0, 0.0, 0, 10, 0.1, exemplo2)
#generate_csv("nystromteste.csv", l)

#Questao 6
for tempo in range(11): 
    print(tempo/10.0)
q6 = rungekutta_2(0,11,0.1,2.0,questao6)
print(q6)
generate_csv("questao6_rk2.csv", q6)
q6_b = rungekutta_1_euler(0,2, 0.1, 2, questao6)
#print(q6_b)
#generate_csv("questao6_rk1_euler.csv", q6_b)

#rungekutta_1_euler(ki, kf, h, xi, func

#Questao 9
l_h = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
for i in l_h:
    q9 = rungekuta_nystrom(0, 0.0, 0.0, 0, 100, i, questao9)
    nome = "questao9_" + str(i) + ".csv"
    #generate_csv(nome, q9)
    a = ""; c = ""; cont = 0; t = 0
    for i in q9:
        if cont!=0:
            a = a + ";" + str(i)
            c = c + ";" + str(t)
            t = t + 1
        else:
            a = str(i)
            c = str(t)
            t = t + 1
            cont=1
    
    file = open(nome, 'w')
    file.write(c)
    file.write(a)
    file.close()

#print(q9)

q9 = rungekuta_nystrom(0, 0.0, 0.0, 0, 400, 0.25, questao9)
file = open('questao9_0.25.txt', 'w')
for i in q9:
    file.write(str(i)+"\n")
file.close()