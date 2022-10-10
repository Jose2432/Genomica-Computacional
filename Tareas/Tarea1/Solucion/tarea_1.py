__author__ = "Guzmán Cortés Dulce Dyliang"
__author__ = "Marin Parra Jose Guadalupe de Jesus"
__author__ = "Rangel Limón Erik"

#!coding:utf8
import math
import random
import re

#Ejercicio 2
cero   = 'ATATATACATACTGGTAATGGGCGCGCGTGTGTTAAGTTCTGTTGTAGGGGTGATTAGGGGCG'#
uno    = 'GGCCCACACCCCACACCAATATATGTGGTGTGGGCTCCACTCTCTCGCGCTCGCGCTGGGGAT'
dos    = 'ATAAGGTGTGTGGGCGCGCCCCGCGCGCGCGTTTTTTCGCGCGCCCCCGCGCGCGCGCGCGCG'
tres   = 'GGCGCGGGACGCGGCGGCGGATCCCGATCCGTGCGTCAATACTATTATGGCCAGATAGAATAA'#
cuatro = 'GTGCTGCTGCGGCGCCCACACCTATTATCTCTCTCTCTCTGCCTCTCCACCTCGGGGCTTAAT'
cinco  = 'GCGCTGCTGCTGGCTCGATGGGCGCGTGCGTCGTAGCTCGATGCTGGCTCGAGCTGTAATCTT'#
seis   = 'GGCGCTCGCTCGGATGCGCGGCCGGGCTCTCTGCTCGCGCTCGCTTCGCGCTCGTGACCGCTG'
siete  = 'AATTGGTGCGCGCTCGCGCACACACAGAGAGAGGGTTTATATAGGATGATATATCCACATTGG'
ocho   = 'ATGCTGCTGCTGGCTCTGCTTGCGCTCTGCTCGCTGGGGTGTGTGTGCCGCGCGCTGCTGCTC'
nueve  = 'GCTGGGCTCGCTCGATGCGCGCGGGCGCGCGACCGCGGACGGCGTCGCTGCTAAATGGGCTTC'
lista_secuencias = [cero,uno,dos,tres,cuatro,cinco,seis,siete,ocho,nueve]
def ejercicio_2(lista_secuencias):
    L = []
    r = re.compile('(ATG|TTG|GTG)([ACTG]{3})+(TAG)')
    L = [i for i, item in enumerate(lista_secuencias) if re.search(r, item)]
    return L

#Auxiliar para el ejercicio_3 que se encarga de leer un archivo
def lee_archivo(nombre_archivo):
    with open(nombre_archivo, "r") as file:
        for line in file:
            yield line.strip()

promotores = ["AGATAG", "TGATAG", "AGATAA", "TGATAA"]

#Ejercicio 3
def ejercicio_3(archivo):
    L = []
    for line in lee_archivo(archivo):
        cuenta = 0
        line = re.split(r'\s+', line)[1]
        for _ in re.finditer('|'.join(promotores), line):
            cuenta += 1
        L.append(cuenta)
    return L

#Ejercicio 4
def ejercicio_4(M):
    x = 0.
    i = 1
    D = 0
    while i < M:
        i = i+1
        x = random.uniform(-1,1)
        y = random.uniform(-1,1)
        d = math.sqrt((x**2)+(y**2))
        if d <= 1:
            D = D+1
    x = 4*D/i
    return x

#------------------------------------------ M E N U ------------------------------------------

#Ejercicio 1 (c)
cadena = 'GATTATATACATAGTAGTATA'
er = 'TA'
print('\nEJERCICIO 1.c \nDada la cadena',cadena)
for match in re.finditer(er, cadena):
    s = match.start()
    e = match.end()
    print('Extiste un "%s" en %d:%d' % (cadena[s:e], s, e))

print('\nEJERCICIO 2 \nLas lineas con genes validos son:',ejercicio_2(lista_secuencias),"\n")

print('EJERCICIO 3 \nPara el archivo promotores.txt, tenemos lo siguiente.\n',ejercicio_3('promotores.txt'),'\n')

while True:
  try:
     inputM = int(input("EJERCICIO 4 \nIngresa un numero entero.\n"))
  except ValueError:
     print("La cadena que ingresaste no es un numero o lo escribiste mal.\n")
     break
  else:
     print("El numero ingresado es " + str(inputM))
     print("La solucion del ejercicio 4 con el numero " + str(inputM) + " ingresado es:", ejercicio_4(inputM),"\n")
     break
