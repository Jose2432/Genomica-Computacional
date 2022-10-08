#!coding:utf8
import math
import random
import re

__author__ = "Marin Parra Jose Guadalupe de Jesus"
__NoCuenta__ = "316264176"

#Ejercicio 1 (c)
#cadena = 'GATTATATACATAGTAGTATA'
#er = 'TA'
#print('Dada la cadena',cadena)
#for match in re.finditer(er, cadena):
#    s = match.start()
#    e = match.end()
#    print('Extiste un "%s" en %d:%d' % (cadena[s:e], s, e))

#Inicio ATG,TTG,GTG
#Final TAG
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
    r = re.compile('[ATG,TTG,GTG].*(TAG)') #Agregar los tripletes de nucleotidos
    L = [i for i, item in enumerate(lista_secuencias) if re.search(r, item)]
    return L

print('Ejercicio 2:',ejercicio_2(lista_secuencias))

def ejercicio_3(archivo):
    L = []
    ...
    return L

#Agregar M como input
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

#print("Ejercicio 4:",ejercicio_4(3))
