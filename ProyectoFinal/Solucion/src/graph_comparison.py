"""

Este script genera el árbol filogenético dada una carpeta con archivos
FASTA.

Jose Guadalupe de Jesús Marín Parra
Erik Rangel Limón

"""
# Paquetería usada para leer archivos
import os
# Paquetería para el uso de REGEX
import re

# Paqueterías para crear árboles binarios.
# Es necesario descargarla con:
# $ pip install binarytree
from binarytree import Node
from binarytree import build
# Para generar la imagen con el árbol filogenético es necesario tener
# instalado el programa "graphviz"
# Fedora:
# $ sudo dnf install graphviz
# Debian (y derivados):
# $ sudo apt install graphviz
# Windows:
# $ winget install graphviz
# o bien, también en Windows pero con chocolatey:
# $ choco install graphviz

"""
Esta clase modela el comportamiento de la Matriz de Adyacencias
reducida de una gráfica bipartita completa con pesos en sus aristas.
"""
class WCompleteBipartite:

    """
    Método constructor que recibe dos listas y genera una matriz de
    adyacencias en 0's
    """    
    def __init__(self,v1,v2):
        self.v1 = v1
        self.indexv1 = {}
        for i in range(len(v1)):
            self.indexv1.update({v1[i] : i})
        self.v2 = v2
        self.indexv2 = {}
        for i in range(len(v2)):
            self.indexv2.update({v2[i] : i})
        self.matrix = [0] * len(v1)
        for i in range(len(v1)):
            self.matrix[i] = [0] * len(v2)

    """
    Método que recibe dós vértices de la gráfica y suma uno a su
    respectiva arista en la matriz de adyacencias.
    """            
    def update(self,u,v):
        if self.indexv1.get(u) != None:
            if self.indexv2.get(v) != None:
                w = self.matrix[self.indexv1[u]][self.indexv2[v]]
                self.matrix[self.indexv1[u]][self.indexv2[v]] = w + 1
                return
            elif self.indexv1.get(v) != None:
                raise Exception("Two elements are in the same independent set")
            else:
                raise Exception(f"No such element {v}")
        elif self.indexv1.get(v) != None:
            if self.indexv2.get(u) != None:
                w = self.matrix[self.indexv1[v]][self.indexv2[u]]
                self.matrix[self.indexv1[v]][self.indexv2[u]] = w + 1
                return
            elif self.indexv1.get(u) != None:
                raise Exception("Two elements are in the same independent set")
            else:
                raise Exception(f"No such element {v}")
        raise Exception("No such elements")

"""
Clase usada para asignarle a cada archivo dado su respectiva matriz de
adyacencias reducida de una gráfica bipartita.
"""
class DataGen:

    """
    Constructor sin parámetros con los respectivos posibles
    dinucluetoidos y nucleotidos.
    """
    def __init__(self):
        self.dinucleotide = ['AA','AC','AT','AG','CA','CC','CT','CG','TA','TC','TT','TG','GA','GC','GT','GG']
        self.nucleotide = ['A','C','T','G']
        self.data = {}

    """
    Constructor con parámetros que crea la asignación de cada gráfica
    para cada archivo dado en el folder especificado.
    """
    def __init__(self,folder):
        self.dinucleotide = ['AA','AC','AT','AG','CA','CC','CT','CG','TA','TC','TT','TG','GA','GC','GT','GG']
        self.nucleotide = ['A','C','T','G']
        self.data = {}
        for filename in os.scandir(folder):
            if filename.is_file():
                self.scanfile(filename.path)

    """
    Método que lee la secuencia y va actualizando los parámetros de la
    gráfica bipartita correspondiente.
    """                
    def scansecuence(self, name, secuence):
        bipartite_graph = WCompleteBipartite(self.dinucleotide,self.nucleotide)
        for i in range(len(secuence)):
            if (i+2 >= len(secuence) - 1):
                break
            if secuence[i] not in self.nucleotide or secuence[i+1] not in self.nucleotide or secuence[i+2] not in self.nucleotide:
                continue
            u = secuence[i:i+2]
            v = secuence[i+2]
            bipartite_graph.update(u,v)
        if name in self.data:
            count = 0
            for k in self.data:
                if re.split(r'\s', k)[0] == name:
                    count = count + 1
            name = name + f' ({count - 1})'
        self.data.update({name : bipartite_graph})


    """
    Hace lo mismo que el anterior, pero con un archivo 
    """
    def scanfile(self,filename):
        print(f'Scanning {filename}...')
        with open(filename, 'r') as file:
            lines = file.readlines()
            firstline = re.split(r'\s', lines[0])
            name = firstline[0][1:] + ' ' + firstline[1] + r'\n' + re.split(r'\.', re.split(r'/', filename)[2])[0]
            secuence = lines[1]
            self.scansecuence(name, secuence)


"""
Clase que crea la matriz de distancias a partir de los datos generados
de la clase DataGen.
"""
class DistanceMatrix:

    """
    Método constructor que recibe como parámetro los datos de un
    DataGen y llena la matriz de distancias.
    """
    def __init__(self, data):
        print('Generating distance matrix')
        self.data = data
        self.files = list(data.keys())
        self.distance = [0] * len(self.files)
        for i in range(len(self.files)):
            self.distance[i] = [0] * len(self.files)
        self.fillDistance()

    """
    Método que llena cada elemento en la matriz de distancias con la
    operación que nos regresa el siguiente método.
    """
    def fillDistance(self):
        for i in range(len(self.files)):
            for j in range(len(self.files)):
                if (i == j):
                    continue
                else:
                    self.distance[i][j] = self.distanceOperation(self.files[i], self.files[j])

    """
    Método que dadas dos secuencias, nos regresa la distancia a la que
    se encuentran. Ésta distancia es el promedio aritmético de la
    diferencia de apariciones de cada trinucleótido entre las dos
    secuencias. En total son 64 posibles combinaciones.
    """
    def distanceOperation(self, x, y):
        x_values = self.data[x].matrix
        y_values = self.data[y].matrix
        addition = 0
        for i in range(len(x_values)):
            for j in range(len(x_values[i])):
                addition = addition + abs(x_values[i][j]-y_values[i][j])
        return (addition / 64)

"""
Clase para crear el árbol filogenético por medio del algoritmo UPGMA
con la matriz de distancias generada en la clase anterior.
"""
class UpgmaMatrix:

    """
    Método constructor que recibe las secuencias y su respectiva
    matriz de distancias.
    """
    def __init__(self, files, distance):
        self.leaves = files
        self.distance = distance

    """
    Método auxiliar para contar las hojas de un subárbol.
    """
    def countLeaves(self, node):
        if node == None:
            return 0
        if node.left == None and node.right == None:
            return 1
        if node.left == None:
            return self.countLeaves(node.right)
        if node.right == None:
            return self.countLeaves(node.left)
        return self.countLeaves(node.left) + self.countLeaves(node.right)

    """
    Método auxiliar para calcular la distancia del nuevo clúster
    (A U B) con otro cluster X.
    """
    def distanceOperation(self, card_a, card_b, d_ax, d_bx):
        numerator = card_a * d_ax + card_b * d_bx
        denominator = card_a + card_b
        return numerator / denominator

    """
    Método que nos regresa el árbol filogenético.
    Especificación: https://en.wikipedia.org/wiki/UPGMA    
    """
    def getFilTree(self):
        print('Creating phylogenetic tree...')
        fil_leaves = []
        for i in range(len(self.leaves)):
            fil_leaves.append(Node(self.leaves[i]))
        tree_distance = self.distance.copy()
        while (True):
            leaf_index = {}
            for i in range(len(fil_leaves)):
                leaf_index.update({fil_leaves[i] : i})
            minvalue = -1
            minleaf_l = -1
            minleaf_r = -1
            for i in range(len(tree_distance)):
                for j in range(len(tree_distance[i])):
                    if i == j:
                        continue
                    if minvalue == -1:
                        minvalue = tree_distance[i][j]
                        minleaf_l = i
                        minleaf_r = j
                        continue
                    if tree_distance[i][j] < minvalue:
                        minvalue = tree_distance[i][j]
                        minleaf_l = i
                        minleaf_r = j
                        continue
            node_l = fil_leaves[minleaf_l]
            node_r = fil_leaves[minleaf_r]
            fil_leaves.remove(node_l)
            fil_leaves.remove(node_r)
            subtree = Node(round(minvalue / 2,2))
            subtree.left = node_l
            subtree.right = node_r
            fil_leaves.insert(0, subtree)
            if (len(fil_leaves) == 1):
                break
            new_distance = [0] * len(fil_leaves)
            for i in range(len(new_distance)):
                new_distance[i] = [0] * len(fil_leaves)
            card_a = self.countLeaves(node_l)
            card_b = self.countLeaves(node_r)
            for i in range(1,len(fil_leaves)):
                d_ax = tree_distance[minleaf_l][leaf_index[fil_leaves[i]]]
                d_bx = tree_distance[minleaf_r][leaf_index[fil_leaves[i]]]
                dist = self.distanceOperation(card_a,card_b,d_ax,d_bx)
                new_distance[0][i] = dist
                new_distance[i][0] = dist
            for i in range(1,len(fil_leaves)):
                for j in range(i+1,len(fil_leaves)):
                    val = tree_distance[leaf_index[fil_leaves[i]]][leaf_index[fil_leaves[j]]]
                    new_distance[i][j] = val
                    new_distance[j][i] = val
            tree_distance = new_distance
        print('¡Phylogenetic tree generated!')
        return fil_leaves[0]

"""
Método general para crear el diagrama del árbol filogenético de todos
los archivos dentro de una carpeta.
"""    
def prueba(carpeta):
    datagen = DataGen(carpeta)
    distance = DistanceMatrix(datagen.data)
    upgma = UpgmaMatrix(distance.files, distance.distance)
    fil_tree = upgma.getFilTree()
    graph = fil_tree.graphviz()
    graph.body
    graph.render(filename='output/' + re.split(r'/', carpeta)[1])

"""
Resultados aquí:
"""
prueba("input/mamiferos")
prueba("input/bacteria")
