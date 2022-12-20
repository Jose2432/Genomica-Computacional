__author__ = "Marin Parra Jose Guadalupe de Jesus"
__NoCuenta__ = "316264176"

import os.path as path

#Lee la secuancia del archivo FASTA
def leerDNA(nombreArchivo):
    cadena = ''
    with open(nombreArchivo, 'r') as f:
        for line in f:
            # ignore header line with cadena information
            if not line[0] == '>':
                cadena += line.rstrip()
    return cadena

#Transcribe la secuencia FASTA
def transcribir(nombreArchivo):
    cadena = ''
    with open(nombreArchivo, 'r') as f:
        for line in f:
            if not line[0] == '>':
                cadena += line.replace('T','U')
    return cadena

#Traduce la el RNA mensajero a proteina
def traducir(secuencia):
    codones_traduccion = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}
    proteina = ""
    if len(secuencia)%3 == 0:
        for i in range(0, len(secuencia), 3):
            codon = secuencia[i : i+3]
            proteina += codones_traduccion[codon]
    return proteina

#Menu
archivoConsultar = raw_input("Asegurate de que el archivo que quieras consultar este en el mismo directorio que este archivo.py\n"
+ "Escribe el nombre del archivo que quieres revisar, incluida la extension(nombreArchivo.extension).\n")
if path.exists(archivoConsultar):
    print("El archivo a consultar es " + str(archivoConsultar))
    cadena = leerDNA(archivoConsultar)
    print("\nLa secuencia de DNA es la siguiente.\n" + str(cadena))
    rnam = transcribir(archivoConsultar)
    print("\nLa secuencia de mRNA es el siguiente.\n" + str(rnam))
    print('\nLa secuencia de aminoacidos es el siguiente.\n'+str(traducir(leerDNA(archivoConsultar)))+"\n")
else:
    print("El nombre del archivo que ingresaste no existe o lo escribiste mal.")
