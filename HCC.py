import sys
import matplotlib.pyplot as plt

#Definición de funciones--------------------------------------------------------

# Funcion que abre el archivo con formato fasta, ignora la linea si empieza con ">"
# lee cada linea como string y la convierte a lista de caracteres en mayuscula
# luego combina todas las listas en una

def leer_seq(imputfile):
    with open(imputfile) as f:
        secuencia = []
        for line in f:
            if line.startswith( ">"):
                continue
            a=line.replace("\n","")
            b=list(a.upper())
            secuencia.extend(b)
    return secuencia

# Funcion que me permite buscar los marcos abiertos de lectura

def buscar_orf (adn):
    list_marcos = []
    marco = ""
    estado = 0
    i=0
    while i < len(adn):
        if estado == 0:
            if adn[i] == "A":
                estado = 1 
                marco += adn[i]
                inicio = i+1
        elif estado == 1:
            if adn[i] == "T":
                estado = 2
                marco += adn[i]
            elif adn[i] == "A":
                estado = 1 
                inicio = i+1
            else:
                estado = 0
                marco = ""
        elif estado == 2:
            if adn[i] == "G":
                estado = 3
                marco += adn[i]
            elif adn[i] == "A":
                estado = 1
                marco = ""
                marco += adn[i]
                inicio = i+1
            else:
                estado = 0
                marco = ""
        elif estado == 3:
            if adn[i] == "T" and (i+1-inicio) %3== 0:
                estado = 4
                marco += adn[i]
            else:
                estado = 3
                marco += adn[i]
        elif estado == 4:
            if adn[i] == "A":
                estado = 5
                marco += adn[i]
            elif adn[i] == "G":
                estado = 6
                marco += adn[i]
            else:
                estado = 3
                marco += adn[i]
        elif estado == 5:
            if adn[i] == "G" or adn[i] == "A":
                marco += adn[i]
                if len(marco) >= 30:
                    fin = i+1
                    list_marcos.append([marco , inicio , fin])
                    marco=""
                    estado = 0
                    i=inicio+1
                else:
                    marco = ""
                    estado = 0
            else:
                estado = 3
                marco += adn[i]
        elif estado == 6:
            if adn[i] == "A":
                marco += adn[i]
                if len(marco) >= 30:
                    fin = i+1
                    list_marcos.append([marco , inicio , fin])
                    marco=""
                    estado = 0
                    i=inicio+1
                else:
                    marco = ""
                    estado = 0
            else:
                estado = 3
                marco += adn[i]
        i +=1
    return list_marcos

#Funcion que corrige los valores de inicio y fin de la hebra antisentido

def corrige(lista_in,cant):
    for i in range(len(lista_in)):
        lista_in[i][1] = cant - lista_in[i][1] + 1
        lista_in[i][2] = cant - lista_in[i][2] + 1  
    return lista_in

# Funcion que calcula el porcentaje de GC 

def porcentajeGC(sequence):
    if len(sequence)==0:
        return
    else:
        cantGC=0
        for nt in sequence:
            if nt in "GC":
                cantGC +=1
        return (cantGC*100)/len(sequence)

# Función que traduce la secuencia de DNA a proteina

def translate(seq):
    """
    tabla de codones
    """
    tabla = {
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
    for i in range(0, len(seq), 3):
        if len(seq)%3 == 0:
            codon = seq[i : i+3]
            proteina += tabla[codon]
    return proteina

# Función que calcula el % aa distintas clases

def clasesaa(prot):
    polares_neutros=0
    no_polares=0
    negativos=0
    positivos=0
    aromaticos=0
    for aa in prot:
        if aa in "STQNC":
            polares_neutros +=1
        elif aa in "AVLIMPG":
            no_polares += 1
        elif aa in "DE":
            negativos += 1
        elif aa in "KRH":
            positivos += 1
        elif aa in "FYW":
            aromaticos += 1
    clases=[polares_neutros*100/(len(prot)-1),no_polares*100/(len(prot)-1),negativos*100/(len(prot)-1),positivos*100/(len(prot)-1),aromaticos*100/(len(prot)-1)]
    return clases

# Funcion que busca fragmentos de más de 17 aa no polares en la proteina

def buscar_transmembrana (prot):
    list_fragmentos = []
    fragmento = ""
    estado = 0
    index = 1
    for aa in prot:
        if estado == 0:
            if aa in "AVLIMPG":
                estado = 1 
                fragmento += aa
                inicio = index
        elif estado == 1:
            if aa in "AVLIMPG":
                fragmento += aa
            else:
                estado = 0
                if index - inicio >= 17:
                    fin = index-1
                    list_fragmentos.append([fragmento , inicio , fin])
                fragmento = ""
        if index == len(prot) and estado==1 and index - inicio >= 17:
            fin = index
            list_fragmentos.append([fragmento , inicio , fin])
        index += 1
    return list_fragmentos


#Programa-----------------------------------------------------------------------

secuenciafasta=sys.argv[1]

if len(sys.argv)==3:
    archivosalida=sys.argv[2]
else:
    archivosalida="salida.txt"

#Leo el archivo de entrada

sec = leer_seq(secuenciafasta)

#Obtengo la secuencia de la hebra complementaria de ADN

revsec = []
reverse = {"A": "T", "C": "G", "T": "A", "G": "C"}
for i in range (len(sec)):
    revsec.append(reverse[sec[-i -1]])

# Busco los marcos en la hebra sentido y la antisentido

list_orfs_d = buscar_orf(sec)
list_orfs_i = buscar_orf(revsec)
list_orfs_i = corrige(list_orfs_i,len(revsec))

#Calculo el porcentaje de GC

for i in range (0, len(list_orfs_i)):
    list_orfs_i[i].append(porcentajeGC(list_orfs_i[i][0]))

for i in range (0, len(list_orfs_d)):
    list_orfs_d[i].append(porcentajeGC(list_orfs_d[i][0]))

# Traduzco las secuencias de DNA a proteina

for i in range (0, len(list_orfs_i)):
    list_orfs_i[i].append(translate(list_orfs_i[i][0]))

for i in range (0, len(list_orfs_d)):
    list_orfs_d[i].append(translate(list_orfs_d[i][0]))

#Calculo %aa clases

clases= ["polares_neutros","no_polares","negativos","positivos","aromaticos"]

for i in range (0, len(list_orfs_i)):
    list_orfs_i[i].append(clasesaa(list_orfs_i[i][4]))
    plt.pie(x=list_orfs_i[i][5], labels=clases)
    plt.title("Clases de aminoacidos - Hebra antisentido - Prot " + str(i+1))
    plt.savefig("Clases de aminoacidos - Hebra sentido - Prot " + str(i+1))
    plt.close()

for i in range (0, len(list_orfs_d)):
    list_orfs_d[i].append(clasesaa(list_orfs_d[i][4]))
    plt.pie(x=list_orfs_d[i][5], labels=clases)
    plt.title("Clases de aminoacidos - Hebra sentido - Prot " + str(i+1))
    plt.savefig("Clases de aminoacidos - Hebra sentido - Prot" + str(i+1))
    plt.close()

# busca fragmentos de más de 17 aa no polares en l prot

for i in range (0, len(list_orfs_i)):
    list_orfs_i[i].append(buscar_transmembrana(list_orfs_i[i][4]))

for i in range (0, len(list_orfs_d)):
    list_orfs_d[i].append(buscar_transmembrana(list_orfs_d[i][4]))


#Salida-------------------------------------------------------------------------

#Para los marcos en la hebra sentido

fsal1 = open("tabla_hebra_sentido.csv",'w')
fseq1 = open("sec_hebra_sentido.fasta", 'w')
fsal1.write("N°ORF,Inicio,Fin,Longitud ADN,Longitud de proteina,% GC,% aa polares neutros,% aa no polares,% aa negativos,% aa positivos,% aa aromáticos\n")
cont_filas=1
for orfs in list_orfs_d:
    longADN = orfs[2]-orfs[1]+1
    longProt = (longADN-1) // 3
    fila = "{:4d},{:4d},{:4d},{:4d},{:4d},{:6.2f}".format(cont_filas,orfs[1],orfs[2],longADN,longProt,orfs[3])
    for p in orfs[5]:
        fila += ",{:6.2f}".format(p)
    fila += '\n'
    fsal1.write(fila)
    fseq1.write(">ORF{:d} DNA sec\n".format(cont_filas))
    fseq1.write(orfs[0]+"\n\n")
    fseq1.write(">ORF{:d} aa sec\n".format(cont_filas))
    fseq1.write(orfs[4]+"\n\n")
    cont_filas +=1

fsal1.close()
fseq1.close()

#Para los marcos en la hebra antisentido

fsal2 = open("tabla_hebra_antisentido.csv",'w')
fseq2 = open("sec_hebra_antisentido.fasta", 'w')
fsal2.write("N°ORF,Inicio,Fin,Longitud ADN,Longitud de proteina,% GC,% aa polares neutros,% aa no polares,% aa negativos,% aa positivos,% aa aromáticos\n")
cont_filas=1
for orfs in list_orfs_i:
    longADN = orfs[1]-orfs[2]+1
    longProt = (longADN-1) // 3
    fila = "{:4d},{:4d},{:4d},{:4d},{:4d},{:6.2f}".format(cont_filas,orfs[1],orfs[2],longADN,longProt,orfs[3])
    for p in orfs[5]:
        fila += ",{:6.2f}".format(p)
    fila += '\n'
    fsal2.write(fila)
    fseq2.write(">ORF{:d} DNA sec\n".format(cont_filas))
    fseq2.write(orfs[0]+"\n\n")
    fseq2.write(">ORF{:d} aa sec\n".format(cont_filas))
    fseq2.write(orfs[4]+"\n\n")
    cont_filas +=1

fsal2.close()
fseq2.close()