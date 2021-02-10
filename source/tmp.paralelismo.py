#!/bin/python 

# %% 
"""
means_rmvs= []
my_range  =list(set(range(251,500))-set(disconnected)-set(energy_cores))#list(set(range(len(G2.nodes)))-set(disconnected))
"Da un rango de nodos que hace?

1000 nodos

4 notebook con rangos de 250 nodos cada uno
cada notebook se ejecuto en una máquina 

notebook_rango(250)_i (internamente paralelo) -> máquina_i(8 nucleos 64 ram)

1- paralelismo de nodos en 4 maquinas
    2- paralelismo de nucleos y threads por cosa 

with concurrent.futures.ProcessPoolExecutor() as executor: # +400 USD
     for r in executor.map(calc_centr_rmv, my_range):      # +400 USD
         means_rmvs.append(r)                              # +400 USD 

overall_centrality_node_i =sum(abs(delta_C_eigenv_i_subsystem_1), abs(delta_C_harmonic_i_subsystem_1),...,... abs(delta_C_eigenv_i_subsystem_2),abs(delta_C_harmonic_i_subsystem_2),...))
"""

# %% --- Importando modulos base

def squaress(): # Función basica de demostración
    for i in range(50000):
        i ** i
    print(i)

import time

num_processes = 1
print('Empezando de nuevo','\n','Procesos paralelos:', num_processes)

from multiprocessing import Process
processes = []

for i in range(num_processes): # Creamos un proceso por CPU
    p = Process( target= squaress ) # Crea un proceso con la función
    processes.append( p ) # Añade el proceso a la lista de procesos

start = time.time()

for p in processes:
    p.start() 

for p in processes:
    p.join() # Espera que un proceso termine 

end = time.time()

print('ended! Elapsed time:', (end - start)*1000 )

num_processes = 4
print('Empezando de nuevo','\n','Procesos paralelos:', num_processes)

from multiprocessing import Process
processes = []

for i in range(num_processes): # Creamos un proceso por CPU
    p = Process( target= squaress ) # Crea un proceso con la función
    processes.append( p ) # Añade el proceso a la lista de procesos

start = time.time()

for p in processes:
    p.start() 

for p in processes:
    p.join() # Espera que un proceso termine 

end = time.time()

print('ended! Elapsed time:', (end - start)*1000 )

# %%
