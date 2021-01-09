
# Module: modulo.py
# Version 0.0.1


def metabolites2reactions(modelo, metabolitos): # Toma metabolitos
  """Crea un agregado de las reacciones segun un listado de metabolitos. 

  Parameters
  ----------
  modelo: model
      El modelo inicial del cual se toman las reacciones
  metabolitos : list
      IDs de un metabolito. La función acepta multiples para crear un joint. 

  Returns
  -------
  df
      crea un dataframe pandas agregada de ids, nombre de reacciones, y reacciones en si
  """
  frozen = modelo.metabolites.get_by_id(metabolitos[0]).reactions # Crea el inicio
  for met in metabolitos: # Añade metabolistos al listado de reacciones
    frozen = frozen.union(modelo.metabolites.get_by_id(met).reactions)
  
  ids   = [list(frozen)[i].id       for i in range(0,len(frozen))] # List Comprehesion
  names = [list(frozen)[i].name     for i in range(0,len(frozen))]
  rxns  = [list(frozen)[i].reaction for i in range(0,len(frozen))]

  list_of_tuples = list(zip(ids, names, rxns))     # DataFrame en Python; crea triplets
  columns        = ["ids","names","rxns"]   # Titulos de columnas de la tripleta
  df             =  pd.DataFrame(list_of_tuples, columns = columns) # Titulos asignados a la lista de tuplas
  
  return df
# 