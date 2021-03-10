#!/bin/bash
### Este script funciona para cualquier numero de nodos; Ray se encarga de manejarlos
### puedes definir este parametro en la ejecución de la tarea como 'sbatch submit-ray... --nodes=5'
#SBATCH --nodes=5

### Asigna todos los recursos de un nodo al proceso Ray, el cual los maneja internamente
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-task=0
#SBATCH --cpus-per-task=48

# Carga modulos aqui
# TODO: revisar como funciona esto en NLHPC
# module load miniconda/ 

################# DON NOT CHANGE THINGS HERE UNLESS YOU KNOW WHAT YOU ARE DOING ###############
# This script is a modification to the implementation suggest by gregSchwartz18 here:
# https://github.com/ray-project/ray/issues/826#issuecomment-522116599

echo "INFO Empieza el setup: $(date)"

redis_password=$(uuidgen)
export redis_password

nodes=$(scontrol show hostnames $SLURM_JOB_NODELIST) # Getting the node names
nodes_array=( $nodes )

node_1=${nodes_array[0]} 
ip=$(srun --nodes=1 --ntasks=1 -w $node_1 hostname) # making redis-address
port=6379
ip_head=$ip:$port
export ip_head
echo "IP Head: $ip_head"

echo "INFO Empezando HEAD en $node_1"
srun --nodes=1 --ntasks=1 -w $node_1 ./source/ray/start-head.sh $ip $redis_password &
sleep 30

worker_num=$(($SLURM_JOB_NUM_NODES - 1)) #number of nodes other than the head node
for ((  i=1; i<=$worker_num; i++ ))
do
  node_i=${nodes_array[$i]}
  echo "INFO Empezando WORKER $i en $node_i"
  srun --nodes=1 --ntasks=1 -w $node_i ./source/ray/start-worker.sh $ip_head $redis_password &
  sleep 5
done

#### SECCIÓN QUE SI EJECUTA CODIGO PYTHON Y DEMAS

echo "INFO Empieza el proceso paralelo: $(date)"

###  Helper que cuenta cuantos cores hay disponibles
n_cores=$(sinfo -o "%C" -h | awk -F '/' '{print $2}')

echo "INFO Corriendo en $n_cores nodos"

echo "INFO Status del cluster Ray"
ray status --address $ip_head --redis_password $redis_password

ray timeline --address $ip_head 

time python -u ./source/delta_centrality.py './data/stimulated_2021.json' $n_cores

### INFORMACIÓN FINAL QUE SALE AL OUTPUT DEL SISTEMA
# Esta parte está fallando, así que mejor imprimir algo como un grep o algo así porsiaca

echo "INFO Termina el proceso paralelo: $(date)"

sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize

exit

