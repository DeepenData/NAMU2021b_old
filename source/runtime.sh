#!/bin/bash
#SBATCH --job-name=test-ray-py  # Nombre del Trabajo para el cluster

#SBATCH --nodes=3               # Cantidad de nodos worker
#SBATCH --exclusive             # Uso exclusivo de esos nodos (?)
#SBATCH --tasks-per-node=1      # Una instancia Ray por nodo

#SBATCH --cpus-per-task=2       # CPUs por nodo
#SBATCH --mem-per-cpu=0.5GB     # RAM por CPU

conda activate human-metnet     # El entorno de Conda con Ray

# --- Define cosas del entorno
nodes=$(scontrol show hostnames "$SLURM_JOB_NODELIST")
nodes_array=($nodes)

head_node=${nodes_array[0]} # Nodo principal con el server de Ray
head_node_ip=$(srun --nodes=1 --ntasks=1 -w "$head_node" hostname --ip-address)

# Paso opcional en caso de que el IP incluya ' '
if [[ "$head_node_ip" == *" "* ]]; then
IFS=' ' read -ra ADDR <<<"$head_node_ip"
if [[ ${#ADDR[0]} -gt 16 ]]; then
  head_node_ip=${ADDR[1]}
else
  head_node_ip=${ADDR[0]}
fi
echo "IPV6 address detected. IPV4 es $head_node_ip"
fi

# --- Inicializa el server de Ray
port=6379
ip_head=$head_node_ip:$port
export ip_head
echo "IP Head: $ip_head"

echo "Inicializando servidor Ray en: $head_node"
srun --nodes=1 --ntasks=1 -w "$head_node" \
    ray start --head --node-ip-address="$head_node_ip" --port=$port \
    --num-cpus "${SLURM_CPUS_PER_TASK}" --num-gpus "${SLURM_GPUS_PER_TASK}" --block &

# --- Inicia los nodos del server de Ray
worker_num=$((SLURM_JOB_NUM_NODES - 1))

for ((i = 1; i <= worker_num; i++)); do
    node_i=${nodes_array[$i]}
    echo "Inicializando WORKER $i at $node_i"
    srun --nodes=1 --ntasks=1 -w "$node_i" \
        ray start --address "$ip_head" \
        --num-cpus "${SLURM_CPUS_PER_TASK}" --num-gpus "${SLURM_GPUS_PER_TASK}" --block &
    sleep 5
done

# --- Invocando el Script Python

time python -u source/delta_centrality.py './data/toy_metabolism_AA.json' 6
