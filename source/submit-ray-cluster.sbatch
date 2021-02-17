#!/bin/bash
### This script works for any number of nodes, Ray will find and manage all resources
#SBATCH --nodes=5

### Give all resources to a single Ray task, ray can manage the resources internally
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-task=0
#SBATCH --cpus-per-task=44

# Load modules or your own conda environment here
# module load pytorch/v1.4.0-gpu 

################# DON NOT CHANGE THINGS HERE UNLESS YOU KNOW WHAT YOU ARE DOING ###############
# This script is a modification to the implementation suggest by gregSchwartz18 here:
# https://github.com/ray-project/ray/issues/826#issuecomment-522116599
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

echo "STARTING HEAD at $node_1"
srun --nodes=1 --ntasks=1 -w $node_1 ./source/ray/start-head.sh $ip $redis_password &
sleep 30

worker_num=$(($SLURM_JOB_NUM_NODES - 1)) #number of nodes other than the head node
for ((  i=1; i<=$worker_num; i++ ))
do
  node_i=${nodes_array[$i]}
  echo "STARTING WORKER $i at $node_i"
  srun --nodes=1 --ntasks=1 -w $node_i ./source/ray/start-worker.sh $ip_head $redis_password &
  sleep 5
done

#### SECCIÓN QUE SI EJECUTA CODIGO PYTHON Y DEMAS
###  Helper que cuenta cuantos cores hay disponibles
n_cores=$(sinfo -o "%C" -h | awk -F '/' '{print $2}')

echo "Status del cluster Ray"
ray status --address $ip_head --redis_password $redis_password

time python -u ./source/delta_centrality.py './data/stimulated_2021.json' $n_cores

exit

