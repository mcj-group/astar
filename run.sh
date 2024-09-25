#!/bin/bash

exec_dir=./build
input_dir=../../inputs_astar
startnode=0
endnode=320970

# MQBucket
echo "Results" > results_astar_1_germany_MQBucket.log
for j in {1..10}
do
	${exec_dir}/astar ${input_dir}/germany.bin $startnode $endnode MQBucket 1 4 4 4 11 64 0 >> results_astar_1_germany_MQBucket.log
done
mv path.txt MQBucket_1_germany_path.txt

echo "Results" > results_astar_48_germany_MQBucket.log
for j in {1..10}
do
        ${exec_dir}/astar ${input_dir}/germany.bin $startnode $endnode MQBucket 48 192 4 4 11 64 0 >> results_astar_48_germany_MQBucket.log
done
mv path.txt MQBucket_48_germany_path.txt

# MQ
echo "Results" > results_astar_1_germany_MQ.log
for j in {1..10}
do
        ${exec_dir}/astar ${input_dir}/germany.bin $startnode $endnode MQ 1 4 2 2 >> results_astar_1_germany_MQ.log
done
mv path.txt MQ_1_germany_path.txt

echo "Results" > results_astar_48_germany_MQ.log
for j in {1..10}
do
        ${exec_dir}/astar ${input_dir}/germany.bin $startnode $endnode MQ 48 192 2 2 >> results_astar_48_germany_MQ.log
done
mv path.txt MQ_48_germany_path.txt


# MQ Plain
echo "Results" > results_astar_1_germany_MQPlain.log
for j in {1..10}
do
        ${exec_dir}/astar ${input_dir}/germany.bin $startnode $endnode MQ 1 4 1 1 >> results_astar_1_germany_MQPlain.log
done
mv path.txt MQPlain_1_germany_path.txt

echo "Results" > results_astar_48_germany_MQPlain.log
for j in {1..10}
do
        ${exec_dir}/astar ${input_dir}/germany.bin $startnode $endnode MQ 48 192 2 2 >> results_astar_48_germany_MQPlain.log
done
mv path.txt MQPlain_48_germany_path.txt

