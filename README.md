# astar
Repo for cpp astar algorithm, based on the implementation from swarm & chronos.

# Build & Run
Need to include the cps repo for MultiQueue and BucketMultiQueue.

`CMakeLists.txt` includes the cps repo by `current_directory/../cps/include`, modify the path if needed.

```
# after cloning the astar repo and cps repo
./build.sh

# in the build folder, run astar like so,
# arguments must follow this order
./astar <inFile> [startNode endNode qType threadNum queueNum batchPop batchPush delta bucketNum printFull]

# to run the experiments
./run.sh

# run astar with default arguments on croatia
./astar croatia.bin 

# run astar serial
./astar croatia.bin 0 320970 Serial

# run astar with MultiQueue and 16 threads
./astar croatia.bin 0 320970 MQ 16

# run astar with BucketMQ, 16 threads, 64 queues, batchsizes of 4, delta of 11 and 64 buckets per bucket queue
./astar croatia.bin 0 320970 MQBucket 16 64 4 4 11 64 
```

# Verification
`astar` outputs a `path.txt` that prints the path from `startNode` to `targetNode`.

To verify the output of a `{map, startNode, targetNode}` configuration:
```
# 1. run the configuration with `Serial` type
#    this outputs a correct path.txt with the final distances of the path
./astar <map> <startNode> <targetNode> Serial

# 2. save the outputted path file with a different name
mv path.txt mapName.txt

# 3. run the configuration with the desired type and threads, etc
./astar <map> <startNode> <targetNode> <desiredType> <threadNum>

# 4. compare the path and distances
diff path.txt mapName.txt

```
  
  
  
  
