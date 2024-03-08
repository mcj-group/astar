# astar
Repo for cpp astar algorithm, based on the implementation from swarm-benchmarks.

# Build & Run
Need to include the CPS repo for MultiQueue and BucketMultiQueue.

`CMakeLists.txt` includes the CPS repo by `current_directory/../CPS/include`, modify the path if needed.

```
# after cloning the astar repo and CPS repo
mkdir astar/build
cd astar/build/
cmake .. -DCMAKE_BUILD_TYPE=Release
make astar -j
make astar_fine_grain -j

# in the build folder, run astar like so,
# arguments must follow this order
./astar <inFile> [startNode endNode qType threadNum bucketNum printFull]

# run astar and astar_fine_grain with default arguments on croatia
./astar croatia.bin 
./astar_fine_grain croatia.bin

# run astar serial
./astar croatia.bin 0 320970 Serial

# run astar with MultiQueue and 16 threads
./astar croatia.bin 0 320970 MQIO 16

# run astar with BucketMQ, 16 threads, 128 buckets, and delta of 10
./astar croatia.bin 0 320970 MQBucket 16 128 10 
```

# Verification
`astar` and `astar_fine_grain` both output a `path.txt` that prints the path from `startNode` to `targetNode`.

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
  
  
  
  
