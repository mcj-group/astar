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

# in the build folder, run astar like so
./astar <input file> <startNode> <targetNode> <type> <numThreads> <numQueues> <numBuckets> <delta> <printFull>

# run astar and astar_fine_grain
./astar croatia.bin 0 320970 Serial 1 4 64 0 0
./astar_fine_grain croatia.bin 0 320970 Serial 1 4 64 0 0
```

# Verification
`astar` and `astar_fine_grain` both output a `path.txt` that prints the path from `startNode` to `targetNode`.

To verify the output of a `{map, startNode, targetNode}` configuration:
  1. run the configuration with `Serial` type, this outputs a correct path.txt with the final distances of the path
  2. change the outputted `path.txt` to `mapName.txt`
  3. run the configuration with the desired type
  4. run the command `diff path.txt mapName.txt` to check if the path and distance differs
