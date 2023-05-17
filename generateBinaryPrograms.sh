
nvcc -std=c++11 -o bin/mainBad2Good src/main_Bad2Good.cu src/mergeClusteringGB.cu -Xcompiler -fopenmp -O3

nvcc -std=c++11 -o bin/mainCluster2DemRes src/main_clusters2demRes.cu src/mergeClusteringGB.cu -Xcompiler -fopenmp -O3

nvcc -std=c++11 -o bin/mainMergeCluster src/main_mergeGoodCluster.cu src/mergeClusteringGB.cu -Xcompiler -fopenmp -O3

nvcc -std=c++11 -o bin/mainRefinOneCluster src/main_refineOneCluster.cu src/mergeClusteringGB.cu -Xcompiler -fopenmp -O3
