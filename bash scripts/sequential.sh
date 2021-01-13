#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

cd ../

module load gcc openblas

make v0

mkdir results

echo "v0 results" > ./results/v0.txt
	
for file in ./matrices/*; do
	echo "" >> ./results/v0.txt
	echo "Running v0 for matrix: "$file >> ./results/v0.txt
    for k in 5 10 20 50 100; do
		./v0 $file $k >> ./results/v0.txt
	done
	echo "" >> ./results/v0.txt
done
