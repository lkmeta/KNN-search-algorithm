#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10

cd ../

module load gcc openblas openmpi

make v1
make v2

mkdir results

echo "v1 results" > ./results/v1.txt
	
echo "" >> ./results/v1.txt
for file in ./matrices/*; do
	echo "" >> ./results/v1.txt
	echo "Running v1 for matrix: "$file >> ./results/v1.txt
	for i in 2 4 8 16 20; do
		for k in 10 20 50 100; do
			echo "-------------------" >> ./results/v1.txt
			srun -n $i ./v1 $file $k >> ./results/v1.txt
		done
	done
	echo "" >> ./results/v1.txt
done
echo "" >> ./results/v1.txt


echo "v2 results" > ./results/v2.txt
	
echo "" >> ./results/v2.txt
for file in ./matrices/*; do
	echo "" >> ./results/v2.txt
	echo "Running v2 for matrix: "$file >> ./results/v2.txt
	for i in 2 4 8 16 20; do
		for k in 10 20 50 100; do
			echo "-------------------" >> ./results/v2.txt
			srun -n $i ./v2 $file $k >> ./results/v2.txt
		done
	done
	echo "" >> ./results/v2.txt
done
echo "" >> ./results/v2.txt
