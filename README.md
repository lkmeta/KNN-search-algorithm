# KNN search-algorithm
# **Parallel and Distributed Computer Systems**  
## **Exercise 2**
## Tsoumplekas Georgios | gktsoump@ece.auth.gr  <br />
## Meta Louis-Kosmas    | louismeta@ece.auth.gr <br />

---
### Contents
   1. [How to build](#1)
   2. [Results from HPC](#2)
   3. [Graphs](#3)
   
<a name="1"></a>
### **How to build**  

**File ```Makefile```:**  
   + Use command line ```make all``` to build v0, v1 and v2 immediately.
   + Use command line ```make testv0``` to test v0 sequential with random array (n = 10000, d = 5, k = 20). 
   + Use command line ```make testv1``` to test v1 Asynchronous with random array (n = 10000, d = 5, k = 20, numOfTasks = 4).
   + Use command line ```make testv2``` to test v2 version with random array (n = 10000, d = 5, k = 20, numOfTasks = 4).
   
**Examples:**  
The functions can either accept 3 or 4 arguments.  
If you want to run a function with a matrix of random values run as shown below:
   + ```./v0 <Number of data points> <Number of dimensions> <Number of nearest neighbors> ```  
   + ```mpirun -np <Number of tasks> ./v1 <Number of data points> <Number of dimensions> <Number of nearest neighbors> ```  
   + ```mpirun -np <Number of tasks> ./v2 <Number of data points> <Number of dimensions> <Number of nearest neighbors> ```  

If you want to run a function with one of the given matrices run as shown below:
   + ```./v0 <Matrix Path> <Number of nearest neighbors> ```   
   + ``` mpirun -np <Number of tasks> ./v1 <Matrix Path> <Number of nearest neighbors> ```
   + ``` mpirun -np <Number of tasks> ./v2 <Matrix Path> <Number of nearest neighbors> ```
   
**Note:** When running the functions with one of the given matrices it is important to not have changed the name of the matrix previously, since the reader can only
recognize the original names.

---
<a name="2"></a>
### **Results from HPC** 

**Note:** We used the following [Scripts](https://github.com/lkmeta/KNN-search-algorithm/tree/main/bash%20scripts) in AUTh High Performance Computing (HPC) infrastructure to test our code.  
We used different time durations for these scripts to compute the [results](https://github.com/lkmeta/KNN-search-algorithm/tree/main/results) for the tested matrices.  
To build and run these scripts you should contain a subfolder named "matrices" which will have the tested matrices.  

**Script [```parallel.sh```](https://github.com/lkmeta/KNN-search-algorithm/blob/main/bash%20scripts/parallel.sh):**
   + calculates KNN using v1 and v2 with the matrices for the following combinations of tasks and k value
   + number of tasks: 2 4 8 16 20 
   + number of nearest neighbors: 10 20 50 100
   + creates the directory results and outputs for v1 and v2 files

**Script [```sequential.sh```](https://github.com/lkmeta/KNN-search-algorithm/blob/main/bash%20scripts/sequential.sh):**
   + calculates KNN using v0 with the matrices for the following combinations of k value
   + number of nearest neighbors: 10 20 50 100
   + creates the directory results and outputs for v0 file


**_Code Validation_**  
The sequential code for the KNN search algorithm we used has been validated succesfully from automated tester on e-learning.  
Furthermore, we can integrate the tester from e-learning in our code by using the functions in [test](https://github.com/lkmeta/KNN-search-algorithm/tree/main/test) subfolder.
 
---
<a name="3"></a>
### **Graphs** 
We used all the data we created from AUTh High Performance Computing (HPC) infrastructure as we mentioned with this [Matlab Script]() to build the Graphs.
For every tested matrix we have 8 graphs:
   + v0 execution time depending on k value
   + v1 execution time depending on number of processes value for each k
   + v1 acceleration depending on number of processes value for each k
   + v1 for 2 and 20 processes depending on k value
   + v2 execution time depending on number of processes value for each k
   + v2 acceleration depending on number of processes value for each k
   + v2 for 2 and 20 processes depending on k value
   + v1 and v2 execution time depending on number of processes

You can find [here](https://github.com/lkmeta/KNN-search-algorithm/tree/main/diagrams) all the graphs as we described before.

Note: We faced some problems with 2 matrices when we tried to read and test them with the Scripts.
Specially, for "features.csv" matrix we couldn't get the needed data for the graphs. 
Also, for "BBC.txt" we couldn't run some of the combinations we wanted to test.

