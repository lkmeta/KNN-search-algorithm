# KNN search-algorithm
# **Parallel and Distributed Computer Systems**  
## **Exercise 2**
## Tsoumplekas Georgios | gktsoump@ece.auth.gr  <br />
## Meta Louis-Kosmas    | louismeta@ece.auth.gr <br />



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
   + ``` mpirun -np <Number of tasks> ./v1 <path to matrix> <number of nearest neighbors> ```
   
**Note:** When running the functions with one of the given matrices it is important to not have changed the name of the matrix previously, since the reader can only
recognize the original names.

### **Process for tested Matrices** 

**Note:** We used the following [Scripts](https://github.com/lkmeta/KNN-search-algorithm/tree/main/bash%20scripts) in AUTh High Performance Computing (HPC) infrastructure to test our code.  
We used different time durations for these scripts to compute the [results](https://github.com/lkmeta/KNN-search-algorithm/tree/main/results) for the tested matrices.  
To build and run these scripts you should contain a subfolder named "matrices" which will have the tested matrices.  

**Script ```parallel.sh```:**
   + calculates KNN using v1 and v2 with the matrices for the following combinations of tasks and k value
   + number of tasks: 2 4 8 16 20 
   + number of nearest neighbors: 5 10 20 50 100
   + creates the directory results and outputs for v1 and v2 files

**Script ```sequential.sh```:**
   + calculates KNN using v0 with the matrices for the following combinations of k value
   + number of nearest neighbors: 10 20 50 100
   + creates the directory results and outputs for v0 file





