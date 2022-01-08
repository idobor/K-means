# K-means
C + Python Implementation of the normalized spectral clustering algorithm (Based on k-means)

Algorithm:
Given a set of n points X:
![image](https://user-images.githubusercontent.com/61387800/148657566-6b8355af-7812-49f5-ae55-f3de3dbee55d.png)


# How to use?
Python spkmeans.py,and C spkmeans.c
Reads user CMD arguments:

(a) k (int, < N): Number of required clusters. If equal 0, the program will find the optimal K.
(b) goal (Enum): 
Can get the following values:
i. spk: Perform full spectral kmeans.
ii. wam: Calculate and output the Weighted Adjacency Matrix. 
iii. ddg: Calculate and output the Diagonal Degree Matrix
iv. lnorm: Calculate and output the Normalized Graph Laplacian.
v. jacobi: Calculate and output the eigenvalues and eigenvectors.
(c) file name (.txt or .csv): The path to the file that will contain N observations, the file
extension is .txt or .csv.

C exmaple:
![image](https://user-images.githubusercontent.com/61387800/148657719-c031e920-3326-4f53-8a30-dd0388dbaf86.png)


# How to compile?
First build module extension:
**python setup.py build_ext --inplace**

Then compile spkmeans.c with GCC
