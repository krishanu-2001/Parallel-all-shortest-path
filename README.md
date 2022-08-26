# COMMUNICATION AVOIDING ALL-PAIRS SHORTEST PATHS ALGORITHM FOR SPARSE GRAPHS

## Introduction
The All Source Shortest Path problem, initially solved by Floyd-Warshall algorithm is a variation of matrix multiplication problem.  
This problem emerges in the field of computer network topology, architechture and biological networks. Utilizing novel algorithms like strassen matrix multiplication, 2.5D approaches one can improve the efficiency and speedup.  
Bandwidth and latency cost is the major scope of optimization for this parallel program. In this implementation we provide an architecture for sparse networks using state of the art research techniques.

## INSTRUCTIONS TO RUN THE CODE:
1. Download all the files in your operating system and run the files in any text editor
2. First run the classicFW.c code 
3. Run blockedFW.c code by writing "mpiexec -n (number of processes)   ./blockedFW"
4. Run cart.c code to test basic cartesian functions
5. Run graph_gen.cpp code

6. Run 2dsparseFW.h code
7. Run 2dsparseFW.c code by writing "mpiexec -n (number of processes)   ./2dsparseFW"

8. Run improvement.cpp code to check if improvement of algorithm can be applied
9. Run test.py code to test the efficiency of our algorithm

## Contribute
First of all, thanks for your contribution! Every small bit of it counts! You can:

1. Create a new issue for bugs, feature requests, and enhancements.
2. Fork the repo, make changes, and submit a pull request, describing the changes made.
3. Discuss new approaches and use cases with the community.
4. Do star the repo if you like it!

## REFRENCES:
1. Lin Zhu, Qiang-Sheng Hua, and Hai Jin. 2021. Communication Avoiding All-Pairs Shortest Paths Algorithm for Sparse Graphs. In 50th International Conference on Parallel Processing (ICPP '21), August 9–12, 2021, Lemont, IL, USA. ACM, New York, NY, USA, 10 pages. https://doi.org/10.1145/3472456.3472524
2. Alfred V. Aho, John E. Hopcroft, and Jeffrey D. Ullman. 1974. The Design and
Analysis of Computer Algorithms. Addison-Wesley
3. Ed Anderson, Zhaojun Bai, Christian H. Bischof, L. Susan Blackford, James
Demmel, Jack J. Dongarra, Jeremy Du Croz, Anne Greenbaum, Sven Hammarling,A. McKenney, and Danny C. Sorensen. 1999. LAPACK Users’ Guide, Third Edition.
SIAM
