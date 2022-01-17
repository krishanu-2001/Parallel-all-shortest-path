# COMMUNICATION AVOIDING ALL-PAIRS SHORTEST PATHS ALGORITHM FOR SPARSE GRAPHS

## TEAM DETAILS:
+ Krishanu Saini  -  190001029
+ Kuldeep Singh  - 190001030
+ Rahul Kumar  -  190001049
+ Keren Tudu   -   190001023

## INSTRUCTIONS TO RUN THE CODE:
1. Download all the files in your operating system and run the files in any text editor
2. First run the classicFW.c code 
3. Run blockedFW.c code by writing "mpiexec -n (number of processes)   ./blockedFW"
4. Run cart.c code to test basic cartesian functions
5. Run graph_gen.cpp code

6. Run 2dsparseFW.h code
7. Run 2dsparseFW.c code by writing "mpiexec -n (number of processes)   ./2dsparseFW"

8.  Run improvement.cpp code to check if improvement of algorithm can be applied
8. Run test.py code to test the efficiency of our algorithm

## REFRENCES:
 1. Alfred V. Aho, John E. Hopcroft, and Jeffrey D. Ullman. 1974. The Design and
Analysis of Computer Algorithms. Addison-Wesley
2. Ed Anderson, Zhaojun Bai, Christian H. Bischof, L. Susan Blackford, James
Demmel, Jack J. Dongarra, Jeremy Du Croz, Anne Greenbaum, Sven Hammarling,A. McKenney, and Danny C. Sorensen. 1999. LAPACK Users’ Guide, Third Edition.
SIAM

## File structure
    
      1. classicFW.c  gcc  
      2. blockedFW.c  MPI  
      3. graph_gen.cpp  g++  
      4. 2dsparseFW.c   MPI  
      5. 2dsparseFW.h   gcc  
      6. improvement.cpp mpic++
        
## Test code     
run in code/ directory  

      1. $ ./make.sh  
      2. $ cd lib  
      3. $ python test.py
