#ifndef MYCLOCK_H
#define MYCLOCK_H

#include <mpi.h>
#include <stdio.h>

class MyClock {
public:
  /** default constructor */
  MyClock(int _nt) {
    tTot = MPI_Wtime();
    nt = _nt;
    t        = new double[nt];
    ttotal = new double[nt];
    ns       = new int[nt];
    active   = new bool[nt];
    for (int i=0; i<nt; i++) {
      t[i] = ttotal[i] = 0;
      ns[i] = 0;
      active[i] = false;
    }
  }
  
  ~MyClock() {
    delete t;
    delete ttotal;
    delete ns;
    delete active;
  }
  
  void start(int i) {
    active[i] = true;
    t[i] = MPI_Wtime();
  }

  void stop(int i) {
    if (active[i]) {
      double t1 = MPI_Wtime();
      t[i] = t1 - t[i];
      ttotal[i] += t[i];
      ns[i]++;
      active[i] = false;
    }
  }
  
  void print(int i) {
    printf("%7.3f", ttotal[i]); 
  }

  double get(int i) { return ttotal[i];}
  

  void printTotal(){
    double t1 = MPI_Wtime();
    printf(" Total time: %.3f \n", t1 - tTot);
  }

  void print() {
    MPI_Barrier(MPI_COMM_WORLD);
    int rnk;
    MPI_Comm_rank (MPI_COMM_WORLD, &rnk);  
    if (!rnk) {
    double tot = 0;
    printf("##############################\n"
           "#  Profiling: \n"
           "#        last      avg     tot\n");
    for(int i=0; i<nt; i++) {print(i); tot+=ttotal[i];}
    printf("\n                          %7.3f\n",tot);
    printTotal();
    printf("\n##############################\n");
    }
  }

private:
  int nt;
  double tTot;
  double *t;
  double *ttotal;
  int    *ns;
  bool   *active;

};

#endif
