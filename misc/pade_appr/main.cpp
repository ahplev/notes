#include <cstdlib>
#include <ginac/ex.h>
#include <ginac/ginac.h>
#include <iostream>
#include <string>
#include <pthread.h>

#define NUM_THREADS 12
#define NUM_X 786432
#define START -3
#define LEN 11.

using namespace ::std;
using namespace ::GiNaC;

struct fps {
  FUNCP_1P double_f;
  FUNCP_1P double_T;
  FUNCP_1P double_pade;
};

struct parg {
  fps fp;
  double* value;
  int thread_id;
  double st;
  double step;
};

fps compile(int m, int n);
void* routine(void *arg);

int main(int args, char* argv[]) {
  if (args != 3) {
    printf("wrong arguments");
    return 1;
  }
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);

  fps fp = compile(m, n);

  char *res = (char*)malloc((4*NUM_X+1) * 27 * sizeof(char));
  double *value = (double*)malloc(3 * NUM_X * sizeof(double));

  pthread_t threads[NUM_THREADS];
  parg pargs[NUM_THREADS];

  for (int i = 0; i < NUM_THREADS; i++) {
    pargs[i].fp = fp;
    pargs[i].value = value;
    pargs[i].thread_id = i;
    pargs[i].st = START;
    pargs[i].step = LEN/NUM_X;

    if (pthread_create(threads + i, NULL, routine, (void *)(pargs + i)) != 0) {
      perror("pthread_create error");
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    if (pthread_join(threads[i], NULL) != 0) {
      perror("pthread_join error");
      exit(EXIT_FAILURE);
    }
  }

  long long len = 0;
  for (int i = 0; i < NUM_X; i++) {
    len += sprintf(res + len, "%f\t", START + i * LEN/NUM_X);
    for (int j = 0; j < 3; j++) {
      len += sprintf(res + len, "%f\t", value[i * 3 + j]);
    }
    len += sprintf(res + len, "\n");
  }
  free(value);

  FILE *fptr = fopen("val.dat", "w");
  fprintf(fptr, "%s", res);
  fclose(fptr);
  free(res);

 return 0;
}

fps compile(int m, int n) {
  symbol x("x");
  ex f = sin(x);

  ex T = series_to_poly(f.series(x ==  0, m + n + 1));

  ex P = 0, Q = 1;
  symbol *a = new symbol[m + 1];
  symbol *b = new symbol[n + 1];

  lst eqs, vars;

  for (int j = 0; j <= m; j++) {
    a[j] = symbol("a" + to_string(j));
    P += a[j] * pow(x, j);
    vars.append(a[j]);
  }

  for (int k = 1; k <= n; k++) {
    b[k] = symbol("b"+to_string(k));
    Q += b[k]*pow(x, k);
    vars.append(b[k]);
  }
  ex product = expand(Q * T);
  ex equation = P - product;

  for (int i = 0; i <= m + n; i++) {
    eqs.append(coeff(equation, x, i) == 0);
  }
  auto solutions = lsolve(eqs, vars);
  ex pade = (P/Q).subs(solutions);
  cout << solutions << '\n';

  FUNCP_1P double_f, double_T, double_pade;
  compile_ex(f, x, double_f);
  compile_ex(T, x, double_T);
  compile_ex(pade, x, double_pade);

  return fps{double_f, double_T, double_pade};
}

void* routine(void* arg) {
  parg thisarg = *(parg*)arg;
  int t = (NUM_X + NUM_THREADS - 1) / NUM_THREADS;
  for (int i = 0; i < t; i++) {
    int column = thisarg.thread_id * t + i; 
    if (column >= NUM_X) continue;
    double x = thisarg.st + column * thisarg.step;
    thisarg.value[3 * column + 0] = thisarg.fp.double_f(x);
    thisarg.value[3 * column + 1] = thisarg.fp.double_T(x);
    thisarg.value[3 * column + 2] = thisarg.fp.double_pade(x);
  }
  return NULL;
}
