#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "util.h"


void timestamp(void)
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf("%s\n", time_buffer);

  return;
#undef TIME_SIZE
}

double potential(double a, double b, double c, double x, double y, double z)
{
  return 2.0 * (pow(x / a / a, 2) + pow(y / b / b, 2) + pow(z / c / c, 2)) + 1.0 / a / a + 1.0 / b / b + 1.0 / c / c;
}

// generator pseudoslučajnih brojeva po uniformnoj raspodeli - svaka nit ima svoj seed, jer kada bi bio shared, uniformnost ne bi bila garantovana
double r8_uniform_01(int *seed)
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * (*seed - k * 127773) - k * 2836;

  if (*seed < 0)
  {
    *seed = *seed + 2147483647;
  }
  r = (double)(*seed) * 4.656612875E-10;

  return r;
}

double feynman_1(const double a, const double b, const double c, const double h, const double stepsz, const int ni, const int nj, const int nk, const int N) 
{
  int seed = 123456789;
  double err = 0.0;
  int n_inside = 0;

#pragma omp parallel default(none) shared(a, b, c, h, stepsz, ni, nj, nk, N) \
                                   firstprivate(seed) \
                                   reduction(+ : err) \
                                   reduction(+ : n_inside)
{
  // seed is private variable, so numbers can generate uniformly
  seed += omp_get_thread_num();

#pragma omp for collapse(3)
  for (int i = 1; i <= ni; i++)
  {
    for (int j = 1; j <= nj; j++)
    {
      for (int k = 1; k <= nk; k++)
      {
        double x = ((double)(ni - i) * (-a) + (double)(i - 1) * a) / (double)(ni - 1);
        double y = ((double)(nj - j) * (-b) + (double)(j - 1) * b) / (double)(nj - 1);
        double z = ((double)(nk - k) * (-c) + (double)(k - 1) * c) / (double)(nk - 1);
        double chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);

        if (1.0 < chk)
        {

// deo oznacen za DEBUG se kompajlira ako se prilikom prevodjenja navede opcija -DDEBUG
#ifdef DEBUG
          printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
                 x, y, z, 1.0, 1.0, 0.0, 0);
#endif

          continue;
        }

        n_inside++;

        double w_exact = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);
        double wt = 0.0;

#ifdef DEBUG
        int steps = 0;
#endif

        for (int trial = 0; trial < N; trial++)
        {
          double x1 = x;
          double x2 = y;
          double x3 = z;
          double w = 1.0;
          chk = 0.0;
          while (chk < 1.0)
          {
            double ut = r8_uniform_01(&seed);
            double us;

            double dx;
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
                dx = -stepsz;
              else
                dx = stepsz;
            }
            else
              dx = 0.0;

            double dy;
            ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
                dy = -stepsz;
              else
                dy = stepsz;
            }
            else
              dy = 0.0;

            double dz;
            ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
                dz = -stepsz;
              else
                dz = stepsz;
            }
            else
              dz = 0.0;

            double vs = potential(a, b, c, x1, x2, x3);
            x1 = x1 + dx;
            x2 = x2 + dy;
            x3 = x3 + dz;

#ifdef DEBUG
            ++steps;
#endif

            double vh = potential(a, b, c, x1, x2, x3);

            double we = (1.0 - h * vs) * w;
            w = w - 0.5 * h * (vh * we + vs * w);

            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
          }
          wt = wt + w;
        }
        wt = wt / (double)(N);

        err += pow(w_exact - wt, 2);

#ifdef DEBUG
        printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
               x, y, z, wt, w_exact, fabs(w_exact - wt), steps / N);
#endif
      }
    }
  }
} // parallel
  return sqrt(err / (double)(n_inside));
  //printf("\n\nRMS absolute error in solution = %e\n", err);
}

double feynman_2(const double a, const double b, const double c, const double h, const double stepsz, const int ni, const int nj, const int nk, const int N) 
{
  static int seed = 123456789;  // set to be static -> to be global (shared) by default
  int n_inside = 0;
  double w_exact[17][12][7];
  double wt[17][12][7];

#pragma omp threadprivate(seed) // then, it is private for every thread
#pragma omp parallel default(none) shared(a, b, c, h, stepsz, ni, nj, nk, N, n_inside, w_exact, wt)
{
  seed += omp_get_thread_num();

#pragma omp single
{
  for (int i = 1; i <= ni; i++)
  {
    for (int j = 1; j <= nj; j++)
    {
      for (int k = 1; k <= nk; k++)
      {
        double x = ((double)(ni - i) * (-a) + (double)(i - 1) * a) / (double)(ni - 1);
        double y = ((double)(nj - j) * (-b) + (double)(j - 1) * b) / (double)(nj - 1);
        double z = ((double)(nk - k) * (-c) + (double)(k - 1) * c) / (double)(nk - 1);
        double chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);
        w_exact[i][j][k] = 0.0;
        wt[i][j][k] = 0.0;

        if (1.0 < chk)
        {
          continue;
        }

        n_inside++;

        w_exact[i][j][k] = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);

        for (int trial = 0; trial < N; trial++)
        {
#pragma omp task shared(wt)
{
          double x1 = x;
          double x2 = y;
          double x3 = z;
          double w = 1.0;
          chk = 0.0;
          while (chk < 1.0)
          {
            double ut = r8_uniform_01(&seed);
            double us;

            double dx;
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
                dx = -stepsz;
              else
                dx = stepsz;
            }
            else
              dx = 0.0;

            double dy;
            ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
                dy = -stepsz;
              else
                dy = stepsz;
            }
            else
              dy = 0.0;

            double dz;
            ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
                dz = -stepsz;
              else
                dz = stepsz;
            }
            else
              dz = 0.0;

            double vs = potential(a, b, c, x1, x2, x3);
            x1 = x1 + dx;
            x2 = x2 + dy;
            x3 = x3 + dz;

            double vh = potential(a, b, c, x1, x2, x3);

            double we = (1.0 - h * vs) * w;
            w = w - 0.5 * h * (vh * we + vs * w);

            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
          }
#pragma omp atomic
          wt[i][j][k] += w;
} // task
        }
      }
    }
  }
} // single
} // parallel
  // na kraju obracunati gresku po svim osama - u matrici wt
  double err = 0.0;
  for (int i = 0; i <= 16; ++i)
  {
      for (int j = 0; j <= 11; ++j)
      {
          for (int k = 0; k <= 6; ++k)
          {
              if (w_exact[i][j][k] == 0.0)
              {
                  continue;
              }
              err += pow(w_exact[i][j][k] - (wt[i][j][k] / (double)(N)), 2);
          }
      }
  }
  return sqrt(err / (double)(n_inside));
}


double (*FUNCS[])(const double, const double, const double, const double, const double, const int, const int, const int, const int) = {feynman_1, feynman_2};

// print na stdout upotrebiti u validaciji paralelnog resenja
int main(int argc, char **argv)
{
  const double a = 3.0;
  const double b = 2.0;
  const double c = 1.0;
  const double h = 0.001;
  // broj tacaka po svakoj dimenziji
  const int ni = 16;
  const int nj = 11;
  const int nk = 6;

  const double stepsz = sqrt(3 * h);;

  if (argc < 3)
  {
    printf("Invalid number of arguments passed.\n");
    return 1;
  }

  const int func = atoi(argv[1]);
  const int N = atoi(argv[2]);

  printf("TEST: func=%d, N=%d, num_threads=%ld\n", func, N, get_num_threads());
  double wtime = omp_get_wtime();
  double err = FUNCS[func](a, b, c, h, stepsz, ni, nj, nk, N);
  wtime = omp_get_wtime() - wtime;
  printf("%d    %lf    %lf\n", N, err, wtime);
  printf("TEST END\n");

  return 0;
}

// commands: 
// make
/*
./feyman 0 1000
./feyman 0 5000
./feyman 0 10000
./feyman 0 20000
 */
