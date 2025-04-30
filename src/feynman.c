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

// potencijalna energiju u nekom centralno simetricnom polju unutar elipsoida
// energija veca sto je tacka dalje od centra i sto su dimenzije elipsoida manje
double potential(double a, double b, double c, double x, double y, double z)
{
  return 2.0 * (pow(x / a / a, 2) + pow(y / b / b, 2) + pow(z / c / c, 2)) + 1.0 / a / a + 1.0 / b / b + 1.0 / c / c;
}

// generator pseudoslučajnih brojeva po uniformnoj raspodeli - svaka nit ima svoj seed, jer kada bi bio shared, uniformnost ne bi bila garantovana
// real 8-byte number in [0,1)
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
  int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)

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
        // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = ni -> x = a
        double x = ((double)(ni - i) * (-a) + (double)(i - 1) * a) / (double)(ni - 1);
        double y = ((double)(nj - j) * (-b) + (double)(j - 1) * b) / (double)(nj - 1);
        double z = ((double)(nk - k) * (-c) + (double)(k - 1) * c) / (double)(nk - 1);

        // check if point is inside elipsoid [a, b, c]
        double chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);

        if (chk > 1.0)
        {
          // tacka nije unutar elipsoida
// deo oznacen za DEBUG se kompajlira ako se prilikom prevodjenja navede opcija -DDEBUG
#ifdef DEBUG
          printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
                 x, y, z, 1.0, 1.0, 0.0, 0);
#endif

          continue;
        }

        // tacka je unutar elipsoida
        n_inside++;

        // analitička vrednost funkcije gustine/potencijala u tački unutar elipsoida - referentna vrednost koju poredimo u odnosu na numericku - wt
        double w_exact = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);
        double wt = 0.0;

#ifdef DEBUG
        int steps = 0;
#endif
        // pustamo N tacaka iz izabrane koordinate - visestruki pokusaji kako bi se dobila bolja aproksimacija
        for (int trial = 0; trial < N; trial++)
        {
          double x1 = x;
          double x2 = y;
          double x3 = z;
          double w = 1.0;
          chk = 0.0;
          // kretanje cestice - dok se nalazi unutar elipsoida
          while (chk < 1.0)
          {
            // da li se pomeramo ili ne (po sve tri ose) - sa verovatnocom 1/3 za svaku osu
            double ut = r8_uniform_01(&seed);

            // da li se pomeramo za +stepsz ili -stepsz
            double us;

            double dx;
            double dy;
            double dz;

            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
                dx = -stepsz;
              else
                dx = stepsz;
            }
            else
            {
              dx = 0.0;
            }
            
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
            {
              dy = 0.0;
            }
           
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
            {
              dz = 0.0;
            }
              
            // potential before moving
            double vs = potential(a, b, c, x1, x2, x3);

            // move
            x1 = x1 + dx;
            x2 = x2 + dy;
            x3 = x3 + dz;

#ifdef DEBUG
            ++steps;
#endif

            // potential after moving
            double vh = potential(a, b, c, x1, x2, x3);

            double we = (1.0 - h * vs) * w;           // Euler-ov korak
            w = w - 0.5 * h * (vh * we + vs * w);     // trapezna aproksimacija

            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
          }
          wt = wt + w;
        }
        // srednja vrenost tezine za N pokusaja
        wt = wt / (double)(N);

        // kvadrat razlike tacne i numericki dobijene vrednosti
        err += pow(w_exact - wt, 2);

#ifdef DEBUG
        printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
               x, y, z, wt, w_exact, fabs(w_exact - wt), steps / N);
#endif
      }
    }
  }
} // parallel
  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}

double feynman_2(const double a, const double b, const double c, const double h, const double stepsz, const int ni, const int nj, const int nk, const int N) 
{
  static int seed = 123456789;  // set to be static -> to be global (shared) by default
  int n_inside = 0;

  // for every point in grid
  double w_exact[17][12][7];
  double wt[17][12][7];

#pragma omp threadprivate(seed) // then, it is private for every thread (no race condition)
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
            double dy;
            double dz;
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
                dx = -stepsz;
              else
                dx = stepsz;
            }
            else
            {
              dx = 0.0;
            }


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
            {
              dy = 0.0;
            }
            
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
            {
              dz = 0.0;
            }

            // potential before moving
            double vs = potential(a, b, c, x1, x2, x3);

            // move
            x1 = x1 + dx;
            x2 = x2 + dy;
            x3 = x3 + dz;

            // potential after moving
            double vh = potential(a, b, c, x1, x2, x3);

            double we = (1.0 - h * vs) * w;
            w = w - 0.5 * h * (vh * we + vs * w);

            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
          }
// jer se menja u kontekstu vise taskova
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
                  // kada tacka nije unutar elipsoida
                  continue;
              }
              err += pow(w_exact[i][j][k] - (wt[i][j][k] / (double)(N)), 2);
          }
      }
  }
  return sqrt(err / (double)(n_inside));
}


double (*FUNCS[])(const double, const double, const double, const double, const double, const int, const int, const int, const int) = {feynman_1, feynman_2};


int main(int argc, char **argv)
{
  // velicina prostora po x, y i z
  const double a = 3.0;
  const double b = 2.0;
  const double c = 1.0;
  const double h = 0.001;       // vremenski interval između dva uzastopna stanja cestice tokom simulacije
  // broj tacaka po svakoj dimenziji
  const int ni = 16;
  const int nj = 11;
  const int nk = 6;

  const double stepsz = sqrt(3 * h);      // velicina koraka kojim se cestica pomera u svakom koraku

  if (argc < 3)
  {
    printf("Invalid number of arguments passed.\n");
    return 1;
  }

  // index of function
  const int func = atoi(argv[1]);
  
  // numer of walks per point
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




 /*
 a, b and represents elipsoid - potential energy is like springs in three directions with different stiffness
 (x * x) / (a * a) + (y * y) / (b * b) + (z * z) / (c * c) <= 1.0

 - cestica se kreci od tacke do tacke - po diskretizovanoj mrezi
 - iz svake tacke resetke, cestica krece N puta

 - w - trenutna tezina cestice u simulaciji
 */