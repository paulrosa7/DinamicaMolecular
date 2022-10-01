#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NormRANu (2.3283063671E-10F)
#define PUNTOSMAX  250000
#define BLOQUESMAX 500
#define MEDIDASMAX 100000

#define PRUEBA_BOX_MULLER

double pi = acos(-1);
double eta = 1;
double kBT = 1;
double k = 1;
double m = 1;
double h = 0.0001;
double const_estocastica = 2*eta*kBT;

double parisi_rapuano ()
{
    unsigned int irr [256];
    unsigned char ind_ran=0, ig1=0, ig2=0, ig3=0;
    double numero_aleatorio;
    int i;

    for (i=0; i<256; i++)
        irr [i] = (rand () << 16) + rand ();

    ig1 = ind_ran - 24;
    ig2 = ind_ran - 55;
    ig3 = ind_ran - 61;
    irr [ind_ran] = irr [ig1] + irr [ig2];

    numero_aleatorio = (irr[ind_ran]^irr [ig3]);
    ind_ran++;

    return numero_aleatorio*NormRANu;
}

void med_var (double *input, int numero, double *media, double *varianza, double *varianza_media)
{
    int i;
    double suma_med = 0, suma_var = 0;
    for (i = 0; i<numero; i++)
    {
        suma_med += input [i];
    }

    *media = suma_med/numero;

    for (i=0; i<numero; i++)
    {
        suma_var += pow (input[i]-*media, 2);
    }

    *varianza = suma_var/(numero-1);
    *varianza_media = *varianza/numero;
}

void box_muller (double *g1, double *g2, double varianza)
{
    double w1, w2;
    w1 = parisi_rapuano ();
    w2 = parisi_rapuano ();

    *g1 = -sqrt(-2*log(w1))*cos(2*pi*w2)*sqrt(varianza);
    *g2 = -sqrt(-2*log(w1))*sin(2*pi*w2)*sqrt(varianza);
}

void min_max (int npuntos, double *puntos, double *minimo, double *maximo)
{
    int i;

    *minimo = *maximo = puntos [0];
    for (i=1; i<npuntos; i++)
    {
        if (puntos[i]<*minimo) *minimo=puntos[i];
        if (puntos[i]>*maximo) *maximo=puntos[i];
    }
}

void construye_histograma (int nbloques, int npuntos, double *puntos, double *hist, double *x)
{
    double minimo, maximo, delta;
    int i, aux;

    min_max (npuntos, puntos, &minimo, &maximo);
    delta = (maximo-minimo)/((double)(nbloques));
    printf ("minimo=%lf\tmaximo=%lf\tdelta=%lf\n", minimo, maximo, delta);

    for (i=0; i<npuntos; i++)
    {
        aux = (int)((puntos[i]-minimo)/delta);
        if (aux==nbloques) aux--;
        hist[aux] ++;
    }

    for (i=0; i<nbloques; i++)
    {
        hist[i]/=npuntos;
        x[i]=minimo+i*delta;
    }


}

void exporta_gaussiana (int nbloques, double *x, double *hist)
{
    int i;
    FILE *f = fopen ("hist_gaussiana.txt", "w");

    for (i=0; i<nbloques; i++)
    {
        fprintf (f, "%lf\t%lf\n", x[i], hist[i]);
    }
    fclose (f);
}

int main ()
{
    int i, npuntos = PUNTOSMAX, nbloques=BLOQUESMAX, nmedidas=MEDIDASMAX;
    double g1, g2, med, var, varmed;
    double puntos [PUNTOSMAX];
    double hist [BLOQUESMAX], hist_ejex[BLOQUESMAX];
    double xaux, x[MEDIDASMAX], histx [BLOQUESMAX], x_ejex[BLOQUESMAX];
    double paux, p[MEDIDASMAX], histp [BLOQUESMAX], p_ejex[BLOQUESMAX];

    srand (time(NULL));

    #ifdef PRUEBA_BOX_MULLER
    for (i=0; i<npuntos/2; i++)
    {
        box_muller (&g1, &g2, const_estocastica);
        puntos [2*i]=g1; puntos [2*i+1]=g2;
        printf ("puntos generados=%d\n", 2*i);
    }

    med_var (puntos, npuntos, &med, &var, &varmed);
    printf ("Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);
    construye_histograma (nbloques, npuntos, puntos, hist, hist_ejex);
    exporta_gaussiana (nbloques, hist_ejex, hist);
    #endif // PRUEBA_BOX_MULLER

    return 0;
}
