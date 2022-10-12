#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NormRANu (2.3283063671E-10F)
#define PUNTOSMAX  250000
#define BLOQUESMAX 500
#define TMAX 100
#define HMIN 0.0001

//#define PRUEBA_BOX_MULLER
#define RUNGE_KUTTA_OSCILADOR
//#define RUNGE_KUTTA_DOBLEPOZO

int flag;
double pi = acos(-1);
double eta;
double kBT;
double k;
double B;
double m;
double h;
double const_estocastica;

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
        hist[i]/=(npuntos*delta);
        x[i]=minimo+i*delta;
    }
}

void exporta_histograma (int nbloques, double *x, double *hist, char *nombre)
{
    int i;
    FILE *f = fopen (nombre, "w");

    for (i=0; i<nbloques; i++)
    {
        fprintf (f, "%lf\t%lf\n", x[i], hist[i]);
    }
    fclose (f);
}

double potencial (int flag, double x)
{
    double aux;
    if (flag == 1) //oscilador armónico
        aux = 0.5*k*x*x;
    if (flag == 2)
        aux = B*(x*x-1)*(x*x-1);

    return aux;
}

double fuerza (int flag, double x)
{
    double aux;
    if (flag == 1) //oscilador armónico
        aux = -k*x;
    else if (flag == 2)
        aux = B*2*(x*x-1)*2*x;
    return aux;
}

double f(int flag, double p)
{
    return p/m;
}

double g(int flag, double x, double p)
{
    double aux;
    if (flag == 1)
        aux = fuerza(flag, x) -eta*p/m;
    return aux;
}

void RK (int flag, double *x, double *p)
{
    double fx1, fx2, gp1, gp2;
    double xn, pn, cte;
    double random1, random2;

    xn=*x; pn=*p;

    box_muller (&random1, &random2, const_estocastica);
    cte = sqrt(2*eta*kBT*h)*random1;
    //printf ("cte=%lf\n", cte);

    fx1 = f(flag, pn+cte);
    gp1 = g(flag, xn, pn+cte);
    //printf ("fx1=%lf\tgp1=%lf\n", fx1, gp1);
    fx2 = f(flag, pn+h*gp1);
    gp2 = g(flag, xn+h*fx1, pn+h*gp1);
    //printf ("fx2=%lf\tgp2=%lf\n", fx2, gp2);

    *x+=h/2.*(fx1+fx2);
    *p+=h/2.*(gp1+gp2)+cte;

    //printf ("x=%lf\tp=%lf\n", *x, *p);
}

double ecinetica (int flag, double p)
{
    return p*p/2./m;
}

double emecanica(int flag, double x, double p)
{
    return ecinetica (flag, p) + potencial (flag, x);
}

int main ()
{
    int i, t=TMAX, npuntos = PUNTOSMAX, nbloques=BLOQUESMAX, nmedidas;
    double g1, g2, med, var, varmed, ecin, epot, etot, sumaecin, sumaepot, sumaetot;
    double puntos [PUNTOSMAX];
    double hist [BLOQUESMAX], hist_ejex[BLOQUESMAX];
    double xaux, x[(int)(TMAX/HMIN)], histx [BLOQUESMAX], x_ejex[BLOQUESMAX];
    double paux, p[(int)(TMAX/HMIN)], histp [BLOQUESMAX], p_ejex[BLOQUESMAX];
    char *nombre;

    srand (time(NULL));

    #ifdef PRUEBA_BOX_MULLER
    const_estocastica = 2.5;
    for (i=0; i<npuntos/2; i++)
    {
        box_muller (&g1, &g2, const_estocastica);
        puntos [2*i]=g1; puntos [2*i+1]=g2;
        printf ("puntos generados=%d\n", 2*i);
    }

    med_var (puntos, npuntos, &med, &var, &varmed);
    printf ("Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);
    construye_histograma (nbloques, npuntos, puntos, hist, hist_ejex);
    nombre = "hist_boxmuller.txt";
    exporta_histograma (nbloques, hist_ejex, hist, nombre);
    #endif // PRUEBA_BOX_MULLER

    #ifdef RUNGE_KUTTA_OSCILADOR
    flag = 1;
    eta = 0; k=1; m=1; kBT = 1;
    const_estocastica = 2*eta*kBT;

    //Sin damping

    h=HMIN;
    nmedidas = t/h;
    xaux = 1; paux = -1;
    sumaecin = sumaepot = sumaetot = 0;
    for (i=0; i<nmedidas; i++)
    {
        if (i%(nmedidas/10)==0) printf ("eta=%lf\th=%lf\ti=%d\n", eta, h, i);
        RK (flag, &xaux, &paux);
        x[i]=xaux; p[i]=paux;
        //ecin = ecinetica (flag, p[i]); epot = potencial(flag, x[i]); etot = emecanica (flag, x[i], p[i]);
        //printf ("Ecin=%lf\tEpot=%lf\tEtot=%lf\n\n", ecin, epot, etot);
        //sumaecin += ecin; sumaepot +=epot; sumaetot +=etot;
    }
    //sumaecin /= (double) nmedidas; sumaepot /= (double) nmedidas; sumaetot /= (double) nmedidas;
    //printf ("<Ecin>=%lf\t<Epot>=%lf\t<Etot>%lf\n", sumaecin, sumaepot, sumaetot);

    construye_histograma (nbloques, nmedidas, x, histx, x_ejex);
    construye_histograma (nbloques, nmedidas, p, histp, p_ejex);

    printf ("Calculando <Ecin>, <Epot> y <Etot> a partir de los valores esperados de <p^2> y <x^2>\n");
    med_var (x, nmedidas, &med, &var, &varmed);
    printf ("Para las x: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);
    printf ("<Epot>=%lf\n", 0.5*k*(var+med*med));
    med_var (p, nmedidas, &med, &var, &varmed);
    printf ("Para las p: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);
    printf ("<Ecin>=%lf\n", 0.5/m*(var+med*med));

    sprintf (nombre, "hist_x_eta=%lf_h=%lf.txt", eta, h);
    exporta_histograma (nbloques, x_ejex, histx, nombre);
    sprintf (nombre, "hist_p_eta=%lf_h=%lf.txt", eta, h);
    exporta_histograma (nbloques, p_ejex, histp, nombre);

    /**Con damping**/
    eta = 0.1;
    while (eta<=10)
    {
        h =HMIN;
        while (h<=0.1)
        {
            nmedidas = t/h;
            xaux = 1; paux = -1;
            sumaecin = sumaepot = sumaetot = 0;
            for (i=0; i<nmedidas; i++)
            {
                if (i%(nmedidas/10)==0) printf ("eta=%lf\th=%lf\ti=%d\n", eta, h, i);
                RK (flag, &xaux, &paux);
                x[i]=xaux; p[i]=paux;
                //ecin = ecinetica (flag, p[i]); epot = potencial(flag, x[i]); etot = emecanica (flag, x[i], p[i]);
                //printf ("Ecin=%lf\tEpot=%lf\tEtot=%lf\n\n", ecin, epot, etot);
                //sumaecin += ecin; sumaepot +=epot; sumaetot +=etot;
            }

            //sumaecin /= (double) nmedidas; sumaepot /= (double) nmedidas; sumaetot /= (double) nmedidas;
            //printf ("<Ecin>=%lf\t<Epot>=%lf\t<Etot>%lf\n", sumaecin, sumaepot, sumaetot);

            construye_histograma (nbloques, nmedidas, x, histx, x_ejex);
            construye_histograma (nbloques, nmedidas, p, histp, p_ejex);

            printf ("Calculando <Ecin>, <Epot> y <Etot> a partir de los valores esperados de <p^2> y <x^2>\n");
            med_var (x, nmedidas, &med, &var, &varmed);
            printf ("Para las x: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);
            printf ("<Epot>=%lf\n", 0.5*k*(var+med*med));
            med_var (p, nmedidas, &med, &var, &varmed);
            printf ("Para las p: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);
            printf ("<Ecin>=%lf\n", 0.5/m*(var+med*med));

            sprintf (nombre, "hist_x_eta=%lf_h=%lf.txt", eta, h);
            exporta_histograma (nbloques, x_ejex, histx, nombre);
            sprintf (nombre, "hist_p_eta=%lf_h=%lf.txt", eta, h);
            exporta_histograma (nbloques, p_ejex, histp, nombre);

            h=10*h;
            printf ("\n\n");
        }
        eta=10*eta;
    }
    #endif // RUNGE_KUTTA_OSCILADOR

    #ifdef RUNGE_KUTTA_DOBLEPOZO /**Estamos trabajandolo aun**/
    flag = 2;
    eta = 0.1; B=0.5; m=1; kBT = 0.2;
    const_estocastica = 2*eta*kBT;

    /**Sin damping**/
    //eta = 0;

    /**Con damping variable**/
    while (eta<=10)
    {
        h =HMIN*10;
        while (B<=2)
        {
            nmedidas = t/h;
            xaux = 0; paux = 2*(parisi_rapuano()-0.5);
            sumaecin = sumaepot = sumaetot = 0;
            for (i=0; i<nmedidas; i++)
            {
                if (i%(nmedidas/10)==0) printf ("eta=%lf\tB=%lf\ti=%d\n", eta, B, i);
                RK (flag, &xaux, &paux);
                x[i]=xaux; p[i]=paux;
                ecin = ecinetica (flag, p[i]); epot = potencial(flag, x[i]); etot = emecanica (flag, x[i], p[i]);
                //printf ("Ecin=%lf\tEpot=%lf\tEtot=%lf\n\n", ecin, epot, etot);
                sumaecin += ecin; sumaepot +=epot; sumaetot +=etot;
            }

            sumaecin /= (double) nmedidas; sumaepot /= (double) nmedidas; sumaetot /= (double) nmedidas;
            printf ("<Ecin>=%lf\t<Epot>=%lf\t<Etot>%lf\n", sumaecin, sumaepot, sumaetot);

            construye_histograma (nbloques, nmedidas, x, histx, x_ejex);
            construye_histograma (nbloques, nmedidas, p, histp, p_ejex);

            med_var (x, nmedidas, &med, &var, &varmed);
            printf ("Para las x: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);
            med_var (p, nmedidas, &med, &var, &varmed);
            printf ("Para las p: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);

            //sprintf (nombre, "hist_x_eta=%lf_h=%lf.txt", eta, h);
            //exporta_histograma (nbloques, x_ejex, histx, nombre);
            //sprintf (nombre, "hist_p_eta=%lf_h=%lf.txt", eta, h);
            //exporta_histograma (nbloques, p_ejex, histp, nombre);

            B=2*B;
            printf ("\n\n");
        }
        eta=10*eta;
    }

    #endif // RUNGE_KUTTA_DOBLEPOZO


    return 0;
}

