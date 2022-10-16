#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NormRANu (2.3283063671E-10F)
#define PUNTOSMAX  250000
#define BLOQUESMAX 500
#define MEDIDASMAX 100000

//#define PRUEBA_BOX_MULLER
#define RUNGE_KUTTA
//#define Verlet
//#define EULER_MARUYAMA

double pi = acos(-1);
double eta = 1;
double kBT = 1;
double k = 1;
double m = 1;
double h = 0.1;
double B = 0.5;
double const_estocastica = 2;

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
        suma_var += input[i]*input[i];
    }

    *varianza = suma_var/numero;
    *varianza_media = *varianza/numero;
}

void box_muller (double *g1, double *g2) //LO cambiamos a generador normal
{
    double w1, w2;
    w1 = parisi_rapuano ();
    w2 = parisi_rapuano ();

    *g1 = -sqrt(-2*log(w1))*cos(2*pi*w2);
    *g2 = -sqrt(-2*log(w1))*sin(2*pi*w2);
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

void mostrarGraficas(char *n1, char *n2){
    FILE *f;
    f=fopen("gnu.plt","w");
    fprintf(f,"set multiplot layout 1,2 rowsfirst \n plot '%s' u 1:2\n plot '%s' u 1:2 \n unset multiplot",n1,n2);
    fclose(f);
    system("gnu.plt");
}

double potencial (int flag, double x)
{
    double aux;
    if (flag == 1)
        aux = 0.5*k*x*x;
    if (flag == 2)
        aux = B*(x*x-1)*(x*x-1);
    return aux;
}

double cinetica (double p)
{
    return p*p/m;
}

double mecanica (int flag, double x, double p)
{
    double pot, cine;
    pot = potencial (flag, x);
    cine = cinetica (p);
    return pot + cine;
}

double fuerza (int flag, double x)
{
    double aux;
    if (flag == 1) //Oscilador armónico
        aux = -k*x;
    if (flag == 2)
        aux = -2*B*2*x*(x*x-1);
    return aux;
}

double f(double p)
{
    return p/m;
}

double g(int flag, double x, double p)
{
    return (fuerza(flag, x) - eta/m*p);
}

void RK (int flag, double *x, double *p, int i, double *random)
{
    double fx1, fx2, gp1, gp2;
    double xn, pn, cte;
    double random1, random2;

    xn=*x; pn=*p;

    if (i==0)
    {
        box_muller (&random1, &random2);
        *random = random2;
    }
    if (i==1)
    {
        random1 = *random;
    }

    cte = sqrt(2*eta*kBT*h)*random1;
    //printf ("cte=%lf\n", cte);

    fx1 = f(pn+cte);
    gp1 = g(flag, xn, pn+cte);
    //printf ("fx1=%lf\tgp1=%lf\n", fx1, gp1);
    fx2 = f(pn+h*gp1);
    gp2 = g(flag, xn+h*fx1, pn+h*gp1);
    //printf ("fx2=%lf\tgp2=%lf\n", fx2, gp2);

    *x+=h/2.*(fx1+fx2);
    *p+=h/2.*(gp1+gp2)+cte;

    //printf ("x=%lf\tp=%lf\n\n", *x, *p);
}

void verlet (int flag, double *x, double *p, double a, double b, int i, double *random)
{

    double random1, random2;
    double xo, po, vo, vn;
    xo=*x;
    po=*p;
    vo=po/m;
    if (i==0)
    {
        box_muller (&random1, &random2);
        *random = random2;
    }
    if (i==1)
    {
        random1 = *random;
    }

    random1 = sqrt(2*eta*kBT*h)*random1;
    *x = xo + b*h*vo + b*h*h/(2*m)*fuerza(flag, xo) + b*h/(2*m)*random1;
    vn=a*vo+h/(2*m)*(a*f(xo)+fuerza(flag, *x))+b/m*random1;
    *p=m*vn;

}

void eulerMaruyama(int flag, double *x,double *p,int i, double *random)
{
    double random1, random2;
    double xn, pn;

    xn=*x; pn=*p;


    if (i==0)
    {
        box_muller (&random1, &random2);
        *random = random2;
    }
    if (i==1)
    {
        random1 = *random;
    }

    random1 = sqrt(2*eta*kBT*h)*random1;

    *x = xn + f(pn)*h;
    *p = pn + g(flag, xn, pn) + random1;

    if(i%1000==0){printf ("x=%lf\tp=%lf\n\n", *x, *p);}

}

int distribucion_tiempos (int n, double *input, double *t)
{
    int i, contador, index_t;
    double signo;

    signo = input[0];
    contador = 1; index_t=0;
    for (i=0; i<n; i++)
    {
        if (input[i]*signo > 0)
            contador ++;
        else
        {
            t[index_t]=contador;
            contador = 1; index_t++;
        }
    }

    return index_t;
}

int main ()
{
    int i, flag, ntermalizacion = 100000, npuntos = PUNTOSMAX, nbloques=100, nmedidas=MEDIDASMAX, ntiempos;
    double g1, g2, med, var, varmed, random;
    double puntos [PUNTOSMAX];
    double hist [BLOQUESMAX], hist_ejex[BLOQUESMAX];
    double xaux, x[MEDIDASMAX], histx [BLOQUESMAX], x_ejex[BLOQUESMAX];
    double paux, p[MEDIDASMAX], histp [BLOQUESMAX], p_ejex[BLOQUESMAX];
    double tiempos [MEDIDASMAX], histt[BLOQUESMAX], t_ejex [BLOQUESMAX];
    char *nombre1,*nombre2, *nombre3;

    flag = 2; //1 para oscilador armónico, 2 para doble pozo

    #ifdef Verlet
    double a,b;
    a=(1-eta*h/(2*m))/(1+eta*h/(2*m));
    b=1/(1+eta*h/(2*m));
    #endif // Verlet


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
    nombre = "hist_boxmuller.txt";
    exporta_histograma (nbloques, hist_ejex, hist, nombre);
    #endif // PRUEBA_BOX_MULLER


    xaux = 0; paux = 0;
    for (i=0; i<ntermalizacion; i++)
    {
        if(i%1000==0){printf ("i=%d\n", i);}
        #ifdef RUNGE_KUTTA
        RK (flag, &xaux, &paux, i%2, &random);
        #endif // RUNGE_KUTTA
        #ifdef EULER_MARUYAMA
        eulerMaruyama(flag, &xaux,&paux,(i%2), &random);
        #endif // EULER_MARUYAMA
        #ifdef Verlet
        verlet (&xaux, &paux, a, b);
        #endif // Verlet
    }
    for (i=0; i<nmedidas; i++)
    {
        if(i%1000==0){printf ("i=%d\n", i);}
        #ifdef RUNGE_KUTTA
        RK (flag, &xaux, &paux, i%2, &random);
        #endif // RUNGE_KUTTA
        #ifdef EULER_MARUYAMA
        eulerMaruyama(flag, &xaux,&paux,(i%2), &random);
        #endif // EULER_MARUYAMA
        #ifdef Verlet
        verlet (&xaux, &paux, a, b);
        #endif // Verlet
        x[i]=xaux; p[i]=paux;
    }
    construye_histograma (nbloques, nmedidas, x, histx, x_ejex);
    construye_histograma (nbloques, nmedidas, p, histp, p_ejex);

    ntiempos = distribucion_tiempos (nmedidas, x, tiempos);
    construye_histograma (nbloques, ntiempos, tiempos, histt, t_ejex);

    med_var (x, nmedidas, &med, &var, &varmed);
    printf ("Para las x: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);
    med_var (p, nmedidas, &med, &var, &varmed);
    printf ("Para las p: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);

    nombre1 = "hist_x.txt";
    exporta_histograma (nbloques, x_ejex, histx, nombre1);
    nombre2 = "hist_p.txt";
    exporta_histograma (nbloques, p_ejex, histp, nombre2);
    nombre3 = "hist_t.txt";
    exporta_histograma (nbloques, t_ejex, histt, nombre3);

    mostrarGraficas(nombre1,nombre3);

}
