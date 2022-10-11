#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NormRANu (2.3283063671E-10F)
#define PUNTOSMAX  250000
#define BLOQUESMAX 500
#define MEDIDASMAX 100000
#define const_estocastica  2*eta*kBT

#define RESONANCIA 0 //fuerza periodica sí o no

//#define PRUEBA_BOX_MULLER
//#define RUNGE_KUTTA
//#define Verlet
#define EULER_MARUYAMA

double pi = acos(-1);
double eta = 1;
double kBT =0.05;
double k = 1;
double m = 1;
double h = 0.1;
double B=100;
double A=1; //constante fuerza periodica
double w=1; //frecuencia fuerza periodica

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

void RK (double *x, double *p)
{
    double fx1, fx2, gp1, gp2;
    double xn, pn, cte;
    double random1, random2;

    xn=*x; pn=*p;

    box_muller (&random1, &random2, const_estocastica);
    cte = sqrt(2*eta*kBT*h)*random1;
    //printf ("cte=%lf\n", cte);

    fx1 = (pn+cte)/m;
    gp1 = (-eta/m*(pn+cte))-k*xn;
    //printf ("fx1=%lf\tgp1=%lf\n", fx1, gp1);
    fx2 = (pn+h*gp1)/m;
    gp2 = -eta/m*(pn+h*gp1)-k*(xn+h*fx1);
    //printf ("fx2=%lf\tgp2=%lf\n", fx2, gp2);

    *x+=h/2.*(fx1+fx2);
    *p+=h/2.*(gp1+gp2)+cte;

    printf ("x=%lf\tp=%lf\n\n", *x, *p);
}

double f(double x){
    return -k*x;
}

void numeroOcupacion(int n[2],double x){
    if(x>0){n[1]++;}else{n[0]++;}
}

void verlet (double *x, double *p, double a, double b){

    double random1, random2;
    double xo, po, vo, vn;
    xo=*x;
    po=*p;
    vo=po/m;
    box_muller (&random1, &random2, const_estocastica);
    *x=xo+b*h*vo+b*h*h/(2*m)*f(xo)+b*h/(2*m)*random1;
    vn=a*vo+h/(2*m)*(a*f(xo)+f(*x))+b/m*random1;
    *p=m*vn;

}

void eulerMaruyama(double *x,double *p,int i,int resonancia, double t){
    double random1, random2;


    box_muller(&random1, &random2, const_estocastica);

    *x=*x+h*(*p)/m;
    if(!resonancia){
        *p=*p+(-eta*(*p)/m-4/B/m*(*x)*((*x)*(*x)-1))*h+sqrt(2*eta*kBT*h)*random1;
    }else{
        *p=*p+(-eta*(*p)/m-4/B/m*(*x)*((*x)*(*x)-1)+A*sin(w*t))*h+sqrt(2*eta*kBT*h)*random1;
    }
    if(i%1000==0){printf ("x=%lf\tp=%lf\n\n", *x, *p);}

}

int main ()
{
    int i, npuntos = PUNTOSMAX, nbloques=BLOQUESMAX, nmedidas=MEDIDASMAX;
    double g1, g2, med, var, varmed;
    double puntos [PUNTOSMAX];
    double hist [BLOQUESMAX], hist_ejex[BLOQUESMAX];
    double xaux, x[MEDIDASMAX], histx [BLOQUESMAX], x_ejex[BLOQUESMAX];
    double paux, p[MEDIDASMAX], histp [BLOQUESMAX], p_ejex[BLOQUESMAX];
    char *nombre1,*nombre2;
    int num_ocupacion[2];
    double t=0; //necesaria para fuerza periodica

    num_ocupacion[0]=num_ocupacion[1]=0;
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


    xaux = 1; paux = -1;
    for (i=0; i<nmedidas; i++)
    {
        if(i%1000==0){printf ("i=%d\n", i);}
        #ifdef RUNGE_KUTTA
        RK (&xaux, &paux);
        #endif // RUNGE_KUTTA
        #ifdef EULER_MARUYAMA
        eulerMaruyama(&xaux,&paux,i,RESONANCIA,t);
        #endif // EULER_MARUYAMA
        #ifdef Verlet
        verlet (&xaux, &paux, a, b);
        #endif // Verlet
        x[i]=xaux; p[i]=paux;
        numeroOcupacion(num_ocupacion,xaux);
        t+=h;
    }
    construye_histograma (nbloques, nmedidas, x, histx, x_ejex);
    construye_histograma (nbloques, nmedidas, p, histp, p_ejex);

    med_var (x, nmedidas, &med, &var, &varmed);
    printf ("Para las x: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);
    med_var (p, nmedidas, &med, &var, &varmed);
    printf ("Para las p: Media:%lf\tVarianza:%lf\tVarianzaMedia:%lf\n", med, var, varmed);

    printf("Num. de ocupacion x<0: %d \tp Num. de ocupacion x>0: %d",num_ocupacion[0],num_ocupacion[1]);

    nombre1 = "hist_x.txt";
    exporta_histograma (nbloques, x_ejex, histx, nombre1);
    nombre2 = "hist_p.txt";
    exporta_histograma (nbloques, p_ejex, histp, nombre2);

    mostrarGraficas(nombre1,nombre2);

}

