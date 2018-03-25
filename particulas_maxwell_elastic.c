#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>


// srand48(time(NULL));

//r = drand48();

struct particulas
{
    double position[2];
    double velocity[2];

};

double mru(struct particulas *particula,double dt,int L);

void main()
{

    FILE *arq;
    arq = fopen("particulas.txt", "w+");
    srand48(time(NULL));

    double u1,u2,ecm,t,dt,dx,dy,i,k;
    int j,L,maxpart;


    L=100;
    dx = 10;

    maxpart = pow((L/dx),2) ;

    dt = 0.1;

    printf("%d\n", maxpart);
    struct particulas particula[maxpart];
//CAIXA 10X10
    j=0;
    for(i=0; i<L; i+=dx)
    {


        for(k=0; k<L; k+=dx)
        {



            // x1 = sqrt(-2*log(drand48()))*cos(2*M_PI*drand48())

            u1 = drand48();
            u2 = drand48();
            particula[j].position[0] = k;
            particula[j].velocity[0] = sqrt(-2*log(u1))*cos(2*M_PI*u1);


            particula[j].position[1] = i;
            particula[j].velocity[1] = sqrt(-2*log(u2))*sin(2*M_PI*u2);
            j++;
        }


    }




    for(t=0; t<100; t+=dt)
    {
        for(j=0; j<maxpart; j++)
        {
            fprintf(arq,"%lf %lf %lf %lf\n", particula[j].position[0],particula[j].velocity[0],particula[j].position[1],particula[j].velocity[1]);
        }
        fprintf(arq,"\n\n");



        for(j=0; j<maxpart; j++)
        {
            mru(&particula[j],dt,L);

        }

    }




    for(j=0; j<maxpart; j++)
    {
        ecm += 1.*pow(particula[j].velocity[0],2)/2. + 1.*pow(particula[j].velocity[1],2)/2.;
    }
    ecm = ecm/maxpart;
    printf("%lf\n", ecm);




    fclose(arq);


}



double mru(struct particulas *particula,double dt,int L)
{
    double x1temp, x2temp;

    x1temp =  particula->position[0] + particula->velocity[0]*dt;
    x2temp = particula->position[1] + particula->velocity[1]*dt;

    if(x1temp > L - 1 || x1temp < 0)
    {
        particula->velocity[0] = - particula->velocity[0];
        x1temp =  particula->position[0] + particula->velocity[0]*dt;
    }

    if(x2temp > L - 1 || x2temp < 0)
    {
        particula->velocity[1] = - particula->velocity[1];
        x2temp = particula->position[1] + particula->velocity[1]*dt;

    }


    particula->position[0] = x1temp;
    particula->position[1] = x2temp;


}
