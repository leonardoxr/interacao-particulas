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
    double new_aceleration[2];
    double old_aceleration[2];

};

double dist(double x_1,double x_2,double y_1,double y_2);
double update(struct particulas particula[],double dt,int L,int maxpart);
double a(double r,double R);

void main()
{

    FILE *arq;
    arq = fopen("particulas.txt", "w+");
    srand48(time(NULL));

    double u1,u2,ecm,t,dt,dx,dy,i,k;
    int j,L,maxpart;

    L=10;
    dx = 1;
    dy = 1;

    maxpart = pow((L/dx),2);

    dt = 0.01;

   // printf("%d\n", maxpart);
    struct particulas particula[maxpart];


//CAIXA 10X10
    j=0;
    for(i=0; i<L; i+=dy)
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


            particula[j].new_aceleration[0] = 0;
            particula[j].new_aceleration[1] = 0;
            particula[j].old_aceleration[0] = 0;
            particula[j].old_aceleration[1] = 0;
            j++;
        }


    }




    for(t=0; t<100; t+=dt)
    {
        for(j=0; j<maxpart; j++)
        {
            fprintf(arq,"%lf %lf\n", particula[j].position[0],particula[j].position[1]);
        }
        fprintf(arq,"\n\n");



        update(particula,dt,L,maxpart);



    }




    for(j=0; j<maxpart; j++)
    {
        ecm += 1.*pow(particula[j].velocity[0],2)/2. + 1.*pow(particula[j].velocity[1],2)/2.;
    }
    ecm = ecm/maxpart;
    //printf("%lf\n", ecm);




    fclose(arq);


}



double update(struct particulas particula[],double dt,int L,int maxpart)
{

    double r[maxpart][2],R[maxpart],x1temp,x2temp;
    int i,j;


    for(i=0; i<maxpart; i++)
    {
        particula[i].old_aceleration[0] = 0;
        particula[i].old_aceleration[1] = 0;
        particula[i].new_aceleration[0] = 0;
        particula[i].new_aceleration[1] = 0;



        for(j = 0; j<maxpart; j++)
        {
            R[j] = dist(particula[i].position[0],particula[j].position[0],particula[i].position[1],particula[j].position[1]);
            r[j][0] = particula[i].position[0] - particula[j].position[0];
            r[j][1] = particula[i].position[1] - particula[j].position[1];


        }



        for(j = 0; j<maxpart; j++)
        {

            particula[i].old_aceleration[0] += a(r[j][0],R[j]);
            particula[i].old_aceleration[1] += a(r[j][1],R[j]);





        }
       // printf("%lf\n", particula[i].old_aceleration[0] );



       x1temp = particula[i].position[0] + particula[i].velocity[0]*dt + 0.5 * particula[i].old_aceleration[0] * pow(dt,2);
       x2temp = particula[i].position[1] + particula[i].velocity[1]*dt + 0.5 * particula[i].old_aceleration[1] * pow(dt,2);


    if(x1temp > L - 1 || x1temp < 0)
    {
        particula[i].velocity[0] = - particula[i].velocity[0];
        particula[i].old_aceleration[0] = - particula[i].old_aceleration[0];

        x1temp =   particula[i].position[0] + particula[i].velocity[0]*dt + 0.5 * particula[i].old_aceleration[0] * pow(dt,2);
    }

    if(x2temp > L - 1 || x2temp < 0)
    {
        particula[i].velocity[1] = - particula[i].velocity[1];
        particula[i].old_aceleration[1] = - particula[i].old_aceleration[1];

        x2temp = particula[i].position[1] + particula[i].velocity[1]*dt + 0.5 * particula[i].old_aceleration[1] * pow(dt,2);

    }


        particula[i].position[0] = x1temp;
        particula[i].position[1] = x2temp;





        for(j = 0; j<maxpart; j++)
        {
            R[j] = dist(particula[i].position[0],particula[j].position[0],particula[i].position[1],particula[j].position[1]);
            r[j][0] = particula[i].position[0] - particula[j].position[0];
            r[j][1] = particula[i].position[1] - particula[j].position[1];


        }

        for(j = 0; j<maxpart; j++)
        {

            particula[i].new_aceleration[0] += a(r[j][0],R[j]);
            particula[i].new_aceleration[1] += a(r[j][1],R[j]);


        }



      particula[i].velocity[0] = particula[i].velocity[0] + 0.5 * (particula[i].new_aceleration[0] + particula[i].old_aceleration[0]) * dt;
      particula[i].velocity[1] = particula[i].velocity[1] + 0.5 * (particula[i].new_aceleration[1] + particula[i].old_aceleration[1]) * dt;


    }
}


double dist(double x_1,double x_2,double y_1,double y_2)
{


    double dist;

    dist = sqrt(pow((x_1 - x_2),2)  + pow((y_1 - y_2),2));
    return dist;


}


double a(double r,double R)
{

    double ac,e,sigma,A;
    e = 1;
    sigma = 0.4 ;
    if(fabs(R) < 10e-7)
    {
        return 0;
    }
    else
    {

    A = 48.*e*(pow((sigma/R),14) -  0.5*pow((sigma/R),8));

    ac = A*r;
     //printf("%lf\n", A);
    return ac;
    }

}
