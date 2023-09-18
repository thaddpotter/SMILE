
#include <smile.h>
//Does the initial time shift for orbits
//Takes a set of ri, r'i and increments them through time up to ti
void shift_time(double *x1_0, double *y1_0, double *vx1_0,double *vy1_0,
            double *x2_0, double *y2_0, double *vx2_0,double *vy2_0,
            float ti, float m1, float m2){
    
    double dt0 = ti*sday/1000;

    //Calculate common denominator for r''[0]
    double denom = pow((*x1_0-*x2_0)*(*x1_0-*x2_0) + (*y1_0-*y2_0)*(*y1_0-*y2_0), 1.5);

    //Get 4 (r'')'s
    double ax1 = -(bigg * m2 * (*x1_0-*x2_0))/denom;
    double ay1 = -(bigg * m2 * (*y1_0-*y2_0))/denom;
    double ax2 = (bigg * m1 * (*x1_0-*x2_0))/denom;
    double ay2 = (bigg * m1 * (*y1_0-*y2_0))/denom;

    //Increment r'(i+1/2)
    *vx1_0 += ax1*dt0/2;
    *vy1_0 += ay1*dt0/2;
    *vx2_0 += ax2*dt0/2;
    *vy2_0 += ay2*dt0/2;

    //Increment r(i+1)
    *x1_0 += (*vx1_0)*dt0;
    *y1_0 += (*vy1_0)*dt0;
    *x2_0 += (*vx2_0)*dt0;
    *y2_0 += (*vy2_0)*dt0;

    //Loop through time (Partially unrolled to prevent excessive r'' calcs)
    //-------------------------------------------
    for(int i=1;i<(1000);i++){
        
        //Calculate common denominator for r''[i]
        denom = pow((*x1_0-*x2_0)*(*x1_0-*x2_0) + (*y1_0-*y2_0)*(*y1_0-*y2_0), 1.5);
        
        //Get 4 (r'')'s
        ax1 = -(bigg * m2 * (*x1_0-*x2_0))/denom;
        ay1 = -(bigg * m2 * (*y1_0-*y2_0))/denom;
        ax2 = (bigg * m1 * (*x1_0-*x2_0))/denom;
        ay2 = (bigg * m1 * (*y1_0-*y2_0))/denom;

        //Increment r'(i-1/2) -> r'(i+1/2)
        *vx1_0 += ax1*dt0;
        *vy1_0 += ay1*dt0;
        *vx2_0 += ax2*dt0;
        *vy2_0 += ay2*dt0;

        //Increment r(i) -> r(i+1)
        *x1_0 += (*vx1_0)*dt0;
        *y1_0 += (*vy1_0)*dt0;
        *x2_0 += (*vx2_0)*dt0;
        *y2_0 += (*vy2_0)*dt0;
    }

    //Get 4 (r'')'s
    ax1 = -(bigg * m2 * (*x1_0-*x2_0))/denom;
    ay1 = -(bigg * m2 * (*y1_0-*y2_0))/denom;
    ax2 = (bigg * m1 * (*x1_0-*x2_0))/denom;
    ay2 = (bigg * m1 * (*y1_0-*y2_0))/denom;

    //Increment r'(i-1/2) -> r'(i)
    *vx1_0 += ax1*dt0/2;
    *vy1_0 += ay1*dt0/2;
    *vx2_0 += ax2*dt0/2;
    *vy2_0 += ay2*dt0/2;

}


//Given some set of initial conditions, return vector pairs for r, v
//Uses Kick Drift Kick method
void orbit(double x1_0, double y1_0, double vx1_0,double vy1_0,
            double x2_0, double y2_0, double vx2_0,double vy2_0,
            double x1[2048], double y1[2048], double vx1[2048],double vy1[2048],
            double x2[2048], double y2[2048], double vx2[2048],double vy2[2048], double dt, float m1, float m2){
    
    //Initial Values
    x1[0] = x1_0;
    y1[0] = y1_0;
    vx1[0] = vx1_0;
    vy1[0] = vy1_0;

    x2[0] = x2_0;
    y2[0] = y2_0;
    vx2[0] = vx2_0;
    vy2[0] = vy2_0;

    //Calculate common denominator for r''[0]
    double denom = pow( pow((x1[0]-x2[0]),2) + pow((y1[0]-y2[0]),2), 1.5);

    //Get 4 (r'')'s
    double ax1 = -(bigg * m2 * (x1[0]-x2[0]))/denom;
    double ay1 = -(bigg * m2 * (y1[0]-y2[0]))/denom;
    double ax2 = (bigg * m1 * (x1[0]-x2[0]))/denom;
    double ay2 = (bigg * m1 * (y1[0]-y2[0]))/denom;

    //Increment r'(i+1/2)
    double vshiftx1 = vx1[0] + ax1*dt/2;
    double vshifty1 = vy1[0] + ay1*dt/2;
    double vshiftx2 = vx2[0] + ax2*dt/2;
    double vshifty2 = vy2[0] + ay2*dt/2;

    //Increment r(i+1)
    x1[1] = x1_0 + vshiftx1*dt;
    y1[1] = y1_0 + vshifty1*dt;
    x2[1] = x2_0 + vshiftx2*dt;
    y2[1] = y2_0 + vshifty2*dt;

    //Loop through time (Partially unrolled to prevent excessive r'' calcs)
    //-------------------------------------------
    int i;
    for(i=1;i<(2048-1);i++){
        
        //Calculate common denominator for r''[i]
        denom = pow( pow((x1[i]-x2[i]),2) + pow((y1[i]-y2[i]),2), 1.5);
        
        //Get 4 (r'')'s
        ax1 = -(bigg * m2 * (x1[i]-x2[i]))/denom;
        ay1 = -(bigg * m2 * (y1[i]-y2[i]))/denom;
        ax2 = (bigg * m1 * (x1[i]-x2[i]))/denom;
        ay2 = (bigg * m1 * (y1[i]-y2[i]))/denom;

        //Increment r'(i-1/2) -> r'(i)
        vx1[i]= vshiftx1 + ax1*dt/2;
        vy1[i]= vshifty1 + ay1*dt/2;
        vx2[i]= vshiftx2 + ax2*dt/2;
        vy2[i]= vshifty2 + ay2*dt/2;

        //Increment r'(i) -> r'(i+1/2)
        vshiftx1 = vx1[i] + ax1*dt/2;
        vshifty1 = vy1[i] + ay1*dt/2;
        vshiftx2 = vx2[i] + ax2*dt/2;
        vshifty2 = vy2[i] + ay2*dt/2;

        //Increment r(i) -> r(i+1)
        x1[i+1] = x1[i] + vshiftx1*dt;
        y1[i+1] = y1[i] + vshifty1*dt;
        x2[i+1] = x2[i] + vshiftx2*dt;
        y2[i+1] = y2[i] + vshifty2*dt;
    }

    //Get 4 (r'')'s
    ax1 = -(bigg * m2 * (x1[i]-x2[i]))/denom;
    ay1 = -(bigg * m2 * (y1[i]-y2[i]))/denom;
    ax2 = (bigg * m1 * (x1[i]-x2[i]))/denom;
    ay2 = (bigg * m1 * (y1[i]-y2[i]))/denom;

    //Increment r'(i-1/2) -> r'(i)
    vx1[i]= vshiftx1 + ax1*dt/2;
    vy1[i]= vshifty1 + ay1*dt/2;
    vx2[i]= vshiftx2 + ax2*dt/2;
    vy2[i]= vshifty2 + ay2*dt/2;
}

//Note that the parameter t here is theta
double integrand(double u, double r, double t, double p, float gamma){
    //Calculate Common Value
    double eps = u*u + r*r - 2*u*r*cos(t);
    //Calculate Integrand
    return r*(eps +2)/(sqrt(eps*(4 + eps)))*(1-gamma*(1-(3./2)*sqrt(1-(r*r/(p*p)))));
}

//Trapezoidal double integral...
void big_integral(float gamma, double pstar[2048], double u[2048], double init_curve[2048]){

    //Init Variables needed for later
    double dr;
    double rarr[128];
    int nsteps = 128;

    //Get stepsize and grid for theta
    double dt = 2*pi/(nsteps-1);
    double tarr[128];
    for (int x=0;x<nsteps;x++){
        tarr[x]=x*dt;
    }

    //Loop over time
    for (int i=0;i<2048;i++){

        //Get stepsize and grid for r, since Pstar is a function of time...
        dr = pstar[i]/(nsteps-1);
        for (int x=0;x<nsteps;x++){
            rarr[x]=x*dr;
        }

        //Do the corners and 2 edges to avoid if's inside the big loop
        for (int j=0;j<nsteps;j++){
            if (j==0 || j==nsteps-1){
                init_curve[i] += dr*dt*integrand(u[i],rarr[j],tarr[0],pstar[i],gamma)/4;
                init_curve[i] += dr*dt*integrand(u[i],rarr[j],tarr[nsteps-1],pstar[i],gamma)/4;
            } 
            else{
                init_curve[i] += dr*dt*integrand(u[i],rarr[j],tarr[0],pstar[i],gamma)/2;
                init_curve[i] += dr*dt*integrand(u[i],rarr[j],tarr[nsteps-1],pstar[i],gamma)/2;
            }
            
        }
        //Do the other 2 sides, not including the corners
        for (int k=1;k<nsteps-1;k++){
            init_curve[i] += dr*dt*integrand(u[i],rarr[0],tarr[k],pstar[i],gamma)/2;
            init_curve[i] += dr*dt*integrand(u[i],rarr[nsteps-1],tarr[k],pstar[i],gamma)/2;
        }

        //Do the insides
        for (int j=1;j<nsteps-1;j++){
            for (int k=1;k<nsteps;k++){
                init_curve[i] += dr*dt*integrand(u[i],rarr[j],tarr[k],pstar[i],gamma);
            }
        }

        //The denominator of eqn 9 actually reduces pretty nicely...
        init_curve[i] /= (pi*pstar[i]*pstar[i]);
    }

}

//Trying it with quadrature now...
void quadquad(double pstar[2048], double u[2048], double init_curve[2048], float gamma){

    //Gaussian quad with n=21, highest order I saw that seems common
    //There are algorithms that calculate this, but can just use a table if just using one static n
    //https://pomax.github.io/bezierinfo/legendre-gauss.html
    int nsteps=21;
    double xi[21] = {0,
                    -0.1455618541608951,
                    0.1455618541608951,
                    -0.2880213168024011,
                    0.2880213168024011,
                    -0.4243421202074388,
                    0.4243421202074388,
                    -0.5516188358872198,
                    0.5516188358872198,
                    -0.6671388041974123,
                    0.6671388041974123,
                    -0.7684399634756779,
                    0.7684399634756779,
                    -0.8533633645833173,
                    0.8533633645833173,
                    -0.9200993341504008,
                    0.9200993341504008,
                    -0.9672268385663063,
                    0.9672268385663063,
                    -0.9937521706203895,
                    0.9937521706203895};

    double wi[21] = {0.1460811336496904,
                    0.1445244039899700,
                    0.1445244039899700,
                    0.1398873947910731,
                    0.1398873947910731,
                    0.1322689386333375,
                    0.1322689386333375,
                    0.1218314160537285,
                    0.1218314160537285,
                    0.1087972991671484,
                    0.1087972991671484,
                    0.0934444234560339,
                    0.0934444234560339,
                    0.0761001136283793,
                    0.0761001136283793,
                    0.0571344254268572,
                    0.0571344254268572,
                    0.0369537897708525,
                    0.0369537897708525,
                    0.0160172282577743,
                    0.0160172282577743};

    //Change of variables to get to [-1,1]:
    //theta -> pi*(eps + 1)
    //r ->  Pstar/2*(eps +1)

    //Loop over time
    for (int i=0;i<2048;i++){

        for (int j=0;j<nsteps;j++){
            for (int k=0;k<nsteps;k++){
                init_curve[i] += wi[j]*wi[k]*integrand(u[i], (pstar[i]/2)*(xi[j]+1), pi*(xi[j]+1), pstar[i], gamma);
            }
        }

        //Divide by denominator
        init_curve[i] /= (pi*pstar[i]*pstar[i]);
    }

}

//Documentation for func
int gen_light_curve(double times[2048], double light_curve[2048], float m1, float m2,
                    float rstar, float ecc, float phi, float omega, float alpha,
                    float gamma, float g, float Periods, float ti){
    
    //Since this is going to be run in a loop, refresh the values in the buffer to 0 before we start
    for (int i=0;i<2048;i++){
        light_curve[i]=0;
    }

    //Calculate some values
    double mu = (m1* m2 * ms) / (m1 + m2); //Reduced Mass
    double P = Periods * sday;  // orbital period in seconds
    double a = pow(((P*P* bigg * ms*(m1+m2)) / (4 * pi*pi)),1/3);  // SMA
    double aph = a*(1+ecc); //Aphelion
    double dt = times[1]-times[0]; //Time step for orbits
    double beta = 0.15 * (1+g) * ((15.+gamma)/(3.-gamma)); //Factor for Ellipsoidal

    //Initial Orbital Positions and Velocities
    double x1_0 = (m2 * aph) / (m1 + m2);
    double y1_0 = 0;
    double x2_0 = -(m1 * aph) / (m1 + m2);
    double y2_0 = 0;

    double vx1_0 = 0;
    double vy1_0 = sqrt((bigg * mu / a) * (m2 / m1) * ((1 - ecc) / (1 + ecc)));
    double vx2_0 = 0;
    double vy2_0 = (-m1 / m2)*sqrt((bigg * mu / a) * (m2 / m1) * ((1 - ecc) / (1 + ecc)));

    //Initialize position and velocity arrays for orbits
    double x1[2048],y1[2048],vx1[2048],vy1[2048],
            x2[2048],y2[2048],vx2[2048],vy2[2048];

    //Apply time shift
    shift_time(&x1_0,  &y1_0,  &vx1_0, &vy1_0,
                &x2_0,  &y2_0,  &vx2_0, &vy2_0, ti, m1, m2);

    //Fill Orbital Arrays
    orbit(x1_0, y1_0, vx1_0,vy1_0,x2_0, y2_0, vx2_0, vy2_0,
            x1, y1, vx1, vy1,x2, y2, vx2, vy2, dt, m1, m2);

    //Get Values at pulse peak
    int j = 0;
    double val = times[0];

    //Find where times=t0
    for (int i=1;i<2048;i++){    
        if((times[i] - t0)< val){
            j = i;
            val = times[i];
        }
    }
    double xp1 = x1[j];
    double yp1 = y1[j];
    double xp2 = x2[j];
    double yp2 = y2[j];

    //Calculate some more values...
    //Get single point values at time center
    double ap0 = sqrt((xp2 - xp1)*(xp2 - xp1) + (yp2 - yp1) *(yp2 - yp1));
    double Rs = (2 * bigg * m1 * ms) / (c*c);       //Schwarzschild Radius
    double Rein0 = sqrt(2 * Rs * ap0);              //Einstein Radius
    double u0 = (ap0 / Rein0) * sin(phi*dtor);      //Closest Angular Approach

    //Get values that are functions of time
    double ap[2048],vtran[2048], Rein[2048], te[2048],u[2048],pstar[2048];
    for(int i=0;i<2048;i++){
        ap[i] = sqrt((x2[i] - x1[i])*(x2[i] - x1[i]) + (y2[i] - y1[i]) *(y2[i] - y1[i]));
        vtran[i] = sqrt(vx1[i]*vx1[i] + vy1[i]*vy1[i]) + 
                    sqrt(vx2[i]*vx2[i] + vy2[i]*vy2[i]);
        Rein[i] = sqrt(2 * Rs * ap[i]);           
        te[i] = (Rein[i] / vtran[i]) / sday;                        //Einstein Time
        u[i] = sqrt(u0*u0 + pow((((times[i]/sday)-t0)/te[i]),2));  //Angular Separation
        pstar[i] = ((rstar * (6.957e8))/Rein[i]);
    }

    //Get the more complex values
    double sev[2048], delf[2048];
    for(int i=0;i<2048;i++){
        sev[i] = beta * ((m1/m2)*
        pow(rstar*(6.957e8)/(sqrt((xp2-xp1)*(xp2-xp1)+(yp2-yp1)*(yp2-yp1))),3)*
        cos(2*((2*pi/P)*(times[i]+(ti*sday))+pi+omega))*pow(sin((pi/2-phi*(dtor))),2));

        delf[i] = ((3-(alpha))*((((1/c)*(-1*(vx2[i]*sin(omega)*sin((pi/2-phi*dtor))+vy2[i]*cos(omega)*sin((pi/2-phi*dtor)))))*sin((pi/2-phi*dtor)))));
    }

    //Do the big numeric integral
    quadquad(pstar,u,light_curve,gamma);

    //Apply ellipsoidal and doppler factors
    for (int i=0;i<2048;i++){
        light_curve[i] *= (1 + sev[i] + delf[i]);
    }

    return 1;
}