#include <smile.h>
//Does the initial time shift for orbits
//Takes a set of ri, r'i and increments them through time up to ti
void shift_time(double *x1_0, double *y1_0, double *vx1_0,double *vy1_0,
            double *x2_0, double *y2_0, double *vx2_0,double *vy2_0,
            float ti, float m1, float m2, double rs){
    
    double dt0 = ti*sday/1000;

    //Calculate common denominator for r''[0]
    double denom = pow( sqrt((*x1_0-*x2_0)*(*x1_0-*x2_0) + (*y1_0-*y2_0)*(*y1_0-*y2_0)) - rs ,2) * sqrt((*x1_0-*x2_0)*(*x1_0-*x2_0) + (*y1_0-*y2_0)*(*y1_0-*y2_0));

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
        denom = pow( sqrt((*x1_0-*x2_0)*(*x1_0-*x2_0) + (*y1_0-*y2_0)*(*y1_0-*y2_0)) - rs ,2) * sqrt((*x1_0-*x2_0)*(*x1_0-*x2_0) + (*y1_0-*y2_0)*(*y1_0-*y2_0));
        
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

    //Calculate common denominator for r''[i]
    denom = pow( sqrt((*x1_0-*x2_0)*(*x1_0-*x2_0) + (*y1_0-*y2_0)*(*y1_0-*y2_0)) - rs ,2) * sqrt((*x1_0-*x2_0)*(*x1_0-*x2_0) + (*y1_0-*y2_0)*(*y1_0-*y2_0));

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
            double x2[2048], double y2[2048], double vx2[2048],double vy2[2048], double dt, float m1, float m2, double rs){
    
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
    double denom = pow( sqrt((x1_0-x2_0)*(x1_0-x2_0) + (y1_0-y2_0)*(y1_0-y2_0)) - rs ,2) * sqrt((x1_0-x2_0)*(x1_0-x2_0) + (y1_0-y2_0)*(y1_0-y2_0));

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
        denom = pow( sqrt((x1[i]-x2[i])*(x1[i]-x2[i]) + (y1[i]-y2[i])*(y1[i]-y2[i])) - rs ,2) * sqrt((x1[i]-x2[i])*(x1[i]-x2[i]) + (y1[i]-y2[i])*(y1[i]-y2[i]));
        
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

    //Calculate common denominator for r''[i]
    denom = pow( sqrt((x1[i]-x2[i])*(x1[i]-x2[i]) + (y1[i]-y2[i])*(y1[i]-y2[i])) - rs ,2) * sqrt((x1[i]-x2[i])*(x1[i]-x2[i]) + (y1[i]-y2[i])*(y1[i]-y2[i]));

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
int big_integral(double pstar[2048], double u[2048], double init_curve[2048],float gamma){

    //Init Variables needed for later
    double dr;
    double rarr[512];
    int nsteps = 512;

    //Get stepsize and grid for theta
    double dt = 2*pi/(nsteps-1);
    double tarr[512];
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
    return 1;
}

//Trying it with quadrature now...
int quadquad(double pstar[2048], double u[2048], double init_curve[2048], float gamma){

    //Gaussian quad with n=64, highest order I found a table for
    //There are algorithms that calculate this, but can just use a table if just using one static n
    //https://pomax.github.io/bezierinfo/legendre-gauss.html
    int nsteps=64;
    double xi[64] = {-0.0243502926634244,0.0243502926634244,
                    -0.0729931217877990,0.0729931217877990,
                    -0.1214628192961206,0.1214628192961206,
                    -0.1696444204239928,0.1696444204239928,
                    -0.2174236437400071,0.2174236437400071,
                    -0.2646871622087674,0.2646871622087674,
                    -0.3113228719902110,0.3113228719902110,
                    -0.3572201583376681,0.3572201583376681,
                    -0.4022701579639916,0.4022701579639916,
                    -0.4463660172534641,0.4463660172534641,
                    -0.4894031457070530,0.4894031457070530,
                    -0.5312794640198946,0.5312794640198946,
                    -0.5718956462026340,0.5718956462026340,
                    -0.6111553551723933,0.6111553551723933,
                    -0.6489654712546573,0.6489654712546573,
                    -0.6852363130542333,0.6852363130542333,
                    -0.7198818501716109,0.7198818501716109,
                    -0.7528199072605319,0.7528199072605319,
                    -0.7839723589433414,0.7839723589433414,
                    -0.8132653151227975,0.8132653151227975,
                    -0.8406292962525803,0.8406292962525803,
                    -0.8659993981540928,0.8659993981540928,
                    -0.8893154459951141,0.8893154459951141,
                    -0.9105221370785028,0.9105221370785028,
                    -0.9295691721319396,0.9295691721319396,
                    -0.9464113748584028,0.9464113748584028,
                    -0.9610087996520538,0.9610087996520538,
                    -0.9733268277899110,0.9733268277899110,
                    -0.9833362538846260,0.9833362538846260,
                    -0.9910133714767443,0.9910133714767443,
                    -0.9963401167719553,0.9963401167719553,
                    -0.9993050417357722,0.9993050417357722};

    double wi[64] = {0.048690957009139,0.048690957009139,
                    0.0485754674415034,0.0485754674415034,
                    0.0483447622348030,0.0483447622348030,
                    0.0479993885964583,0.0479993885964583,
                    0.0475401657148303,0.0475401657148303,
                    0.0469681828162100,0.0469681828162100,
                    0.0462847965813144,0.0462847965813144,
                    0.0454916279274181,0.0454916279274181,
                    0.0445905581637566,0.0445905581637566,
                    0.0435837245293235,0.0435837245293235,
                    0.0424735151236536,0.0424735151236536,
                    0.0412625632426235,0.0412625632426235,
                    0.0399537411327203,0.0399537411327203,
                    0.0385501531786156,0.0385501531786156,
                    0.0370551285402400,0.0370551285402400,
                    0.0354722132568824,0.0354722132568824,
                    0.0338051618371416,0.0338051618371416,
                    0.0320579283548516,0.0320579283548516,
                    0.0302346570724025,0.0302346570724025,
                    0.0283396726142595,0.0283396726142595,
                    0.0263774697150547,0.0263774697150547,
                    0.0243527025687109,0.0243527025687109,
                    0.0222701738083833,0.0222701738083833,
                    0.0201348231535302,0.0201348231535302,
                    0.0179517157756973,0.0179517157756973,
                    0.0157260304760247,0.0157260304760247,
                    0.0134630478967186,0.0134630478967186,
                    0.0111681394601311,0.0111681394601311,
                    0.0088467598263639,0.0088467598263639,
                    0.0065044579689784,0.0065044579689784,
                    0.0041470332605625,0.0041470332605625,
                    0.0017832807216964,0.0017832807216964};

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
    return 1;
}

int piecequad(double pstar[2048], double u[2048], double init_curve[2048], float gamma){
    
    int nsteps=64;

    double xi[64] = {-0.0243502926634244,0.0243502926634244,
                    -0.0729931217877990,0.0729931217877990,
                    -0.1214628192961206,0.1214628192961206,
                    -0.1696444204239928,0.1696444204239928,
                    -0.2174236437400071,0.2174236437400071,
                    -0.2646871622087674,0.2646871622087674,
                    -0.3113228719902110,0.3113228719902110,
                    -0.3572201583376681,0.3572201583376681,
                    -0.4022701579639916,0.4022701579639916,
                    -0.4463660172534641,0.4463660172534641,
                    -0.4894031457070530,0.4894031457070530,
                    -0.5312794640198946,0.5312794640198946,
                    -0.5718956462026340,0.5718956462026340,
                    -0.6111553551723933,0.6111553551723933,
                    -0.6489654712546573,0.6489654712546573,
                    -0.6852363130542333,0.6852363130542333,
                    -0.7198818501716109,0.7198818501716109,
                    -0.7528199072605319,0.7528199072605319,
                    -0.7839723589433414,0.7839723589433414,
                    -0.8132653151227975,0.8132653151227975,
                    -0.8406292962525803,0.8406292962525803,
                    -0.8659993981540928,0.8659993981540928,
                    -0.8893154459951141,0.8893154459951141,
                    -0.9105221370785028,0.9105221370785028,
                    -0.9295691721319396,0.9295691721319396,
                    -0.9464113748584028,0.9464113748584028,
                    -0.9610087996520538,0.9610087996520538,
                    -0.9733268277899110,0.9733268277899110,
                    -0.9833362538846260,0.9833362538846260,
                    -0.9910133714767443,0.9910133714767443,
                    -0.9963401167719553,0.9963401167719553,
                    -0.9993050417357722,0.9993050417357722};

    double wi[64] = {0.048690957009139,0.048690957009139,
                    0.0485754674415034,0.0485754674415034,
                    0.0483447622348030,0.0483447622348030,
                    0.0479993885964583,0.0479993885964583,
                    0.0475401657148303,0.0475401657148303,
                    0.0469681828162100,0.0469681828162100,
                    0.0462847965813144,0.0462847965813144,
                    0.0454916279274181,0.0454916279274181,
                    0.0445905581637566,0.0445905581637566,
                    0.0435837245293235,0.0435837245293235,
                    0.0424735151236536,0.0424735151236536,
                    0.0412625632426235,0.0412625632426235,
                    0.0399537411327203,0.0399537411327203,
                    0.0385501531786156,0.0385501531786156,
                    0.0370551285402400,0.0370551285402400,
                    0.0354722132568824,0.0354722132568824,
                    0.0338051618371416,0.0338051618371416,
                    0.0320579283548516,0.0320579283548516,
                    0.0302346570724025,0.0302346570724025,
                    0.0283396726142595,0.0283396726142595,
                    0.0263774697150547,0.0263774697150547,
                    0.0243527025687109,0.0243527025687109,
                    0.0222701738083833,0.0222701738083833,
                    0.0201348231535302,0.0201348231535302,
                    0.0179517157756973,0.0179517157756973,
                    0.0157260304760247,0.0157260304760247,
                    0.0134630478967186,0.0134630478967186,
                    0.0111681394601311,0.0111681394601311,
                    0.0088467598263639,0.0088467598263639,
                    0.0065044579689784,0.0065044579689784,
                    0.0041470332605625,0.0041470332605625,
                    0.0017832807216964,0.0017832807216964};
    
    //Loop over time
    for (int i=0;i<2048;i++){

        //Case 1: u < p (Needed for peak)
        if (u[i] < pstar[i]){
            //Region 1: Pi/4 < t < 7Pi/4, 0 < r < u
            for (int j=0;j<nsteps;j++){
                for (int k=0;k<nsteps;k++){
                    init_curve[i] += wi[j]*wi[k]*integrand(u[i],(u[i]/2)*(xi[j]+1),
                    (3*pi/4)*xi[k]+pi,pstar[i],gamma);
                }
            }
            //Region 2: Pi/4 < t < 7Pi/4, u < r < P
            for (int j=0;j<nsteps;j++){
                for (int k=0;k<nsteps;k++){
                    init_curve[i] += wi[j]*wi[k]*integrand(u[i],(pstar[i]-u[i])/2*xi[j]+
                    (pstar[i]+u[i])/2, (3*pi/4)*xi[k]+pi,pstar[i],gamma);
                }
            }
        //Case 2: u >=p
        } else {
            //Region 1: Pi/4 < t < 7Pi/4, 0 < r < 9P/10
            for (int j=0;j<nsteps;j++){
                for (int k=0;k<nsteps;k++){
                    init_curve[i] += wi[j]*wi[k]*integrand(u[i],(9*pstar[i]/10)*(xi[j]+1), (3*pi/4)*xi[k]+pi,pstar[i],gamma);
                }
            }
            //Region 2: Pi/4 < t < 7Pi/4, 9P/10 < r < P
            for (int j=0;j<nsteps;j++){
                for (int k=0;k<nsteps;k++){
                    init_curve[i] += wi[j]*wi[k]*integrand(u[i],(pstar[i]/20)*xi[j]+
                    19*pstar[i]/20, (3*pi/4)*xi[k]+pi,pstar[i],gamma);
                }
            }
            //Region 1: 7*Pi/4 < t < 9Pi/4, 0 < r < 9P/10
            for (int j=0;j<nsteps;j++){
                for (int k=0;k<nsteps;k++){
                    init_curve[i] += wi[j]*wi[k]*integrand(u[i],(9*pstar[i]/10)*(xi[j]+1), (pi/4)*xi[k]+4*pi,pstar[i],gamma);
                }
            }
            //Region 2: 7*Pi/4 < t < 9Pi/4, 9P/10 < r < P
            for (int j=0;j<nsteps;j++){
                for (int k=0;k<nsteps;k++){
                    init_curve[i] += wi[j]*wi[k]*integrand(u[i],(pstar[i]/20)*xi[j]+
                    19*pstar[i]/20, (pi/4)*xi[k]+4*pi,pstar[i],gamma);
                }
            }
        }
        

        //Divide by denominator
        init_curve[i] /= (pi*pstar[i]*pstar[i]);
    }
    return 1;
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
    double rs = 2 * bigg * ms * (m1+m2)/(c*c);    //Schwartzchild Radius
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
                &x2_0,  &y2_0,  &vx2_0, &vy2_0, ti, m1, m2, rs);

    //Fill Orbital Arrays
    orbit(x1_0, y1_0, vx1_0,vy1_0,x2_0, y2_0, vx2_0, vy2_0,
            x1, y1, vx1, vy1,x2, y2, vx2, vy2, dt, m1, m2, rs);

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
    a = quadquad(pstar,u,light_curve,gamma);

    //Apply ellipsoidal and doppler factors
    for (int i=0;i<2048;i++){
        light_curve[i] *= (1 + sev[i] + delf[i]);
    }

    return a;
}

//Calling functions and timing things
int main(void){

    clock_t start_t, end_t;

    //Set input parameters
    float m1 = 2.4;
    float m2 = 0.5;
    float rstar = 0.95;
    float ecc = 0.1;
    float phi = 0.001;
    float omega = 1;
    float alpha = -1.3;
    float gamma = 0.4;
    float g = 0.4;
    float Periods = 14.1;
    float ti = 12.5;

    //Setup arrays
    double times[2048];
    double light_curve[2048];

    for(int i=0;i>2048;i++){
        times[i]= i * (14 * sday/2048-1);
    }
    
    
    start_t=clock();
    int a = gen_light_curve(times, light_curve, m1, m2, rstar,ecc, phi,omega,
                            alpha, gamma, g, Periods, ti);
    end_t = clock();

    printf("Time for Eval %f\n",(double)(end_t - start_t)/CLOCKS_PER_SEC);
    printf("Return code %d\n",a);
    return 0;
}
