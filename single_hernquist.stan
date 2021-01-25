//Units: distance in kpc, velocity in km/s

functions {

    real df_hernquist_lpdf(real[] y, real x, real yv, real zv, real rho0, real a){ //does this need to be divided by total mass??
        //y[1]=r, y[2]=rv

        //transform into 3d coordinates
        real pos = sqrt(square(y[1])+square(x)); //values cluster around max allowed value! (whole function not running)
        real vel = sqrt(square(y[2])+square(yv)+square(zv)); //seem to do the same here (whole function not running)

        //define parameters for the distribution function
        real M = 2.*pi()*rho0*pow(a,3.); //get the right values but then there's a line of inf's (not whole function )
        real pot = -4.3009173e-6*M/(pos+a) ; //seems to be 1e2 larger than it should be, line of zeros at end (not whole function)
        real epsilon = -square(vel)/2. - pot;
        real epsilon_tilde = epsilon*a / (4.3009173e-6*M); 
        //all epsilon_tilde values are 0.999___ (not whole function)

        //distribution function terms //CHECK UNITS
        real mult_one = pow((4.3009173e-6*M*a), -1.5) / (sqrt2()*pow((2.*pi()),3.));
        real mult_two = sqrt(epsilon_tilde) / (square(1.-epsilon_tilde)); 
        real bracket_one = (1.-2.*epsilon_tilde)* (8.*square(epsilon_tilde)-8.*epsilon_tilde-3); 
        real bracket_two = 3*asin(sqrt(epsilon_tilde)) / sqrt(epsilon_tilde*(1-epsilon_tilde));

        //put it all together!
        real f_eps = mult_one * mult_two * (bracket_one+bracket_two);

        return f_eps;

        //return ;
    }

}

data {
    int<lower=1> N; //total length of data

    vector[N] rv_obs; //radial velocity, or the x component of velocity
    vector[N] r_obs; //this is the coordinates perpendicular to the line of sight, r^2=y^2+z^2

    //uncertainties
    vector[N] rv_err;

}

transformed data {
    //combine data into y[]
    real y[N,2]; 

    for (i in 1:N) {
        y[i,1] = r_obs[i];
        y[i,2] = rv_obs[i];
    }

}

parameters {
    real<lower=0> rho0; //density scale 
    real<lower=0> a; //scale factor

    //missing position and velocity componants 
    vector<lower=-20, upper=20>[N] x; //x component of position
    vector<lower=-1000., upper=1000.>[N] yv; //y component of velocity
    vector <lower=-1000., upper=1000.>[N] zv; //z component of velocity

    //this is used for uncertainties on velocity 
    //vector<lower=min(rv_obs - 5 * rv_err), upper=max(rv_obs + 5 * rv_err)>[N] xv; 

}

transformed parameters {
}

model {
    //Wasserman et al 2018 used uniform priors over the log of the parameters
    //unclear if these are the log uniform or just the uniform
    rho0 ~ uniform(0., pow(10.,8.)); //M_sun/kpc^3
    print("rho0");
    a ~ uniform(pow(10.,-2.), pow(10., 7.)); //in kpc
    //print("a");

    //think some more about these priors, not at all sure if they're accurate
    x ~ normal(0, 40); 
    //print("x");
    yv ~ normal(0, 800); 
    //print("yv");
    zv ~ normal(0, 800);
    //print("zv");

    //rv_obs ~ normal(xv, rv_err);

    //likelihood
    for (i in 1:N) {
        print(df_hernquist_lpdf(y[i] | x[i], yv[i], zv[i], rho0, a));
        y[i] ~ df_hernquist(x[i], yv[i], zv[i], rho0, a);
    }
//note: do we want the likelihood, or the total mass?? will we calculate the total mass after getting parameter values?

}