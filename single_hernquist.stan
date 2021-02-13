//Units: distance in kpc, velocity in km/s

functions {

    real df_hernquist_lpdf(real[] y, real x, real yv, real zv, real logM, real loga){
        //y[1]=r, y[2]=rv

        //exponentiate parameters
        real M = pow(10.,logM);
        real a = pow(10, loga);

        //transform into 3d coordinates
        real pos = sqrt(square(y[1])+square(x)); 
        real velsq = square(y[2])+square(yv)+square(zv); 

        //define parameters for the distribution function
        //real M = 2.*pi()*rho0*pow(a,3.); 
        real pot = -4.3009173e-6*M/(pos+a) ; 
        real epsilon = -velsq/2. - pot;
        real epsilon_tilde = epsilon*a / (4.3009173e-6*M); 

        //distribution function terms 
        real mult_one = pow((4.3009173e-6*M*a), -1.5) / (sqrt2()*pow((2.*pi()),3.));
        real mult_two = sqrt(epsilon_tilde) / (square(1.-epsilon_tilde)); 
        real bracket_one = (1.-2.*epsilon_tilde)* (8.*square(epsilon_tilde)-8.*epsilon_tilde-3); 
        real bracket_two = 3*asin(sqrt(epsilon_tilde)) / sqrt(epsilon_tilde*(1-epsilon_tilde));

        //put it all together!
        real f_eps = mult_one * mult_two * (bracket_one+bracket_two) / M;

        if (epsilon_tilde < 0.) {
            return 0;
        }  

        if (epsilon_tilde > 1.) {
            return 0;
        }         

        return f_eps;

    }

}

data {
    int<lower=1> N; //total length of data

    vector[N] rv_obs; //radial velocity, or the x component of velocity
    vector[N] r_obs; //coordinates perpendicular to the line of sight, r^2=y^2+z^2

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
    //real<lower=0> rho0; //density scale 
    real<lower=10, upper=16> logM; //total mass
    real<lower=0> loga; //scale factor

    //missing position and velocity componants 
    vector<lower=0, upper=100.>[N] x; //x component of position
    vector<lower=-0, upper=1000.>[N] yv; //y component of velocity
    vector <lower=-0, upper=1000.>[N] zv; //z component of velocity

    //this is used for uncertainties on velocity, comment for troubleshooting
    vector<lower=min(rv_obs - 5 * rv_err), upper=max(rv_obs + 5 * rv_err)>[N] xv; 

}

transformed parameters {
}

model {
    //Wasserman et al 2018 used uniform priors with these values over the log of the parameters
    //unclear if these are the log uniform or just the uniform
    //rho0 ~ uniform(1e2, 1e8); //M_sun/kpc^3
    logM ~ uniform(10,15);
    loga ~ uniform(1, 7);

    //think some more about these priors
    x ~ normal(0, 40); 
    yv ~ normal(0, 200); 
    zv ~ normal(0, 200);

    rv_obs ~ normal(xv, rv_err); //adding uncertainty, comment for troubleshooting 

    //likelihood
    for (i in 1:N) {
        //print(df_hernquist_lpdf(y[i] | x[i], yv[i], zv[i], rho0, a));
        y[i] ~ df_hernquist(x[i], yv[i], zv[i], logM, loga);
    }
//note: do we want the likelihood, or the total mass?? will we calculate the total mass after getting parameter values?

}