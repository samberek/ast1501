//Units: distance in kpc, velocity in km/s

functions {

    real df_hernquist_lpdf(real[] y, real logM, real loga){
        //y[1]=r, y[2]=rv

        //exponentiate parameters
        real M = pow(10.,logM);
        real a = pow(10, loga);

        //define parameters for the distribution function
        //real pot = -(4.3009173e-6*M)/(y[1]+a) ; 
        real epsilon = -square(y[2])/2. + (4.3009173e-6*M)/(y[1]+a);
        real epsilon_tilde = epsilon*a / (4.3009173e-6*M); 

        //distribution function terms 
        real mult_one = pow((4.3009173e-6*M*a), -1.5) / (sqrt2()*pow((2.*pi()),3.));
        real mult_two = sqrt(epsilon_tilde) / (square(1.-epsilon_tilde)); 
        real bracket_one = (1.-2.*epsilon_tilde) * (8.*square(epsilon_tilde)-8.*epsilon_tilde-3); 
        real bracket_two = 3*asin(sqrt(epsilon_tilde)) / sqrt(epsilon_tilde*(1-epsilon_tilde));

        //put it all together!
        real f_eps = mult_one * mult_two * (bracket_one+bracket_two);

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
    real<lower=11, upper=13.5> logM; //total mass
    real<lower=0> loga; //scale factor
}

transformed parameters {
}

model {
    //Wasserman et al 2018 used uniform priors with these values over the log of the parameters
    //unclear if these are the log uniform or just the uniform
    //rho0 ~ uniform(1e2, 1e8); //M_sun/kpc^3
    logM ~ normal(12.5,0.5); 
    loga ~ uniform(1, 3);

    //likelihood
    for (i in 1:N) {
        //print(df_hernquist_lpdf(y[i] | x[i], yv[i], zv[i], rho0, a));
        y[i] ~ df_hernquist(logM, loga);
    }
//note: do we want the likelihood, or the total mass?? will we calculate the total mass after getting parameter values?

}