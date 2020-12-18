functions {

    real df_hernquist(real[] y, real rho0, real a){
        //y[1]=r, y[2]=rv
        //don't really understand the syntax real[] y, ask Jeff

        //define parameters for the distribution function
        real M = 2.*pi()*rho0*pow(a,3.);
        real pot = -G*M/(y[1]+a) //FIX: need to define G
        real epsilon = -square(y[2])/2. - pot;
        real epsilon_tilde = epsilon*a / (G*M);

        //distribution function terms
        real mult_one = pow((G*M*a), -1.5) / (sqrt(2.)*pow((2.*pi()),3.));
        real mult_two = sqrt(epsilon_tilde) / (square(1.-epsilon_tilde));
        real bracket_one = (1.-2.*epsilon_tilde) * (8.*square(epsilon_tilde)-8.*epsilon_tilde-3);
        real bracket_two = 3*asin(sqrt(epsilon_tilde)) / sqrt(epsilon_tilde*(1-epsilon_tilde));

        //put it all together!
        real f_eps = mult_one * mult_two * (bracket_one+bracket_two);

        return f_eps;

    }

}

data {
    int<lower=1> N; //total length of data

    //I don't really know if I need these two
    vector[N] ra_obs; 
    vector[N] dec_obs; 

    // definitely need these
    vector[N] rv_obs; 
    vector[N] r_obs; 

    //uncertainties
    vector[N] rv_err

}

transformed data {
    //combine data into y[]
    real y[N,3]; 

    for (i in 1:N) {
        y[i,1] = r_obs[i]
        y[i,2] = rv_obs[i]
    }

}

parameters {
    real<lower=0> rho0; //density scale 
    real<lower=0> a; //scale factor

    //taken from Jeff, this is used for uncertainties but not sure if I need
    vector<lower=min(rv_obs - 5 * rv_err), upper=max(rv_obs + 5 * rv_err)>[N] rv;
    //also don't have any error on the radii

}

transformed parameters {
}

model {
    //hyperpriors (what are these)
    rho0 ~ //normal? gamma? uniform?
    a ~ 

    //priors (what are these?)
    rv_obs ~ normal(rv, rv_err)

    //likelihood
    for (i in 1:N) {
        y[i] ~ df_hernquist(rho0, a);
    }

}