functions {

    real df_hernquist(real[] y, real x, real yv, real zv, real rho0, real a){
        //y[1]=r, y[2]=xv
        //don't really understand the syntax real[] y

        //transform into 3d coordinates
        real pos = sqrt(square(y[1])+square(x))
        real vel = sqrt(square(y[2])+square(yv)+square(zv))

        //define parameters for the distribution function
        real M = 2.*pi()*rho0*pow(a,3.);
        real pot = -G*M/(pos+a) //FIX: need to define G
        real epsilon = -square(vel)/2. - pot;
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
    vector[N] rv_obs; //radial velocity, or the x component of velocity
    vector[N] r_obs; //this is the coordinates perpendicular to the line of sight, r^2=y^2+z^2

    //uncertainties
    vector[N] rv_err

}

transformed data {
    //combine data into y[]
    real y[N,2]; 

    for (i in 1:N) {
        y[i,1] = r_obs[i]
        y[i,2] = xv_obs[i]
    }

}

parameters {
    real<lower=0> rho0; //density scale 
    real<lower=0> a; //scale factor

    //missing position and velocity componants 
    real x; //x component of position
    real yv; //y component of velocity
    real zv; //z component of velocity

    //taken from Jeff, this is used for uncertainties 
    vector<lower=min(rv_obs - 5 * rv_err), upper=max(rv_obs + 5 * rv_err)>[N] xv;
    //also don't have any error on the radii

}

transformed parameters {
}

model {
    //hyperpriors (what are these)
    rho0 ~ //normal? gamma? uniform?
    a ~ 

    //priors (what are these?)
    rv_obs ~ normal(xv, rv_err)

    //likelihood
    for (i in 1:N) {
        y[i] ~ df_hernquist(rho0, a);
    }

}