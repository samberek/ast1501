//Units: distance in kpc, velocity in km/s

functions {

    real df_hernquist_lpdf(real[] y, real x, real v_yz, real logM, real a){
        //y[1]=r, y[2]=rv

        //exponentiate parameters
        real M = pow(10.,logM);
        //real a = pow(10, loga);

        //transform into 3d coordinates
        real pos = sqrt(square(y[1])+square(x)); 
        real velsq = square(y[2])+square(v_yz); 

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
        real f_eps = mult_one * mult_two * (bracket_one+bracket_two);

        if (epsilon_tilde < 0.) {
            return 0;
        }  

        if (epsilon_tilde > 1.) {
            return 0;
        }         

        return f_eps;

    }

    vector max_position(vector r_obs, int N){
        //set the size of the galaxy
        real R = 120;

        //y[1] is the position data, r, where r^2 = y^2+z^2

        vector[N] x_max_square = square(R) - square(r_obs);
        vector[N] x_max = sqrt(x_max_square);

        return x_max;
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
    real<lower=12, upper=13.7> logM; //total mass
    real<lower=0> a; //scale factor

    //missing position and velocity componants 
    vector<lower=0, upper=1>[N] x_raw; //x component of position
    vector<lower=-0, upper=1500.>[N] v_yz; //combined y and z component of velocity

    //this is used for uncertainties on velocity, comment for troubleshooting
    //vector<lower=min(rv_obs - 5 * rv_err), upper=max(rv_obs + 5 * rv_err)>[N] xv; 

}

transformed parameters {
    vector[N] x = max_position(r_obs, N) .* x_raw;
}

model {
    //Wasserman et al 2018 used uniform priors over the log of the parameters
    logM ~ normal(13.2,0.1); 
    a ~ normal(5, 3);

    x ~ normal(0, 0.25 * max_position(r_obs, N)); 
    v_yz ~ normal(0, 300); 

    //rv_obs ~ normal(xv, rv_err); //adding uncertainty, comment for troubleshooting 

    //likelihood
    for (i in 1:N) {
        //print(df_hernquist_lpdf(y[i] | x[i], yv[i], zv[i], rho0, a));
        y[i] ~ df_hernquist(x[i], v_yz[i], logM, a);
    }

}