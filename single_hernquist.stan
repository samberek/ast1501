functions {

    
}

data {
    int<lower=1> N; //total length of data

    vector[N] ra; //I don't really know if I need these two
    vector[N] dec; 

    // measured data
    vector[N] rv_obs; 
    vector[N] radius_obs; 

    //uncertainties
    vector[N] rv_err_obs

}

transformed data {

}

parameters {

}

transformed parameters {


}

model {


}