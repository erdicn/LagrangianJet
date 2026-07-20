#ifndef DROPLET_STATS_HPP
#define DROPLET_STATS_HPP

#include <vector>
#include <random>

static std::vector<double> d11mm   = {13.09325314   ,  14.24700768  ,  14.67278033 ,  15.21765021 ,  15.34865718 ,  16.85672602 ,  19.29970821 , 20.46537248   ,  22.35157506 ,  23.63038171 ,  27.04996129 ,  28.66372893 ,  31.39403323 ,  35.18876913  , 39.51051033  ,  44.91752516 ,  49.79753469 ,  54.39617698 ,  60.01607813 ,  65.65979873  ,  70.97153576 , 77.18841184  ,  83.53629488 ,  91.95646996 ,  97.98874531 ,  102.3953433  , 104.19966653 , 107.47781814 , 109.23152504 , 109.82998869 ,  111.28595248, 112.22086584  , 112.91907342 , 114.19192521 , 116.22699934 , 118.62680879 , 122.37837194  , 124.42088966 , 128.41809087 , 131.08884654 , 131.95081284 , 133.5854225  , 135.42100875  ,  138.08283213, 141.55005062 , 145.59340201 , 146.89007324 , 149.01149288  , 151.51104627 , 153.86768296 , 155.13457988 , 158.35169416 , 160.93014947  , 164.2678497  , 169.27588876 , 174.08592866 , 183.96951111 , 199.09039481 , 203.2677306   , 207.82915501 ,212.72851783  };
static std::vector<double> pdf11mm = {1.35905823e-06, 3.56256379e-07,3.19752661e-06,8.82414627e-06,2.44000409e-05,3.05447903e-05,1.62310523e-04, 4.74616453e-04,9.65072753e-04,1.56940845e-03,2.50936352e-03,4.28042777e-03,5.23820066e-03, 4.51735104e-03,3.40251084e-03,1.05183956e-03,3.74332614e-04,1.69465537e-04,6.12330808e-05, 2.49485275e-05,1.28064950e-05,7.15621394e-06,4.76257228e-06,6.39886990e-06,1.93992694e-05, 5.61269255e-05,1.48913944e-04,1.75647324e-03,7.96209420e-03,1.07174765e-02,1.71425059e-02, 2.20386024e-02,3.28771256e-02,4.02628853e-02,4.70227639e-02,4.83260523e-02, 4.50827041e-02,3.29853411e-02,1.72098620e-02,1.04599941e-02,8.12572287e-03,6.89797509e-03, 8.20371195e-03,8.73166990e-03,7.67036781e-03,6.47695069e-03,5.21111440e-03, 3.79113647e-03,3.08170932e-03,3.05279567e-03,3.46640142e-03,3.53208434e-03, 2.78714693e-03,2.06564951e-03,1.51502535e-03,8.73645246e-04,1.99036272e-04,1.24280302e-05, 4.60972073e-06,1.48597155e-06,4.37988533e-07};
static std::random_device rd; // seed
static std::mt19937 gen(rd()); // TODO convert to 64 and also seed it from a set of random devices
static std::piecewise_linear_distribution<double> distribution11mm(d11mm.begin(), d11mm.end(), pdf11mm.begin());
static std::normal_distribution normal_dist11mm{10.15, 0.22};
// static std::linear_congruential_engine<uint32_t, 65539, 0, 2147483648> gen(rd());
// static std::piecewise_linear_distribution<double> distribution11mm_satelite(d11mm.begin(), d11mm.end(), pdf11mm.begin());
// static std::normal_distribution normal_distfreq{91000., 2*942.6869706212264};
// static std::uniform_real_distribution<double> uniform_dist_vel_9_11(9.0, 11.0);
// static std::uniform_real_distribution<double> uniform_dist_dia_100_165um(100, 165);
// std::normal_distribution normal_distfreq{43492.14887184294, 942.6869706212264};

static std::normal_distribution<double> normal_dist_droplet_enlet_time{29.25, 1};
// normal_dist_droplet_enlet_time.param(std::normal_distribution<double>::param_type(10.0, 25.0)); // to redo the distribution differently
static void updateNormalDropletDistributionInletTime(double new_mean, double new_stddev) {
    normal_dist_droplet_enlet_time.param(
        std::normal_distribution<double>::param_type(new_mean, new_stddev)
    );
}

static const double sat_mu = 10.155; // Mean velocity (center of the peak)
static const double sat_b  = 0.58;   // Scale parameter (adjust this to make the wings wider/narrower)

template <typename Generator>
double getLaplaceRandom(double mu, double b, Generator& gener) {
    std::exponential_distribution<double> exp_dist(1.0 / b);
    std::bernoulli_distribution coin_flip(0.5); // 50% chance for + or -
    
    // Draw an exponential value
    double exp_val = exp_dist(gener);
    
    // Randomly add or subtract it from the mean
    if (coin_flip(gener)) {
        return mu + exp_val;
    } else {
        return mu - exp_val;
    }
}


#endif /* DROPLET_STATS_HPP */
