#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <functional>
#include <fstream>

using namespace std;

random_device rd;
mt19937 gen(rd());

ostream& operator << (ostream& o, vector<double>v) {
    for (auto e : v) {
        o << e << ",";

    }
    return o;
}


vector<double> hill_climbing(function<double(vector<double>)> f, function<bool(vector<double>)> f_domain, ofstream &out, vector<double> p0, int iterations)
{
    auto p = p0;
    uniform_int_distribution<> distrib(0, p.size() - 1);
    uniform_real_distribution<> distrib_r(-0.1, 0.1);
    if (!f_domain(p)) throw std::invalid_argument("The p0 point must be in domain");
    for (int i = 0; i < iterations; i++) {
        auto p2 = p;
        p2[distrib(gen)] += distrib_r(gen);
        if (f_domain(p2)) {
            double y2 = f(p2);
            if (y2 < f(p)) {
                p = p2;
            }
        }
        out << i <<" "<<f(p)<<endl;

    }
    return p;
}
bool holder_table_domain (std::vector<double> v){
    if(abs(v[0]) <= 10 && abs(v[1]) <= 10) {
        return true;
    }
    return false;
}

double holder_table_function(std::vector<double> v){
    if(!holder_table_domain(v)) {
        return 9999999;
    }

    double result = 0.0;

    std::cout << "{" << v[0] << " " << v[1] << "} = ";

    result = -abs(sin(v[0])* cos(v[1])* exp(abs(1-(sqrt(v[0]*v[0]+v[1]*v[1])/M_PI))));
    std::cout << result << std::endl;

    return result;

}

bool goldstein_domain (std::vector<double> v){
    if(abs(v[0]) <= 2 && abs(v[1]) <= 2) {
        return true;
    }
    return false;
}

double goldstein_function(std::vector<double> v){
    if(!goldstein_domain(v)) {
        return 9999999;
    }

    double result = 0.0;

    std::cout << "{" << v[0] << " " << v[1] << "} = ";

    result = (1+(v[0]+v[1]+1)*(v[0]+v[1]+1)*(19-14*v[0]+3*v[0]*v[0]-14*v[1]+6*v[0]*v[1]+3*v[1]*v[1]))*
             (30+(2*v[0]-3*v[1])*(2*v[0]-3*v[1])*(18-32*v[0]+12*v[0]*v[0]+48*v[1]-36*v[0]*v[1]+27*v[1]*v[1]));

    std::cout << result << std::endl;

    return result;

}

vector<double> simulated_annealing(function<double(vector<double>)> f, function<bool(vector<double>)> f_domain, vector<double> p0, int iterations, ofstream &out,function<vector<double>(vector<double>)> N, function<double(int)> T)
{
    auto currentPoint = p0;
    auto best = p0;


    uniform_real_distribution<> u_k(0.0, 1.0);

    if (!f_domain(currentPoint)) throw std::invalid_argument("The p0 point must be in domain");
    for (int k = 0; k < iterations; k++) {
        auto nextPoint = N(currentPoint);
        if (!f_domain(nextPoint)) continue;
        if (f(nextPoint) < f(currentPoint)) {
            currentPoint = nextPoint;
        }
        else {
            double u = u_k(gen);
            if (u < exp(-abs(f(nextPoint) - f(currentPoint)) / T(k))) {
                currentPoint = nextPoint;
            }
            else {

            }
        }
        if (f(currentPoint) < f(best)) {
            best = currentPoint;
        }
        out << k <<" "<<f(currentPoint)<<endl;
    }
    return best;
}

int main() {
    ofstream outfile("holder_hill.txt");
    ofstream outfile2("holder_annealing.txt");
    ofstream outfile3("gold_hill.txt");
    ofstream outfile4("gold_annealing.txt");


    int iterations;
    cout << "Give number of iterations: " << endl;
    cin >> iterations;
    uniform_real_distribution<> distrib_h(-10, 10);
    vector<double> hp0 = {distrib_h(gen), distrib_h(gen) };

    auto result1 = hill_climbing(holder_table_function, holder_table_domain, outfile, hp0, iterations);
    auto result2 = simulated_annealing(holder_table_function, holder_table_domain, hp0, iterations, outfile2, [](auto p) {
                                           normal_distribution<double> n(0.0, 0.3);
                                           for (auto& e : p) {
                                               e = e + n(gen);
                                           }return p;},
                                       [](int k) { return 1000.0 / k; });

    uniform_real_distribution<> distrib_g(-2, 2);
    vector<double> gp0 = {distrib_g(gen), distrib_g(gen) };

    auto result3 = hill_climbing(goldstein_function, goldstein_domain,outfile3,gp0, iterations);
    auto result4 = simulated_annealing(goldstein_function, goldstein_domain, gp0, iterations, outfile4,[](auto p) {
                                           normal_distribution<double> n(0.0, 0.3);
                                           for (auto& e : p) {
                                               e = e + n(gen);
                                           }return p;},
                                       [](int k) { return 1000.0 / k; });
    return 0;
}