#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <chrono>

using namespace std;

int n;   // n Zeilen 
int stick_length_mm=1000;   // Länge Stab in mm
double omega;

int do_step(Eigen::VectorXd & A) {
    Eigen::VectorXd A_cpy(n);
    A_cpy=A;
    for(int i = 1; i<(n-1); i++) {
        A_cpy(i) = (1-omega)*A(i)+0.5*omega*(A(i-1)+(A(i+1)));
    }
    A = A_cpy;
    return 0;
}

int main(int argc, char *argv[]) {
    double step_size=0.01;    
    double step_size_mm=stick_length_mm*step_size;
    double time_step_size=0.00004; // in Sekunden 

    if(2*time_step_size/(step_size*step_size) >= 1) {
        cout<<"Choose different time/space step length. Will not converge.\n";
        cout<<"time_step_size has to be < "<<(step_size*step_size)/2<<endl;
        return -1;
    }

    n=(int)((stick_length_mm)/step_size_mm);
    omega=2*time_step_size/(step_size*step_size);
    cout<<omega<<endl;
    Eigen::VectorXd x(n);

    x(0)=0.0;
    for(int i = 1; i<n-1; i++) {
        x(i) = sin(M_PI*((double)i/(double)n));
    }
    x(n-1)=0.0;

    for(double i = time_step_size; i<=1.00001; i+=time_step_size) {
        do_step(x);
        //cout<<"Vec:\n"<<x<<endl;
    }

    cout<<"Vec:\n"<<x<<endl;
    Eigen::VectorXd exact(n); 

    exact(0)=0;
    for(int i = 1; i<n-1; i++) {
        exact(i)=exp(-M_PI*M_PI)*sin(M_PI*((double)i/(double)n));
    }
    exact(n-1)=0;   
    cout<<"Exact: \n"<<exact<<endl;


    double sum=0;
    double cnt=0;
    for(int i = 1; i<(n-1); i++) {
        sum+=abs(exact(i)-x(i));
        cnt++;
    }
    cout<<"Mean of local error: "<<sum/cnt<<endl;

    return 0;
}