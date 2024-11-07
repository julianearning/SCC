#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <chrono>

using namespace std;

int n;   // n Zeilen 
int gitarrensaite_mm=1000;   // Länge Stab in mm
double lambda;
double omega;


int do_first_step(Eigen::VectorXd & A, Eigen::VectorXd & B) {
    Eigen::VectorXd A_cpy(n);
    A_cpy=A;
    for(int i = 1; i<(n-1); i++) {
        A_cpy(i) = (1-omega)*A(i)+0.5*omega*(A(i-1)+(A(i+1)));
    }
    B = A_cpy;
    return 0;
}

int do_step(Eigen::VectorXd & x, Eigen::VectorXd & x1) {
    Eigen::VectorXd new_x(n);
    
    double lambda_sq=lambda*lambda;
    new_x(0)=0;
    for(int i = 1; i<n-1; i++) {
        new_x(i)=2*(1-lambda_sq)*x1(i)+lambda_sq*(x1(i-1)+x1(i+1))-x(i);
    }
    new_x(n-1)=0;

    x=x1;
    x1=new_x;
    return 0;
}


int main(int argc, char *argv[]) {
    double step_size=0.1;    
    double step_size_mm=gitarrensaite_mm*step_size;
    double time_step_size=0.001; // in Sekunden

    if(2*time_step_size/(step_size*step_size) >= 1) {
        cout<<"Choose different time/space step length. Will not converge.\n";
        cout<<"time_step_size has to be < "<<(step_size*step_size)/2<<endl;
        return -1;
    }

    n=(int)((gitarrensaite_mm)/step_size_mm);
    omega=2*time_step_size/(step_size*step_size);
    lambda=time_step_size/step_size;
    Eigen::VectorXd x(n);
    Eigen::VectorXd x1(n);

    x(0)=0.0;
    for(int i = 1; i<n-1; i++) {
        x(i) = 0.1*sin(M_PI*((double)i/(double)n));
    }
    x(n-1)=0.0;

    //do_first_step(x, x1);
    x1=x;
    //cout<<"first step:"<<endl;

    for(double i = 2*time_step_size; i<=3.00001; i+=time_step_size) {
        //cout<<"t: "<<i<<"\n";
        do_step(x, x1);
        //cout<<"Vec:\n"<<x<<endl;
        for(int c = 0; c<n-1; c++) {
            cout<<x1(c)<<", ";
        }
        cout<<x1(n-1);
        cout<<"\n";
    }

    Eigen::VectorXd exact(n); 

    exact(0)=0;
    for(int i = 1; i<n-1; i++) {
        exact(i)=0.1*cos(M_PI)*sin(((double)i/(double)n)*M_PI);
    }
    exact(n-1)=0;   
    //cout<<"Exact: \n"<<exact<<endl;

    double sum=0;
    double cnt=0;
    for(int i = 1; i<(n-1); i++) {
        sum+=abs(exact(i)-x(i));
        cnt++;
    }
    //cout<<"Mean of local error: "<<sum/cnt<<endl;


}