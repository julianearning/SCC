#include <iostream>
#include <math.h>
#include <vector>
#include <string>

using namespace std;

int explicit_euler(double a, double b, double c, double d, double h) {
    double factor = 1.2;
    double b0 = factor*(d/c);
    double r0 = factor*(a/b);
    double curr_b=0.0;
    double curr_r=0.0;

    cout<<"$data <<EOD\n";
    for(int i=0; i<100; i++) {
        curr_b = b0 + h * ((a*b0)-(b*r0*b0));
        curr_r = r0 + h * ((c*b0*r0)-(d*r0));
        b0 = curr_b;
        r0 = curr_r;
        cout<<i<<" "<<curr_b<<" "<<curr_r<<"\n";
    }
    cout<<"EOD\n";
    cout<<"set xrange [0:100]\nset yrange [0:300]\nplot $data using 1:2 with linespoints title 'B', '' using 1:3 with linespoints t 'R'\n";

    return 0;
}

int heun(double a, double b, double c, double d, double h) {
    double factor = 1.2;
    double b0 = factor*(d/c);
    double r0 = factor*(a/b);
    double curr_b=b0;
    double curr_r=r0;
    double Bk1=0.0;
    double Rk1=0.0;
    double Bestimate=0.0;
    double Restimate=0.0;
    double Bk2=0.0;
    double Rk2=0.0;

    cout<<"$data <<EOD\n";
    for(int i=0; i<100; i++) {
        Bk1=(a*curr_b)-(b*curr_r*curr_b);
        Rk1=(c*curr_b*curr_r)-d*curr_r;
        Bestimate=curr_b+h*Bk1;
        Restimate=curr_r+h*Rk1;
        Bk2=(a*Bestimate)-(b*Restimate*Bestimate);
        Rk2=(c*Bestimate*Restimate)-d*Restimate;
        curr_b=curr_b+(h/2)*(Bk1+Bk2);
        curr_r=curr_r+(h/2)*(Rk1+Rk2);
        cout<<i<<" "<<curr_b<<" "<<curr_r<<"\n";
    }
    cout<<"EOD\n";
    cout<<"set xrange [0:100]\nset yrange [0:300]\nplot $data using 1:2 with linespoints title 'B', '' using 1:3 with linespoints t 'R'\n";

    return 0;
}


int mittelwert(double a, double b, double c, double d, double h) {
    double factor = 1.2;
    double b0 = factor*(d/c);
    double r0 = factor*(a/b);
    double curr_b=b0;
    double curr_r=r0;
    double Bk1=0.0;
    double Rk1=0.0;
    double Bestimate=0.0;
    double Restimate=0.0;
    double Bk2=0.0;
    double Rk2=0.0;

    cout<<"$data <<EOD\n";
    for(int i=0; i<100; i++) {
        Bk1=(a*curr_b)-(b*curr_r*curr_b);
        Rk1=(c*curr_b*curr_r)-d*curr_r;
        Bestimate=curr_b+(h/2)*Bk1;
        Restimate=curr_r+(h/2)*Rk1;
        Bk2=(a*Bestimate)-(b*Restimate*Bestimate);
        Rk2=(c*Bestimate*Restimate)-d*Restimate;
        curr_b=curr_b+h*Bk2;
        curr_r=curr_r+h*Rk2;
        cout<<i<<" "<<curr_b<<" "<<curr_r<<"\n";
    }
    cout<<"EOD\n";
    cout<<"set xrange [0:100]\nset yrange [0:300]\nplot $data using 1:2 with linespoints title 'B', '' using 1:3 with linespoints t 'R'\n";
    return 0;
}



int runge_kutta(double a, double b, double c, double d, double h) {
    double factor = 1.2;
    double b0 = factor*(d/c);
    double r0 = factor*(a/b);
    double curr_b=b0;
    double curr_r=r0;
    double Bk1=0.0;
    double Rk1=0.0;
    double Bestimate=0.0;
    double Restimate=0.0;
    double Bestimate2=0.0;
    double Restimate2=0.0;
    double Bestimate3=0.0;
    double Restimate3=0.0;
    double Bk2=0.0;
    double Rk2=0.0;
    double Bk3=0.0;
    double Rk3=0.0;
    double Bk4=0.0;
    double Rk4=0.0;

    cout<<"$data <<EOD\n";
    for(int i=0; i<100; i++) {
        Bk1=(a*curr_b)-(b*curr_r*curr_b);
        Rk1=(c*curr_b*curr_r)-d*curr_r;
        Bestimate=curr_b+(h/2)*Bk1;
        Restimate=curr_r+(h/2)*Rk1;
        Bk2=(a*Bestimate)-(b*Restimate*Bestimate);
        Rk2=(c*Bestimate*Restimate)-d*Restimate;
        Bestimate2=curr_b+(h/2)*Bk2;
        Restimate2=curr_r+(h/2)*Rk2;
        Bk3=(a*Bestimate2)-(b*Restimate2*Bestimate2);
        Rk3=(c*Bestimate2*Restimate2)-d*Restimate2;
        Bestimate3=curr_b+h*Bk3;
        Restimate3=curr_r+h*Rk3;
        Bk4=(a*Bestimate3)-(b*Restimate3*Bestimate3);
        Rk4=(c*Bestimate3*Restimate3)-d*Restimate3;
        curr_b=curr_b+h/6*(Bk1+2*Bk2+2*Bk3+Bk4);
        curr_r=curr_r+h/6*(Rk1+2*Rk2+2*Rk3+Rk4);
        cout<<i<<" "<<curr_b<<" "<<curr_r<<"\n";
    }
    cout<<"EOD\n";
    cout<<"set xrange [0:100]\nset yrange [0:300]\nplot $data using 1:2 with linespoints title 'B', '' using 1:3 with linespoints t 'R'\n";

    return 0;
}


int main(int argc, char *argv[]) {
    if(argc <= 2) {
        cout<<"Usage: ./uebung2.cpp [h] [method]"<<endl;
        return -1;
    }
    std::string method(argv[2]);
    if(method=="explicit_euler") {
        explicit_euler(0.3, 0.025, 0.0015, 0.2, atof(argv[1]));
    } else if(method=="heun") {
        heun(0.3, 0.025, 0.0015, 0.2, atof(argv[1]));
    } else if(method=="mittelwert") {
        mittelwert(0.3, 0.025, 0.0015, 0.2, atof(argv[1]));
    } else if(method=="runge_kutta") {
        runge_kutta(0.3, 0.025, 0.0015, 0.2, atof(argv[1]));
    }
    return 0;
}