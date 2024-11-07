#include <omp.h>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <chrono>

using namespace std;


struct Range {
    double start;
    double end;
};
 
struct Range columns {0,0.5};
struct Range rows {0,0.5};

int n;   // n Zeilen 
int m;   // m Spalten
double omega;


double calc_residuum(Eigen::MatrixXd & A) {
    // M*x-b   M: 4 auf der Diagonalen -1 auf der Nebendiagonalen und über
    Eigen::VectorXd x(n*m);
    Eigen::VectorXd Mx(n*m);
    for(int i = 0; i<n; i++) {
        for(int j=0; j<m; j++) {
            x(i*n+j)=A(i,j);
        }
    }
    // M*x
    Mx(0)=0.0;
    double val;
    for(int i = 1; i<(n-1); i++) {
        for(int j = 1; j<(m-1); j++) {
            Mx(i*n+j)=4*x(i*n+j)-x(i*n+j-1)-x(i*n+j+1)-x(i*n+j-n)-x(i*n+j+n);
            //Mx(i*n+j)=4*x(i*n+j)-x(i*n+j-1)-x(i*n+j+1)-x(i*n+j-n)-x(i*n+j+n);
        }
    }
    Mx((n*m)-1)=0.0;

    // Mx - b
    double sum_mx = 0.0;
    double sum_x = 0.0;
    for(int i = 1; i<(n-1); i++) {
        for(int j = 1; j<(m-1); j++) {
            sum_mx+=abs(Mx(i*n+j));
            sum_x+=abs(x(i*n+j));
        }
    }
    return sum_mx/sum_x;
}


int jacobi(Eigen::MatrixXd & init_A, Eigen::MatrixXd * result) {
    int n_iterations=0;
    double val;
    Eigen::MatrixXd A2(m,n);
    Eigen::MatrixXd A(m,n);
    A=init_A;
    A2=init_A;
    while(calc_residuum(A) > 10e-10) {
        for(int i = 1; i<n-1; i++) {
            for(int j=1; j<m-1; j++) {
                val = 0.25*(A(i-1,j)+A(i+1,j)+A(i,j+1)+A(i,j-1));
                A2(i,j) = val;
            }
        }
        A = A2;
        n_iterations++;
    }
    *result = A;
    return n_iterations;
}


int gauss_seidel(Eigen::MatrixXd & init_A, Eigen::MatrixXd * result) {
    int n_iterations=0;
    Eigen::MatrixXd A(m,n);
    A=init_A;
    while(calc_residuum(A) > 10e-10) {
        for(int i = 1; i<(n-1); i++) {
            for(int j=1; j<(m-1); j++) {
                A(i,j) = 0.25*(A(i-1,j)+A(i+1,j)+A(i,j+1)+A(i,j-1));
            }
        }
        n_iterations++;
    }
    *result = A;
    return n_iterations;
}



int SOR(Eigen::MatrixXd & init_A, Eigen::MatrixXd * result) {
    int n_iterations=0;
    Eigen::MatrixXd A(m,n);
    A=init_A;
    while(calc_residuum(A) > 10e-10) {
        for(int i = 1; i<(n-1); i++) {
            for(int j=1; j<(m-1); j++) {
                //A(i,j) = 0.25*(A(i-1,j)+A(i+1,j)+A(i,j+1)+A(i,j-1));
                A(i,j)=(1-omega)*A(i,j)+0.25*omega*(A(i-1,j)+A(i+1,j)+A(i,j-1)+A(i,j+1));
            }
        }
        n_iterations++;
    }
    *result = A;
    return n_iterations;
}


double check_local_error(Eigen::MatrixXd & B, int i, int j) {
    return abs((400*(((double)i/(double)n)*(columns.end-columns.start)+columns.start)*(((double)j/(double)m)*(rows.end-rows.start)+rows.start))- B(i,j));
}


int main(int argc, char *argv[]) {

    if(argc < 2) {
        cout<<"Usage: ./uebung5 [step_size (for example 0.001)]\n";
        return -1;
    }

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    int n_iterations;
    double step_size=atof(argv[1]);

    // m=Anzahl Zeilen
    m=(int)((columns.end-columns.start)/step_size);

    // n=Anzahl Spalten
    n=(int)((rows.end-rows.start)/step_size);

    omega=2/(1+sqrt(1-cos(M_PI/n)*cos(M_PI/n)));

    cout<<"Matrix: "<<m<<"x"<<n<<endl;

    Eigen::MatrixXd A(m,n);

    // für alle Spalten: 
    for(int i=0; i<n; i++) {
        A(0,i)=0;
        A(n-1,i) = 200*((double)i/(double)n)*(columns.end-columns.start)+columns.start;
    }
    for(int i=0; i<m; i++) {
        A(i,0)=0;
        A(i,m-1) = 200*((double)i/(double)m)*(rows.end-rows.start)+rows.start;
    }


    //double init_value=100*(rows.end-rows.start)+rows.start;
    double init_value=0.0;
    for(int i=1; i<(n-1); i++) {
        for(int j=1; j<(m-1); j++) {
            A(i,j) = init_value;
            //A(i,j)=400*(((double)i/(double)n)*(columns.end-columns.start)+columns.start)*(((double)j/(double)m)*(rows.end-rows.start)+rows.start);
        }
    }
    //cout<<A<<endl;


    Eigen::MatrixXd B(m,n);
    

    begin = std::chrono::steady_clock::now();
    n_iterations=jacobi(A, &B);
    end = std::chrono::steady_clock::now();
    std::cout << "t = " << chrono::duration_cast<chrono::seconds>(end - begin).count() << "[s]"<< endl;
    
    cout<<"Jacobi used n_iterations: "<<n_iterations<<endl;
    //cout<<B<<endl;

    double sum=0;
    int cnt=0;
    for(int i = 1; i<(n-1); i++) {
        for(int j=1; j<(m-1); j++) {
            sum+=check_local_error(B, i, j);
            cnt++;
        }
    }
    cout<<"Mean of local error: "<<sum/cnt<<endl;


    n_iterations=gauss_seidel(A, &B);

    cout<<"Gauss-Seidel used n_iterations: "<<n_iterations<<endl;
    //cout<<B<<endl;

    sum=0;
    cnt=0;
    for(int i = 1; i<(n-1); i++) {
        for(int j=1; j<(m-1); j++) {
            sum+=check_local_error(B, i, j);
            cnt++;
        }
    }
    cout<<"Mean of local error: "<<sum/cnt<<endl;

    n_iterations=SOR(A,&B);

    cout<<"SOR used n_iterations: "<<n_iterations<<endl;
    //cout<<B<<endl;


    sum=0;
    cnt=0;
    for(int i = 1; i<(n-1); i++) {
        for(int j=1; j<(m-1); j++) {
            sum+=check_local_error(B, i, j);
            cnt++;
        }
    }
    cout<<"Mean of local error: "<<sum/cnt<<endl;


    return 0;
}