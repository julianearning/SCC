#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <chrono>

using namespace std;

int gauss_elimination(Eigen::MatrixXd & old_A, Eigen::VectorXd & old_b, Eigen::VectorXd * result) {
    // abfangen, wenn es keine Lösung gibt, Pivotsuche fehlt
    double l=0.0;
    double sum=0.0;
    int n = old_A.rows();
    Eigen::MatrixXd A = old_A;  
    Eigen::VectorXd b = old_b;
    for(int i = 0; i<n; i++) {    
        for(int j = i+1; j<n; j++) {
            l = A(j,i) / A(i,i);
            for(int k = 0; k <n; k++) {
                A(j,k) = A(j,k) - l * A(i,k);    
            }
            b(j) = b(j) - l * b(i);
        }
    }
    
    for(int i = n-1; i>=0; i--) {
        sum = 0.0;
        for(int j = i+1; j < n; j++) {
            sum += (A(i,j)/A(i,i))*((*result)(j));
        }
        (*result)(i) = b(i)/A(i,i) - sum;
    }
    return 0;
}

int qr_zerlegung(Eigen::MatrixXd & old_A, Eigen::VectorXd & old_b, Eigen::VectorXd * result) {
    double temp=0.0;
    double rho=0.0;
    double cosphi=0.0;
    double minussinphi=0.0;
    double sum=0.0;
    int n = old_A.rows();
    Eigen::MatrixXd A = old_A;  
    Eigen::VectorXd b = old_b;
    Eigen::VectorXd z(n);
    Eigen::MatrixXd Q(n,n); 
    Eigen::MatrixXd R(n,n);
    for(int i = 0; i<n; i++) {   // initialize Q
        Q(i,i) = 1;
    }
    for(int i = 0; i<n; i++) {     // spalten
        for(int j = i+1; j<n; j++) {
            if(A(j,i) != 0) {      // if A(j,i)=0 do nothing
                if(A(i,i) == 0) {        // if A(i,i)=0 switch row j and i
                    for(int idx = 0; idx<n; idx++) {
                        temp=A(i, idx);
                        A(i,idx) = A(j,idx);
                        A(j,idx) = temp;
                    }
                } 

                rho=sqrt(pow(A(i,i),2) + pow(A(j,i),2));
                cosphi=A(i,i)/rho;
                minussinphi=A(j,i)/rho;
                

                if((i == 0) && (j==1)) { 
                    Q(i,i) = cosphi;
                    Q(j,j) = cosphi;
                    Q(i,j) = minussinphi; 
                    Q(j,i) = -minussinphi;
                } else {
                    for(int spalte=0; spalte<n; spalte++) {
                        Q(i,spalte) = cosphi * Q(i,spalte) + minussinphi * Q(j, spalte);
                        Q(j, spalte) = -minussinphi * Q(i, spalte) + cosphi * Q(j,spalte);
                    }
                }
            } 
        }
    }

    R = Q * A;

    z = Q * b;


    for(int i = n-1; i>=0; i--) {
        sum = 0.0;
        for(int j = i+1; j < n; j++) {
            sum += (R(i,j)/R(i,i))*((*result)(j));
        }
        (*result)(i) = z(i)/R(i,i) - sum;
    }


    return 0;
}


int cholesky_decomposition(Eigen::MatrixXd & old_A, Eigen::VectorXd & old_b, Eigen::VectorXd * result) {

    double sum=0.0;
    int n = old_A.rows();
    Eigen::MatrixXd A = old_A;  
    Eigen::VectorXd b = old_b;
    Eigen::MatrixXd L(n,n);
    Eigen::MatrixXd D(n,n);
    Eigen::MatrixXd R(n,n);
    Eigen::VectorXd y(n);


    for(int k = 0; k<n; k++) {
        sum=0.0;
        for(int j = 0; j<=(k-1); j++) {
            sum+=(pow(L(k,j),2) * D(j,j));
        }
        D(k,k)=A(k,k)-sum;
        L(k,k)=1;
        for(int i = (k+1); i<n; i++) {
            sum=0.0;
            for(int j=0; j<k; j++) {
                sum += (L(i,j) * D(j,j) * L(k,j));
            }
            L(i,k) = (A(i,k) - sum)/D(k,k);
        }
    } 


    // solve Ly=b    
    for(int i = 0; i<n; i++) {
        sum = 0.0;
        for(int j = 0; j < i; j++) {
            sum += (L(i,j)/L(i,i))*(y(j));
        }
        y(i) = b(i)/L(i,i) - sum;

    }


    R = D*L.transpose();


    // solve Rx=y
    for(int i = n-1; i>=0; i--) {
        sum = 0.0;
        for(int j = i+1; j < n; j++) {
            sum += (R(i,j)/R(i,i))*((*result)(j));
        }
        (*result)(i) = y(i)/R(i,i) - sum;
    }

    return 0;
}


bool residuum_almost_zero(Eigen::VectorXd & residuum) {
    double sum=0.0;
    for(int i=0; i<residuum.rows(); i++) {
        sum+=(double)residuum(i)*(double)residuum(i);
    }
    if(sum < 0.00000000001) {
        return true;
    } else {
        return false;
    }
}


int jacobi(Eigen::MatrixXd & old_A, Eigen::VectorXd & old_b, Eigen::VectorXd * result) {

    double sum=0.0;
    int k = 0;
    int n = old_A.rows();
    Eigen::MatrixXd A = old_A;  
    Eigen::VectorXd b = old_b;
    Eigen::VectorXd x(n);
    Eigen::VectorXd x2(n);
    Eigen::VectorXd residuum(n);

    for(int i = 0; i<n; i++) {
        x(i) = (lrand48()%10);
    }
    residuum = A*x-b; 

    while(!residuum_almost_zero(residuum)) {
        for(int j = 0; j<n; j++) {
            sum=0.0;
            for(int k = 0; k<j; k++) {
                sum += A(j,k)*x(k);
            }
            for(int k = (j+1); k<n; k++) {
                sum += A(j,k)*x(k);
            }
            x2(j) = (1/A(j,j)) * (b(j) - sum);
        }
        x = x2;
        residuum = A*x-b;
        k++;
    }

    *result = x; 

    return k;
}


int gauss_seidel(Eigen::MatrixXd & old_A, Eigen::VectorXd & old_b, Eigen::VectorXd * result) {
    double sum1=0.0;
    double sum2=0.0;
    int k=0;
    int n = old_A.rows();
    Eigen::MatrixXd A = old_A;  
    Eigen::VectorXd b = old_b;
    Eigen::VectorXd x(n);
    Eigen::VectorXd residuum(n);

    for(int i = 0; i<n; i++) {
        x(i) = (lrand48()%10);
    }

    residuum = A*x-b;

    while(!residuum_almost_zero(residuum)) {
        for(int j=0;j<n;j++) {
            sum1=0.0;
            sum2=0.0;
            for(int k = j+1; k<n;k++) {
                sum1+=(A(j,k)*x(k));
            }
            for(int k = 0; k<j;k++) {
                sum2+=(A(j,k)*x(k));
            }
            x(j) = ((b(j) - sum1) - sum2)/A(j,j);
        }
        residuum = A*x-b;
        k++;
    }

    *result = x; 

    return k;
}


int wikipedia_conjugate_gradient(Eigen::MatrixXd & old_A, Eigen::VectorXd & old_b, Eigen::VectorXd * result) {
    int n = old_A.rows();
    int k = 0;
    double alpha;
    double beta;
    Eigen::MatrixXd A = old_A;  
    Eigen::VectorXd b = old_b;
    Eigen::VectorXd x(n);
    Eigen::VectorXd p = old_b;
    Eigen::VectorXd r = old_b;
    Eigen::VectorXd r_old = old_b;

    for(int i = 0; i<n; i++) {
        x(i) = 0;
    }

    r_old = b-(A*x);
    if(residuum_almost_zero(r_old)) {
        for(int i = 0; i<n; i++) {
            (*result)(i) = x(i);
        }
        return k;
    }

    p=r;

    while(true) {
        alpha=r_old.dot(r_old)/p.dot(A*p);
        //std::cout<<(A*p)(3)<<std::endl;
        x=x+alpha*p;
        r=r-alpha*(A*p);
        if(residuum_almost_zero(r)) {
            for(int i = 0; i<n; i++) {
                (*result)(i) = x(i);
            }
            return k;
        }
        beta=r.dot(r)/r_old.dot(r_old);
        p=r+(beta*p);

        k++;

        r_old=r;
    }

}


int conjugate_gradient(Eigen::MatrixXd & old_A, Eigen::VectorXd & old_b, Eigen::VectorXd * result) {
    int n=old_A.rows();
    int k=0;
    Eigen::MatrixXd A = old_A;
    Eigen::VectorXd b = old_b;
    Eigen::VectorXd x(n);
    Eigen::VectorXd p = old_b;
    Eigen::VectorXd r = old_b;
    Eigen::VectorXd Ap_i_minus_one(n);
    double alpha=0.0;
    double beta=0.0; 
    double norm_r_i_minus_one_squared=0.0;
    double norm_r_i_squared=0.0;
    double scalar_product=0.0;

    for(int i = 0; i<n; i++) {
        x(i) = 0;
    }

    while(!residuum_almost_zero(r)) {
        norm_r_i_minus_one_squared = 0.0;
        norm_r_i_squared = 0.0;    
        Ap_i_minus_one = A * p;    // hier ist die Matrix mal Vektor Mult.
        scalar_product = 0.0;
        for(int i = 0; i<n; i++) {   
            norm_r_i_minus_one_squared += (r(i)*r(i));
            scalar_product += (p(i) * Ap_i_minus_one(i));
        }
        alpha = (norm_r_i_minus_one_squared) / scalar_product;

        for(int i = 0; i<n; i++) {
            x(i) = x(i) + (alpha * p(i));
            r(i) = r(i) - (alpha * Ap_i_minus_one(i));   
            norm_r_i_squared += (r(i)* r(i));
        }
        beta=norm_r_i_squared/norm_r_i_minus_one_squared;

        for(int i = 0; i<n; i++) {
            p(i) = r(i) + (beta * p(i));
        }
        k++;
    }
    for(int i = 0; i<n; i++) {
        (*result)(i) = x(i);
    }

    return k;
}



void uebung1(int n, bool debug) {
    int n_iterations=0;
    Eigen::MatrixXd A(n,n);
    Eigen::VectorXd b(n);
    Eigen::VectorXd solution1(n);
    Eigen::VectorXd solution2(n);
    Eigen::VectorXd solution3(n);
    Eigen::VectorXd solution4(n);
    Eigen::VectorXd solution5(n);
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    for(int i = 0; i<n; i++) {
        for(int j = 0; j<n; j++) {
            A(j,i) = (lrand48()%11)/(10.0 * (double)n);
        }
    }
    for(int i = 0; i<n; i++) {
        A(i,i) = (lrand48()%10 + 1);
        b(i) = (lrand48()%10 + 1);
    }
    

    begin = std::chrono::steady_clock::now();
    gauss_elimination(A,b,&solution1);
    end = std::chrono::steady_clock::now();

    std::cout << "Gauss Elimination = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s], " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<"[ms]"<< std::endl;

    begin = std::chrono::steady_clock::now();
    qr_zerlegung(A,b,&solution2);
    end = std::chrono::steady_clock::now();

    std::cout << "QR Decomp = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s], " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<"[ms]"<< std::endl;

    begin = std::chrono::steady_clock::now();
    n_iterations=jacobi(A,b,&solution3);
    end = std::chrono::steady_clock::now();

    std::cout << "Jacobi = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s], " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<"[ms], Iterations="<<n_iterations<<std::endl;

    begin = std::chrono::steady_clock::now();
    n_iterations=gauss_seidel(A,b,&solution4);
    end = std::chrono::steady_clock::now();

    std::cout << "Gauss-Seidel = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s], " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"[ms], Iterations="<<n_iterations<<std::endl;


    begin = std::chrono::steady_clock::now();
    n_iterations=conjugate_gradient(A,b,&solution5);
    end = std::chrono::steady_clock::now();

    std::cout << "Conjugate Gradient = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s], " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<"[ms], Iterations="<<n_iterations<< std::endl;

    if(debug) {
        //cout<<"-------\nA\n-------\n"<<A<<"\n";
        //cout<<"-------\nb\n-------\n"<<b<<"\n";
        cout<<"-------\nSolution Gauss-Elimination: \n"<<solution1((int)(n/2))<<"\n";
        cout<<"-------\nSolution QR-Decomposition: \n"<<solution2((int)(n/2))<<"\n";
        cout<<"-------\nSolution Jacobi-Iteration: \n"<<solution3((int)(n/2))<<"\n";
        cout<<"-------\nSolution Gauss-Seidel: \n"<<solution4((int)(n/2))<<"\n";
        cout<<"-------\nSolution Conjugate Gradient: \n"<<solution5((int)(n/2))<<"\n";
    }

}


void uebung2(int n, bool debug) {
    int n_iterations=0;
    double sum=0.0;
    Eigen::MatrixXd A(n,n);
    Eigen::VectorXd b(n);
    Eigen::VectorXd solution1(n);
    Eigen::VectorXd solution2(n);
    Eigen::VectorXd solution3(n);
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    for(int i = 0; i<n; i++) {
        for(int j=0; j<n; j++) {
            A(i,j) = (double)1/((double)(i+1)+(double)(j+1)-1);
        }
    }

    for(int i = 0; i<n; i++) {
        sum=0.0;
        for(int j=0; j<n; j++) {
            sum+=A(i,j);
        }
        b(i) = sum;
    }

    begin = std::chrono::steady_clock::now();
    gauss_elimination(A,b,&solution1);
    end = std::chrono::steady_clock::now();
    
    std::cout << "Gauss-Elimination = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s], " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<"[ms]"<< std::endl;

    begin = std::chrono::steady_clock::now();
    cholesky_decomposition(A,b,&solution2);
    end = std::chrono::steady_clock::now();

    std::cout << "Cholesky Decomposition = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s], " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<"[ms]"<< std::endl;
    
    begin = std::chrono::steady_clock::now();
    n_iterations=conjugate_gradient(A,b,&solution3);
    end = std::chrono::steady_clock::now();

    std::cout << "Conjugate Gradients = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s], " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<"[ms], Iterations="<<n_iterations<< std::endl;


    if(debug) {
        //cout<<"-------\nA\n-------\n"<<A<<"\n";
        //cout<<"-------\nb\n-------\n"<<b<<"\n";
        cout<<"-------\nSolution Gauss-Elimination: \n"<<solution1((int)(n/2))<<"\n";
        cout<<"-------\nSolution Cholesky Decomposition: \n"<<solution2((int)(n/2))<<"\n";
        cout<<"-------\nSolution Conjugate Gradients: \n"<<solution3((int)(n/2))<<"\n";
    }
}


int main(int argc, char *argv[]) {
    if(argc <= 2) {
        cout<<"Usage: ./uebung1.cpp [uebungs_nr] [n] [debug true/" "]"<<endl;
        return -1;
    }
    int uebung=stoi(argv[1]);
    if((uebung!=1) && (uebung!=2)) {
        cout<<"uebung "<<uebung<<" existiert nicht"<<endl;
    }
    int n = stoi(argv[2]);
    bool debug=false;
    
    if(argc == 4) {
        debug=true;
    }

    if(uebung==1) {
        cout<<"uebung1"<<" n="<<n<<endl;
        uebung1(n, debug);
    } else if(uebung==2) {
        cout<<"uebung2"<<" n="<<n<<endl;
        uebung2(n,debug);
    }

    return 0;
}  