#include "typedef.h"

using namespace std;
using namespace std::chrono;


struct poly{
    int n;
    vd p;
    poly(int _n, vd _p){
        n = _n;
        p = _p;
    }

    poly(int _n){
        n = _n;
        p.assign(_n, 0);
    }
};


/// Operadores de Matrices y vectores///
double sum(vd v);
double norma(vd v);
void normalize(vd &x);
vvd Transpose(vvd M);
double Det(vvd A);
vvd Inverse(vvd A);
vvd Inverse_PentaD(vvd A);
vvd Identity(int n);
void Identity(vvd &M);
double factorial(int n);

///Factorizaciones
void LU_Croud(vvd &A, vvd &L, vvd &U);
void Cholesky(vvd &A, vvd &H);
void Cholesky_TriD(vvd &A, vvd &H);
void Cholesky_PentaD(vvd &A, vvd &H);
void Factorizacion_QR(vvd &A, vvd &Q, vvd &R);


///Resuelve Ecuaciones///
void Diagonal(vvd &A, vd &b, vd &x);
void Triangular_Superior(vvd &U, vd &b, vd &x);
void Triangular_Inferior(vvd &L, vd &b, vd &x);
void Eliminacion_Gaussiana(vvd &A, vd &b, vd &x);
void Partial_pivot(vvd &A, int k, vd &b);
void solveEC(vd &x, double Q, double K, double phi_0, double phi_n, double n, double L);
void Jacobi(vvd &A, vd &b, vd &x, double TOL, int max_it);
void Gauss_Seidel(vvd &A, vd &b, vd &x, double TOL, int max_it);
void metodo_Potencia(vvd A, vvd &eigenvectors, vd &eigenvalues, double error, int range, int max_it);
void metodo_Potencia_Inv(vvd &A, vvd &eigenvectors, vd &eigenvalues, double tol, int range, int max_it);
void Jacobi_Eigen_Values(vvd &A, vd &eigenvalues, double tol, int max_it);
void metodo_Subespacio_Pot(vvd &A, vvd &eigenvectors, vd &eigenvalues, double tol, int max_it, int m);
void metodo_Subespacio_Pot_Inv(vvd &A, vvd &eigenvectors, vd &eigenvalues, double tol, int max_it, int m);
void cociente_Rayleigh(vvd &A, vd &v, double &lambda, double tol, int max_it);
void Gradiente_Conjudado(vvd &A, vd &b, vd &x, vvd &M, bool flag, double tol, int max_it);

/// Interpolacion
void Minimos_Cuadrados(vd x, vd y, int m, vd &coef);
void Minimos_Cuadrados_Generalizado(vvd X, vd y, vd &coef);
void Matricial_Interpolation(vd x, vd y, vd &coef);
poly Lagrange_Interpolation(vd x, vd y);
double Lagrange_Interpolation(vd x, vd y, double z);
double Newton_Interpolation(vd x, vd y, double z);
void cubic_Splines(vd x, vd y, vvd &coef);
void quadratic_Splines(vd x, vd y, vvd &coef, double diff);
void lineal_Splines(vd x, vd y, vvd &coef);

/// Integracion Numerica
double mcmc_integration_2d(double a, double b);
double mcmc_integration_3d(double a, double b, double c, double d, double (*f)(double, double));
double Regla_Trapecios(double a, double b, int n, double (*f)(double));
double Cuadratura_Gaussiana(double a, double b, int nodos, double (*f)(double));
double Romberg(double a, double b, double eps, double (*f)(double));
double Extra_Richardson(double a, double b, int n, double (*f)(double));
double Newton_Cotes(double a, double b, int n, double (*f)(double));

///Diferenciacion Numerica///
double Diferencias_Finitas(double x , double (*f)(double), int n, double h);
double Lagrange_Derivative(vd x, vd y, double z, int m);
double Extra_Richardson_diff(double x, double (*f)(double), int n, double h);

/// Funciones Proyecto Final ///
void PCA(vvd X, vvd &Y, int n_comp, double &contribution);
double Softmax(vd x, int k, vvd beta);
void Softmax_Estimation(vvd X, vd y, vvd &beta, double alpha, int max_it);
void Softmax_Predictions(vvd new_x, vvd beta, vd &pred);
double Accuracy(vd pred, vd real);
double Balanced_Accuracy(vd pred, vd real);


///Operadores de Matrices y vectores///

poly operator+ (poly a, poly b){
    if(a.n < b.n)swap(a,b);
    int m = b.n;
    for(int i = 0; i<m; ++i)
        a.p[i] += b.p[i];
    return a;
}

poly operator* (poly a, poly b){
    int n = a.n, m = b.n;
    poly c(n+m,vd(n+m,0));
    for(int i = 0; i<n; ++i)
        for(int j = 0; j<m; ++j)
            c.p[i+j] += a.p[i]*b.p[j];
    return c;
}

poly operator- (poly a, ld x){
    a.p[0]-=x;
    return a;
}

poly operator/ (poly a, ld x){
    int n = a.n;
    for(int i = 0; i<n; ++i)
        a.p[i]/=x;
    return a;
}

poly operator* (poly a, ld x){
    int n = a.n;
    for(int i = 0; i<n; ++i)
        a.p[i]*=x;
    return a;
}

double factorial(int n){

    if (n == 0)
        return 1;
    else
        return n*factorial(n-1);
}

double sum(vd v){

    int n = v.size();

    double s = 0;
    for(int i = 0; i < n; i++){
        s += v[i];
    }

    return s;
}

vvd Identity(int n) {

    vvd M(n, vd(n, 0));

    for (int i = 0; i < n; i++) {
        M[i][i] = 1;
    }

    return M;
}

void Identity(vvd &M) {
    int n = M.size();
    int m = M[0].size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                M[i][i] = 1;
            else
                M[i][j] = 0;
        }
    }
}

double norma(vd v) {

    int n = v.size();
    double suma = 0;

    for (int i = 0; i < n; i++)
        suma += v[i] * v[i];

    return sqrt(suma);
}

double norma(vvd M) {

    int n = M.size();
    int m = M[0].size();
    double suma = 0;

    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            suma += M[i][j] * M[i][j];
        }
    }

    return sqrt(suma);
}

void normalize(vd &x){

    double _x = norma(x);
    for(int i = 0; i < (int)x.size(); ++i)
        x[i] /= _x;
}

double operator* (vd a, vd b){

    int n = a.size();
    double ans = 0;

    for(int i = 0; i < n; ++i)
        ans += a[i]*b[i];

    return ans;
}

vd operator* (double l, vd a){

    int n = a.size();
    for(int i = 0; i < n; ++i)
        a[i] *= l;

    return a;
}

vd operator* (vd a, double l){

    int n = a.size();
    for(int i = 0; i < n; ++i)
        a[i] *= l;

    return a;
}

vvd operator* (double k, vvd A){

    int n = A.size();
    int m = A[0].size();

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++i){
            A[i][j] *= k;
        }
    }

    return A;
}

vd operator* (vvd A, vd a){

    int n = A.size();
    int m = a.size();
    double acc;
    vd v_r(m, 0);

    if (n != m) {
        cout << "ERROR: size of matrix and vector dont match" << endl;
        return v_r;
    }
    for (int i = 0; i < n; i++) {

        acc = 0;
        for (int j = 0; j < n; j++)
            acc += A[i][j] * a[j];

        v_r[i] = acc;
    }

    return v_r;
}

vd operator* (vd a, vvd A){

    int n = A.size();
    int m = a.size();
    double acc;
    vd v_r(m, 0);

    if (n != m) {
        cout << "ERROR: size of matrix and vector dont match" << endl;
        return v_r;
    }
    for (int i = 0; i < n; i++) {

        acc = 0;
        for (int j = 0; j < n; j++)
            acc += a[j]*A[j][i];

        v_r[i] = acc;
    }

    return v_r;
}

vvd operator* (vvd A, vvd B){

    int n = A.size();
    int m = B[0].size();
    int p = A[0].size();

    vvd C(n, vd(m, 0));

    //omp_set_num_threads(omp_get_num_procs());

    int i, j, k;

    //#pragma omp parallel for private(i,j,k) shared(A,B,C,n, m, p) default(none)
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            double acc = 0;

            for (k = 0; k < p; k++) {
                acc += A[i][k] * B[k][j];
            }
            C[i][j] = acc;
        }
    }

    return C;
}

vd operator- (vd a, vd b){

    int n = a.size();

    for(int i = 0; i < n; ++i)
        a[i] -= b[i];

    return a;
}

vvd operator- (vvd A, vvd B){

    int n = A.size();
    int m = A[0].size();

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            A[i][j] -= B[i][j];
        }
    }

    return A;
}

vd operator- (vd a, double l){

    int n = a.size();

    for(int i = 0; i < n; ++i)
        a[i] -= l;

    return a;
}

vd operator+ (vd a, vd b){

    int n = a.size();

    for(int i = 0; i < n; ++i)
        a[i] += b[i];

    return a;
}

vvd operator+ (vvd A, vvd B){

    int n = A.size();
    int m = A[0].size();

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            A[i][j] += B[i][j];
        }
    }

    return A;
}

///Matriz vector multiplicacion suponiendo que A es tridiagonal
vd operator/ (vvd A, vd a) {

    int n = A.size();
    int m = a.size();
    double acc;
    vd v_r(m, 0);

    if (n != m) {
        cout << "ERROR: size of matrix and vector dont match" << endl;
        return v_r;
    }
    for (int i = 0; i < n; i++) {

        acc = 0;
        for (int j = max(0, i - 2); j < min(n, i + 3); j++){
            acc += A[i][j] * a[j];
        }
        v_r[i] = acc;
    }

    return v_r;
}

vd operator/ (vd a, vvd A){

    int n = A.size();
    int m = a.size();
    double acc;
    vd v_r(m, 0);

    if (n != m) {
        cout << "ERROR: size of matrix and vector dont match" << endl;
        return v_r;
    }
    for (int i = 0; i < n; i++) {

        acc = 0;
        for (int j = max(0, i - 2); j < min(n, i + 3); j++)
            acc += a[j]*A[j][i];

        v_r[i] = acc;
    }

    return v_r;
}


///Matriz vector multiplicacion suponiendo que A es tridiagonal


vvd Transpose(vvd M) {

    int n = M.size(), m = M[0].size();

    vvd M_t(m, vd(n, 0));

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            M_t[j][i] = M[i][j];
        }
    }

    return M_t;
}

double Det(vvd A){

    int n = A.size();
    double det = 1;

    for (int i = 0; i < n; i++) {

        double pivotElement = A[i][i];

        int pivotRow = i;
        for (int row = i + 1; row < n; row++) {
            if (abs(A[row][i]) > abs(pivotElement)) {
                pivotElement = A[row][i];
                pivotRow = row;
            }
        }

        if (fabs(pivotElement) < 0.0000000001) {
            return 0.0;
        }

        if (pivotRow != i) {
            A[i].swap(A[pivotRow]);
            det *= -1.0;
        }

        det *= pivotElement;

        for (int row = i + 1; row < n; row++) {
            for (int col = i + 1; col < n; col++) {
                A[row][col] -= A[row][i] * A[i][col] / pivotElement;
            }
        }

    }

    return det;
}

vvd Inverse(vvd A){

    int n = A.size();

    vvd I(n, vd(n, 0));

    double aux;
    double pivotElem;

    for(int i = 0; i < n; i++){
        I[i][i] = 1;
        pivotElem = A[i][i];


        for(int k = 0; k < n; k++){
            A[i][k] /= pivotElem;
            I[i][k] /= pivotElem;
        }

        for(int j = 0; j < n; j++){
            if(i != j){
                aux = A[j][i];
                for(int k = 0; k < n; k++){
                    A[j][k] -= aux * A[i][k];
                    I[j][k] -= aux * I[i][k];
                }
            }
        }

    }

    return I;
}

vvd Inverse_PentaD(vvd A){

    int n = A.size();

    vvd Idem(n, vd(n, 0));
    vvd Ch(n, vd(n, 0));
    vvd I(n, vd(n, 0));

    for(int i = 0; i < n; i++){
        Idem[i][i] = 1;
    }

    Idem = Transpose(Idem);

    Cholesky_PentaD(A, Ch);

    for(int i = 0; i < n; i++){
        Triangular_Inferior(Ch, Idem[i], I[i]);
        Triangular_Superior(Ch, I[i], I[i]);
    }

    return I;
}

double max_element(vvd &A, int &i_max, int &j_max){

    int n = A.size();
    int i, j;
    vd max(n, 0.0);
    vector<int> max_j(n, 0);

    //omp_set_num_threads(omp_get_num_procs());
    //#pragma omp parallel for private(i, j) shared(A, max, max_j, n) default(none)

    for (i = 0; i < n; i++){
        for (j=0; j < n; j++){
            if (i == j)
                continue;

            if(max[i] < fabs(A[i][j])){
                max[i] = fabs(A[i][j]);
                max_j[i] = j;
            }
        }
    }

    double realMax = 0;
    for (int k = 0; k < n; k++){
        if (max[k] > realMax){
            realMax = max[k];
            i_max = k;
            j_max = max_j[k];
        }
    }
    return realMax;
}

void grand_schmidt(vvd &vv, vd &v, int nv) {

    int n = v.size();

    for (int j = 0; j < nv; j++) {
        double a_k = vv[j]*v;
        for (int k = 0; k < n; k++){
            v[k] = v[k] - a_k*vv[j][k];
        }
    }

    return;
}


///Factorizaciones///

void factorizarLDLT(vvd &A, vvd &L, vvd &DLT) {
    int n = A.size();

    vd D;
    D.assign(n, 0.0);

    for (int i = 0; i < n; i++) {


        L[i][i] = 1;
        //L_ij
        for (int j = 0; j <= i; j++) {

            double acc = 0;
            // D_jj
            for (int k = 0; k < j; k++) {
                acc += L[j][k] * L[j][k] * D[k];
            }
            D[j] = A[j][j] - acc;


            acc = 0;

            if (i > j) {
                //h_ij
                for (int k = 0; k < j; k++) {
                    acc += L[i][k] * L[j][k] * D[k];
                }
                L[i][j] = (A[i][j] - acc) / D[j];
            }
        }
    }
    //calcular DLT
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            DLT[i][j] = L[j][i];
            L[j][i] = L[j][i] * D[i];
        }
    }
}

void LU_Croud(vvd &A, vvd &L, vvd &U) {

    int n = A.size();

    for (int i = 0; i < n; i++) {
        U[i][i] = 1;
    }

    for (int k = 0; k < n; k++) {
        double acc;

        for (int i = k; i < n; i++) {
            acc = 0;

            for (int r = 0; r <= k - 1; r++) {
                acc += L[i][r] * U[r][k];
            }
            L[i][k] = A[i][k] - acc;
        }

        if (L[k][k] == 0)
            cout << "LU Decomposition  Msg WARNING: L_[" << k << ", " << k << "] es nulo, se dividira por 0" << endl;

        for (int j = k + 1; j < n; j++) {
            acc = 0;
            for (int r = 0; r <= k - 1; r++) {
                acc += L[k][r] * U[r][j];
            }
            U[k][j] = (A[k][j] - acc) / L[k][k];
        }
    }

    return;
}


void Cholesky(vvd &A, vvd &H) {

    int n = A.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double acc = 0;

            if (i != j) {

                for (int k = 0; k < j; k++) {
                    acc += H[i][k] * H[j][k];
                }
                H[i][j] = (A[i][j] - acc) / H[j][j];
            }
        }

        double acc = 0;

        for (int k = 0; k < i; k++) {
            acc += H[i][k] * H[i][k];
        }

        H[i][i] = sqrt(A[i][i] - acc);

    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < min(n, i + 3); j++) {
            H[i][j] = H[j][i];
        }
    }

    return;
}


void Cholesky_TriD(vvd &A, vvd &H) {

    int n = A.size();

    for (int i = 0; i < n; i++) {

        if (i != 0)
            H[i][i - 1] = A[i][i - 1] / H[i - 1][i - 1];

        H[i][i] = sqrt(A[i][i] - H[i][i - 1] * H[i][i - 1]);
    }

    for (int i = 0; i < n-1; i++) {
        H[i][i + 1] = H[i + 1][i];
    }

    return;
}

void Cholesky_PentaD(vvd &A, vvd &H) {

    int n = A.size();

    for (int i = 0; i < n; i++) {
        for (int j = max(0, i - 2); j < i; j++) {
            double acc = 0;

            if (i != j) {

                for (int k = max(0, i - 2); k < j; k++) {
                    acc += H[i][k] * H[j][k];
                }
                H[i][j] = (A[i][j] - acc) / H[j][j];
            }
        }

        double acc = 0;

        for (int k = max(0, i - 2); k < i; k++) {
            acc += H[i][k] * H[i][k];
        }

        H[i][i] = sqrt(A[i][i] - acc);

    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            H[i][j] = H[j][i];
        }
    }

    return;
}


void Factorizacion_QR(vvd &A, vvd &Q, vvd &R){

    int n = A.size();
    int m = A[0].size();
    vvd V = Transpose(A);

    //inicializar Q, R
    Q.assign(m, vd(n ,0.0));
    R.assign(n, vd(n ,0.0));

    // inicializar Q

    for (int i = 0; i < n; i++){
        R[i][i] = norma(V[i]);
        Q[i] = V[i]*(1/R[i][i]);

        for(int j = i; j < n; j++){
            R[i][j] = Q[i]*V[j];
            V[j] = V[j] - R[i][j]*Q[i];
        }
    }

    Q = Transpose(Q);

    return;
}


///Resuelve Ecuaciones///


void Diagonal(vvd &A, vd &b, vd &x) {

    int n = A.size();

    for (int i = 0; i < n; i++) {
        if (A[i][i] == 0)
            cout << "Diag WARNING: Elemento nulo en posicion [" << i << "][" << i << "]" << endl;

        x[i] = b[i] / A[i][i];
    }

    return;
}

void Triangular_Superior(vvd &U, vd &b, vd &x) {

    int n = U.size();
    double X[n];

    for (int i = n - 1; i >= 0; i--) {
        double acc = 0;
        for (int j = min(n - 1, i + 2); j > i; j--) {
            acc += U[i][j] * x[j];
        }
        if (U[i][i] == 0){
            cout << "Triangular Sup  Msg WARNING: A_[" << i << ", " << i << "] es nulo, se dividira por 0" << endl;
            exit(0);
        }

        x[i] = (b[i] - acc) / U[i][i];
    }

    return;
}

void Triangular_Inferior(vvd &L, vd &b, vd &x) {

    int n = L.size();

    for (int i = 0; i < n; i++) {
        double acc = 0;
        for (int j = max(0, i - 2); j < i; j++) {
            acc += L[i][j] * x[j];
        }
        if (L[i][i] == 0){
            cout << "Triangular Inf  Msg WARNING: A_[" << i << ", " << i << "] es nulo, se dividira por 0" << endl;
            exit(0);
        }

        x[i] = (b[i] - acc) / L[i][i];
    }

    return;
}

void Partial_pivot(vvd &A, int k, vd &b) {

    int n = A.size();
    int i_max = k;
    double maxim = fabs(A[k][k]);

    for (int i = k + 1; i < n; i++) {
        if (fabs(A[i][k]) > maxim) {
            i_max = i;
            maxim = fabs(A[i][k]);
        }
    }

    double temp;

    if (i_max != k) {
        for (int j = 0; j < n; j++) {
            temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
        temp = b[k];
        b[k] = b[i_max];
        b[i_max] = temp;
    }

    return;
}

void Eliminacion_Gaussiana(vvd &A, vd &b, vd &x) {

    int n = A.size();

    for (int k = 0; k < n - 1; k++) {

        if (A[k][k] == 0)
            Partial_pivot(A, k, b);

        if (A[k][k] == 0)
            cout << "GEM  Msg WARNING: A_[" << k << ", " << k << "] es nulo, se dividira por 0" << endl;

        for (int i = k + 1; i < n; i++) {
            double m_ik = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - m_ik * A[k][j];
            }
            b[i] = b[i] - m_ik * b[k];
        }
    }

    Triangular_Superior(A, b, x);

    return;
}


void solveEC_Cholesky(vd &x, double Q, double K, double phi_0, double phi_n, double n, double L ){

    double cf = (Q * (L/n)*(L/n))/K;

    vvd M(n, vd(n, 0));
    vvd H(n, vd (n, 0));
    vd b(n, cf);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j)
                M[i][j] = 2;

            if (i == j - 1 || i == j + 1)
                M[i][j] = -1;
        }
    }

    b[0] = cf + phi_0;
    b[n - 1] = cf + phi_n;

    Cholesky_TriD(M, H);

    Triangular_Inferior(H, b, x);
    Triangular_Superior(H, x, x);

    return;
}

void Jacobi(vvd &A, vd &b, vd &x, double TOL, int max_it) {

    int n = b.size();
    vd xn(n, 0);

    for (int i = 0; i < n; i++) {
        x[i] = b[i] / A[i][i];
    }

    for (int it = 0; it < max_it; it++) {
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                if (i == j)
                    continue;

                sum += A[i][j] * x[j];
            }
            xn[i] = (b[i] - sum) / A[i][i];
        }

        double num = 0;
        double den = 0;
        double fact;

        for (int i = 0; i < n; i++) {
            double dif_v_num = (xn[i] - x[i]);
            double dif_v_den = xn[i] * xn[i];
            num += dif_v_num * dif_v_num;
            den += dif_v_den;
            x[i] = xn[i];
        }

        fact = sqrt(num)/sqrt(den);

        if (fact < TOL){
            cout << "iteracion: " << it << endl;
            return;
        }

    }

    cout << "Sistema no pudo converger con " << max_it << " iteraciones " << endl;
    return;
}

void Jacobi_TriD(vvd &A, vd &b, vd &x, double TOL, int max_it) {

    int n = b.size();
    vd xn(n, 0);

    for (int i = 0; i < n; i++) {
        x[i] = b[i] / A[i][i];
    }

    for (int it = 0; it < max_it; it++) {

        for (int i = 0; i < n; i++) {

            xn[i] = (b[i] - (i != 0) * A[i][i - 1] * x[i - 1] - (i != n - 1) * A[i][i + 1] * x[i + 1]) / A[i][i];

        }

        double num = 0;
        double den = 0;
        double fact;

        for (int i = 0; i < n; i++) {
            double dif_v_num = (xn[i] - x[i]);
            double dif_v_den = xn[i] * xn[i];
            num += dif_v_num*dif_v_num;
            den += dif_v_den;
            x[i] = xn[i];
        }

        fact = sqrt(num)/sqrt(den);

        if (fact < TOL){
            cout << "iteracion: " << it << endl;
            return;
        }
    }

    return;
}


void Gauss_Seidel(vvd &A, vd &b, vd &x, double TOL, int max_it) {

    int n = b.size();
    vd xn(n, 0);

    // inicializacion
    for (int i = 0; i < n; i++) {
        x[i] = b[i] / A[i][i];
    }

    // calculo de nuevo vector
    for (int it = 0; it < max_it; it++) {
        for (int i = 0; i < n; i++) {
            double sum = 0;
            //valores nuevos
            for (int j = 0; j < i; j++) {
                sum += A[i][j] * xn[j];
            }
            // valores viejos
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }

            xn[i] = (b[i] - sum) / A[i][i];
        }


        double num = 0;
        double den = 0;
        double fact;

        //actualizacion del vectors y calculo de error
        for (int i = 0; i < n; i++) {
            double dif_v_num = (xn[i] - x[i]);
            double dif_v_den = xn[i] * xn[i];
            num += dif_v_num * dif_v_num;
            den += dif_v_den;
            x[i] = xn[i];
        }

        fact = sqrt(num)/sqrt(den);

        if (fact < TOL){
            cout << "iteracion: " << it << endl;
            return;
        }

    }

    cout << "Sistema no pudo converger con " << max_it << " iteraciones " << endl;

}

void Gauss_Seidel_TriD(vvd &A, vd &b, vd &x, double TOL, int max_it) {

    int n = b.size();
    vd xn(n, 0);

    // inicializacion
    for (int i = 0; i < n; i++) {
        x[i] = b[i] / A[i][i];
    }

    // calculo de nuevo vector
    for (int it = 0; it < max_it; it++) {

        for (int i = 0; i < n; i++) {

            xn[i] = (b[i] - (i != 0) * A[i][i - 1] * xn[i - 1] - (i != n - 1) * A[i][i + 1] * x[i + 1]) / A[i][i];

        }

        double num = 0;
        double den = 0;
        double fact;

        //actualizacion del vectors y calculo de error
        for (int i = 0; i < n; i++) {
            double dif_v_num = (xn[i] - x[i]);
            double dif_v_den = xn[i] * xn[i];
            num += dif_v_num * dif_v_num;
            den += dif_v_den;
            x[i] = xn[i];
        }

        fact = sqrt(num)/sqrt(den);

        if (fact < TOL){
            cout << "iteracion: " << it << endl;
            return;
        }

    }

    cout << "Sistema no pudo converger con " << max_it << " iteraciones " << endl;

}

void solveEC_Iterativo(vd &x, double Q, double K, double phi_0, double phi_n, double n, double L ){

    int des;
    double cf = (Q * (L/n)*(L/n))/K;

    vvd M(n, vd(n, 0));
    vvd H(n, vd (n, 0));
    vd b(n, cf);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j)
                M[i][j] = 2;

            if (i == j - 1 || i == j + 1)
                M[i][j] = -1;
        }
    }

    b[0] = cf + phi_0;
    b[n - 1] = cf + phi_n;


    printf("Elija que metodo quiere usar \n");
    printf("1. Gauss_Seidel \n");
    printf("2. Jacobi \n");

    cin >> des;

    if(des == 1){
        Gauss_Seidel(M, b, x, 0.0001, 50000);
    }else{
        Jacobi(M, b, x, 0.0001, 50000);
    }

    return;
}


void metodo_Potencia(vvd A, vvd &eigenvectors, vd &eigenvalues, double error, int range, int max_it){

    int n = A.size();
    double a_k;
    double lambda = 0;
    double lambda_err = 0;
    int total_iterations;
    vd v0;
    eigenvalues.assign(range, 0.0);
    eigenvectors.assign(range, vd(n, 0.0));


    for (int e_i = 0; e_i < range; e_i++) {
        v0.assign(n, 1);

        for (int k = 0; k < e_i; k++) {
            a_k = eigenvectors[k]*v0;
            for (int j = 0; j < n; j++) {
                v0[j] = v0[j] - eigenvectors[k][j] * a_k;
            }
        }

        eigenvectors[e_i] = A*v0;
        normalize(eigenvectors[e_i]);

        for (int it = 0; it < max_it; it++) {

            lambda_err = lambda;

            v0 = A*eigenvectors[e_i];
            lambda = eigenvectors[e_i]*v0;

            for (int k = 0; k < e_i; k++) {
                a_k = eigenvectors[k]*v0;
                for (int j = 0; j < n; j++) {
                    v0[j] = v0[j] - eigenvectors[k][j] * a_k;
                }
            }

            double num = 0;
            for (int i = 0; i < n; i++) {
                eigenvectors[e_i][i] = v0[i];
            }

            normalize(eigenvectors[e_i]);

            if (norma(A*eigenvectors[e_i] - lambda*eigenvectors[e_i]) <= error){
                //cout << "Convergencia en: " << it + 1 << endl;
                break;
            }
        }

        eigenvalues[e_i] = lambda;
    }

    return;
}

void metodo_Potencia_Inv(vvd &A, vvd &eigenvectors, vd &eigenvalues, double tol, int range, int max_it) {

    int n = A.size();

    vvd Ch(n, vd(n, 0.0));
    vd v0(n, 1.0);
    vd Av(n, 0.0);

    double a_k;
    double lambda = 0;
    double lambda_err = 0;
    eigenvalues.assign(range, 0.0);
    eigenvectors.assign(range, vd(n, 0.0));

    Cholesky_PentaD(A, Ch);

    int total_iterations = 0;

    for (int e_i = 0; e_i < range; e_i++) {

        v0.assign(n, 1);

        for (int k = 0; k < e_i; k++) {

            a_k = eigenvectors[k]*v0;

            for (int j = 0; j < n; j++) {
                v0[j] = v0[j] - eigenvectors[k][j] * a_k;
            }
        }

        normalize(v0);

        for (int it = 0; it < max_it; it++) {

            lambda_err = lambda;

            Triangular_Inferior(Ch, v0, eigenvectors[e_i]);
            Triangular_Superior(Ch, eigenvectors[e_i], eigenvectors[e_i]);

            normalize(eigenvectors[e_i]);
            Av = A*eigenvectors[e_i];
            lambda = eigenvectors[e_i]*Av;

            //actualizar iteracion
            for (int i = 0; i < n; i++) {
                v0[i] = eigenvectors[e_i][i];
            }

            if (fabs(lambda - lambda_err) < tol) {
                cout << "Convergencia en: " << it + 1 << endl;
                break;
            }

            for (int k = 0; k < e_i; k++) {
                a_k = eigenvectors[k]*v0;
                for (int j = 0; j < n; j++) {
                    v0[j] = v0[j] - eigenvectors[k][j] * a_k;
                }
            }
        }

        eigenvalues[e_i] = lambda;

    }

    return;
}


int JacobiKrylov(vvd &A, vd &eigenvalues, vvd &G, double tol, int max_it) {

    int n = A.size();
    int i_max = 0, j_max = 0;
    int i_max_prev = 1, j_max_prev = 0;

    vvd A_tmp;

    //inicializar G
    G.assign(n, vd(n, 0.0));
    A_tmp.assign(n, vd(n, 0.0));

    for (int i = 0; i < n; i++)
        G[i][i] = 1;


    for (int it = 0; it < max_it; it++) {

        //find MAX
        double max = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j)
                    continue;
                if (fabs(A[i][j]) > max) {
                    i_max = i;
                    j_max = j;
                    max = fabs(A[i][j]);
                }
            }
        }

        // criterio de paro
        if (max < tol) {
            for (int i = 0; i < n; i++)
                eigenvalues[i] = A[i][i];
            return it;
        }

        G[i_max_prev][i_max_prev] = G[j_max_prev][j_max_prev] = 1;
        G[i_max_prev][j_max_prev] = G[j_max_prev][i_max_prev] = 0;

        i_max_prev = i_max;
        j_max_prev = j_max;

        double delta = (A[j_max][j_max] - A[i_max][i_max]) / (2 * A[i_max][j_max]);
        double sign = delta < 0 ? -1 : delta > 0 ? 1 : 0;
        double t = sign / (fabs(delta) + sqrt(1 + delta * delta));
        double c = 1 / (sqrt(1 + t * t));
        double s = c * t;

        G[i_max][i_max] = G[j_max][j_max] = c;
        G[i_max][j_max] = s;
        G[j_max][i_max] = -s;

        A_tmp = A*G;

        double tmp = G[i_max][j_max];
        G[i_max][j_max] = G[j_max][i_max];
        G[j_max][i_max] = tmp;

        A = G*A_tmp;
    }

    cout << "WARNING JACOBI :Metodo no pudo converger " << endl;
    return max_it;
}

void Jacobi_Eigen_Values(vvd &A, vd &eigenvalues, vvd &eigenvectors, double tol, int max_it) {

    int n = A.size();
    int i_max = 0, j_max = 0, it;
    double max;

    vd G1(n,0), G2(n,0);
    vvd A_tmp(n, vd(n, 0.0));
    vvd evect2(n, vd(n, 0.0));
    bool flag = 0;

    for(int i = 0; i<n; ++i)
        for(int j = 0; j<n; ++j)
            eigenvectors[i][j] = evect2[i][j] = (i==j);

    for (it = 0; it < max_it; it++) {

        max = max_element(A, i_max, j_max);

        if (max < tol) {

            for (int i = 0; i < n; i++)
                eigenvalues[i] = A[i][i];
            cout << "Convergencia en " << it + 1 << " iteraciones" << endl;
            return;
        }

        double theta = atan2(2*A[i_max][j_max], A[j_max][j_max] - A[i_max][i_max])/2;
        double c = cos(theta);
        double s = sin(theta);

        G1[i_max] = G2[j_max] = c;
        G1[j_max] = -s; G2[i_max] = s;

        A_tmp = A;
        for(int i = 0; i<n; ++i)
            A_tmp[i][i_max] = A[i][i_max]*G1[i_max] + A[i][j_max]*G1[j_max];
        for(int i = 0; i<n; ++i)
            A_tmp[i][j_max] = A[i][i_max]*G2[i_max] + A[i][j_max]*G2[j_max];

        A = A_tmp;
        for(int j = 0; j<n; ++j)
            A[i_max][j] = G1[i_max]*A_tmp[i_max][j] + G1[j_max]*A_tmp[j_max][j];
        for(int j = 0; j<n; ++j)
            A[j_max][j] = G2[i_max]*A_tmp[i_max][j] + G2[j_max]*A_tmp[j_max][j];


        for(int i = 0; i < n; ++i)
            evect2[i_max][i] = eigenvectors[i_max][i]*G1[i_max] + eigenvectors[j_max][i]*G1[j_max];
        for(int i = 0; i<n; ++i)
            evect2[j_max][i] = eigenvectors[i_max][i]*G2[i_max] + eigenvectors[j_max][i]*G2[j_max];

        eigenvectors = evect2;

        G1[i_max] = G2[j_max] = 0;
        G1[j_max] = G2[i_max] = 0;
    }

    for (int i = 0; i < n; i++)
        eigenvalues[i] = A[i][i];

    cout << "Numero de iteraciones -> " << it << endl;
    cout << "Metodo no pudo converger " << endl;

    return;
}

void metodoSubespacio(vvd &A, vvd &eigenvectors, vd &eigenvalues, double tol, bool powerMethodType, int max_it) {

    int n = A.size();
    int s = eigenvectors.size();

    vd Av;

    vvd L, U, phiT_A, phi_A, eigenvectors_T, G, J;

    L.assign(n, vd(n, 0.0));
    U.assign(n, vd(n, 0.0));

    phiT_A.assign(s, vd(n, 0.0));
    phi_A.assign(n, vd(s, 0.0));

    G.assign(s, vd(s, 0.0));
    J.assign(s, vd(s, 0.0));

    Identity(eigenvectors);

    if (powerMethodType)
        factorizarLDLT(A, L, U);

    for (int it = 0; it < max_it; it++) {

        for (int ev_pos = 0; ev_pos < s; ev_pos++) {
            grand_schmidt(eigenvectors, eigenvectors[ev_pos], ev_pos);
            if (powerMethodType) {
                Triangular_Inferior(L, eigenvectors[ev_pos], eigenvectors[ev_pos]);
                Triangular_Superior(U, eigenvectors[ev_pos], eigenvectors[ev_pos]);
            } else {
                eigenvectors[ev_pos] = A*eigenvectors[ev_pos];
            }
            normalize(eigenvectors[ev_pos]);
        }

        eigenvectors_T = Transpose(eigenvectors);

        phi_A = A*eigenvectors_T;
        J = eigenvectors*phi_A;

        it += JacobiKrylov(J, eigenvalues, G, tol, max_it) + 1;

        phiT_A = G*eigenvectors;

        double err = 0;
        for (int i = 0; i < s; i++) {
            double v_err = 0;
            //v_err = (G[i][i] - 1)*(G[i][i]-1);
            for (int j = 0; j < n; j++) {
                eigenvectors[i][j] = phiT_A[i][j];
                double dif = eigenvalues[i] * eigenvectors[i][j] - phi_A[j][i];
                v_err += dif * dif;
            }

            if (sqrt(v_err) > err)
                err = sqrt(v_err);
        }

        if (err < tol) {
            cout << " Convergencia en iteracion: " << it << endl;
            return;
        }
    }

    cout << "WARNING METODO SUBESPACIO: metodo no pudo converger" << endl;

}

void metodo_Subespacio_Pot(vvd &A, vvd &eigenvectors, vd &eigenvalues, double tol, int max_it, int m) {

    int n = A.size();
    vvd phi_0(m, vd(n, 0));

    for(int i = 0; i<m; ++i){
        phi_0[i][i] = 1;
    }

    bool conv;
    for(int it = 0; it < max_it; ++it){

        for (int j = 0; j < m; j++) {
            grand_schmidt(phi_0, phi_0[j], j);
            phi_0[j] = phi_0[j]*A;
            normalize(phi_0[j]);
        }

        vvd B = phi_0*(A*Transpose(phi_0));

        vvd evect(m, vd(m, 0));
        vd eval(m, 0);

        Jacobi_Eigen_Values(B, eval, evect, tol, max_it);
        phi_0 = evect*phi_0;

        for(int j = 0; j < m; ++j){
            normalize(phi_0[j]);
        }

    }

    eigenvectors = phi_0;

    phi_0 = phi_0*(A*Transpose(phi_0));

    for(int i = 0; i < m; i++){
        eigenvalues[i] = phi_0[i][i];
    }

    return;
}


void metodo_Subespacio_Pot_Inv(vvd &A, vvd &eigenvectors, vd &eigenvalues, double tol, int max_it, int m) {

    int n = A.size();
    vvd phi_0(m, vd(n, 0));

    vvd Ch(n, vd(n));

    Cholesky_PentaD(A, Ch);

    for(int i = 0; i < m; ++i){
        phi_0[i][i] = 1;
    }

    bool conv;
    for(int it = 0; it < max_it; ++it){

        for (int j = 0; j < m; j++) {

            Triangular_Inferior(Ch, phi_0[j], phi_0[j]);
            Triangular_Superior(Ch, phi_0[j], phi_0[j]);
            normalize(phi_0[j]);
            grand_schmidt(phi_0, phi_0[j], j);
            normalize(phi_0[j]);
        }

        vvd B = phi_0*(A*Transpose(phi_0));

        conv = true;
        for(int i = 0; i < m; i++){
            for(int j = 0; j < m; j++){
                if(i != j && fabs(B[i][j]) > tol){
                    conv = false;
                }
            }
        }

        if(conv){
            cout << "\nConvergencia en Iteracion: " << it + 1 << "\n" << endl;
            break;
        }

        vvd evect(m, vd(m, 0));
        vd eval(m, 0);

        Jacobi_Eigen_Values(B, eval, evect, tol, max_it);
        phi_0 = evect*phi_0;

        for(int j = 0; j < m; ++j){
            normalize(phi_0[j]);
        }

    }

    if(!conv){
        cout << "\nWARNING: Metodo no pudo converger \n" << endl;
    }

    eigenvectors = phi_0;

    phi_0 = phi_0*(A*Transpose(phi_0));

    for(int j = 0; j < m; ++j){
        eigenvalues[j] = phi_0[j][j];
    }


    return;
}



void cociente_Rayleigh(vvd &A, vd &v, double &lambda, double tol, int max_it) {

    int n = A.size();

    auto clock_start = high_resolution_clock::now();

    //inicializar
    double rho, lambda2;
    vvd Q, R, A_rho;
    vd Av, v_temp, v_b;

    Q.assign(n, vd(n, 0.0));
    R.assign(n, vd(n, 0.0));
    A_rho.assign(n, vd(n, 0.0));
    Av.assign(n, 0.0);
    v_b.assign(n, 0.0);

    v_temp = v;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A_rho[i][j] = A[i][j];

    normalize(v);

    Av = A*v;
    rho = v*Av;

    for (int it = 0; it < max_it; it++) {

        for (int i = 0; i < n; i++)
            A_rho[i][i] = A[i][i] - rho;


        Factorizacion_QR(A_rho, Q, R);

        v_b = Transpose(Q)*v;
        Triangular_Superior(R, v_b, v);

        normalize(v);

        Av = A*v;
        rho = v*Av;

        lambda = rho;
        if(it >= 1){
            lambda2 = lambda;
        }else{
            lambda2 = 0;
        }

        double err = 0;
        for (int i = 0; i < n; i++) {
            double dif = Av[i] - lambda*v[i];
            err += dif*dif;
        }
        if (sqrt(err) < tol){
            cout << "Convergencia en Iteracion: " << it + 1 << endl;
            cout << "Error del valor propio: " << abs(lambda2 - lambda) << endl;
            cout << "Error del vector propio: " << norma(v_temp-v) << endl;

            auto clock_stop = high_resolution_clock::now();
            auto duracion = duration_cast<microseconds>(clock_stop - clock_start);
            cout << "\n El tiempo de ejecucion del programa es: " << duracion.count()*1e-6 << "\n" << endl;

            return;
        }

        v_temp = v;

    }
    cout << "WARNING RAYLEIGH: Metodo no pudo converger " << endl;

    return;
}

void Gradiente_Conjudado(vvd &A, vd &b, vd &x, vvd &M, bool flag, double tol, int max_it){

    int n = A.size();
    double alpha, beta;

    if(flag){
        M = Inverse(M);
    }

    vd r_0, r_1, p_0, z_0, z_1, Ap, alph_p, x2;
    x.assign(n, 0.0);
    p_0.assign(n, 0.0);
    z_0.assign(n, 0.0);

    r_0 = b - A*x;
    z_0 = M*r_0;
    p_0 = z_0;

    x2 = x;

    for (int it = 0; it < max_it; it++){

        alpha = (r_0*z_0)/((p_0*A)*p_0);

        x = x + alpha*p_0;

        r_1 = r_0;
        r_0 = r_0 - alpha*(A*p_0);

        if (norma(r_0) < tol) {
            cout << "Convergencia en Iteracion: " << it + 1 << endl;
            cout << "Error del vector x" << norma(x2 - x) << endl;
            return;
        }

        z_1 = z_0;
        z_0 = M*r_0;

        beta  = (z_0*r_0)/(z_1*r_1);
        p_0 = z_0 + beta*p_0;

        x2 = x;
    }

    cout << "Metodo no pudo converger "<< endl;

    return;
}

void Minimos_Cuadrados(vd x, vd y, int m, vd &coef){

    m += 1;
    coef.assign(m, 0);
    int n = y.size();
    vvd X(n, vd(m, 0));

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            X[i][j] = pow(x[i], j);
        }
    }

    vd c = Inverse(Transpose(X)*X)*Transpose(X)*y;

    for(int i = 0; i < m; i++){
        coef[i] = c[i];
    }

    return;
}

void Minimos_Cuadrados_Generalizado(vvd X, vd y, vd &coef){

    int m = X[0].size();
    coef.assign(m, 0);

    vd c = Inverse(Transpose(X)*X)*Transpose(X)*y;

    for(int i = 0; i < m; i++){
        coef[i] = c[i];
    }

    return;
}

void Matricial_Interpolation(vd x, vd y, vd &coef){

    int n = y.size();
    coef.assign(n, 0);
    vvd X(n, vd(n, 0));

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            X[i][j] = pow(x[i], j);
        }
    }

    Eliminacion_Gaussiana(X, y, coef);

    return;
}

poly Lagrange_Interpolation(vd x, vd y){

    int n = x.size();
    poly pz(1,{0}), z(2, {0,1});
    for(int i = 0; i < n; ++i){
        poly l(1,{y[i]});
        for(int j = 0; j < n; ++j)
            if(i != j)
                l = l*((z - x[j])/(x[i] - x[j]));
       pz = pz + l;
    }

    return pz;
}

long double Lagrange_Interpolation(vld x, vld y, long double z){

    int n = x.size();
    long double px = 0, l;

    for(int i = 0; i < n; i++){
        l = y[i];
        for(int j = 0; j < n; j++){
            if(i != j){
                l *= (z - x[j])/(x[i] - x[j]);
            }
        }
       px += l;
    }

    return px;
}

long double Newton_Interpolation(vld x, vld y, long double z){

    int i, j;
    long double px = 0, mult;

    int n = y.size();

    for(j = 0; j < n - 1; j++) {
        for(i = n - 1; i > j; i--) {
            y[i] = (y[i] - y[i - 1]) / (x[i] - x[i - j-1]);
        }
    }
    for(i = n - 1; i >= 0; i--) {
        mult = 1;
        for(j = 0; j < i; j++){
            mult *= (z - x[j]);
        }

        mult *= y[j];
        px += mult;
    }

    return px;
}

void lineal_Splines(vd x, vd y, vvd &coef){

    int n = x.size();

    for(int i = 0; i < n - 1; ++i){
        coef[0][i] = (y[i + 1] - y[i])/(x[i + 1] - x[i]);
        coef[1][i] = y[i];
    }

    return;
}

void quadratic_Splines(vd x, vd y, vvd &coef, double diff){

    int n = x.size();

    vd z(n, 0);
    z[0] = diff;

    for(int i = 1; i < n; ++i){
        z[i] = 2*(y[i] - y[i - 1])/(x[i] - x[i - 1]) - z[i - 1];

        coef[0][i - 1] = 0.5*(z[i] - z[i - 1])/(x[i] - x[i - 1]);
        coef[1][i - 1] = z[i - 1];
        coef[2][i - 1] = y[i - 1];
    }

    return;
}


void cubic_Splines(vd x, vd y, vvd &coef){

    int n = x.size();

    vd h(n - 1, 0);
    vd b_temp(n - 1, 0);
    vd b(n - 2, 0);
    vd z(n - 2, 0);
    vvd A(n - 2, vd(n - 2, 0));
    vvd H(n - 2, vd(n - 2, 0));

    for(int i = 1; i < n; ++i){
        h[i - 1] = x[i] - x[i - 1];
        b_temp[i - 1] = (y[i] - y[i - 1])/h[i - 1];
    }

    for(int i = 0; i < n - 2; ++i){

        b[i] = 6*(b_temp[i + 1] - b_temp[i]);

        if(i == 0){
            A[i][i] = 2*(h[0] + h[1]);
            A[i][i + 1] = h[1];
        }else if(i == n - 3){
            A[i][i - 1] = h[i];
            A[i][i] = 2*(h[i] + h[i + 1]);
        }else{
            A[i][i - 1] = h[i];
            A[i][i] = 2*(h[i] + h[i + 1]);
            A[i][i + 1] = h[i + 1];
        }
    }

    Cholesky_TriD(A, H);

    Triangular_Inferior(H, b, z);
    Triangular_Superior(H, z, z);

    z.insert(z.begin(), 0);
    z.insert(z.end(), 0);

    for(int i = 0; i < n - 1; ++i){
        coef[0][i] = z[i + 1]/(6*h[i]);
        coef[1][i] = z[i]/(6*h[i]);
        coef[2][i] = y[i + 1]/h[i] - h[i]*z[i + 1]/6;
        coef[3][i] = y[i]/h[i] - h[i]*z[i]/6;
    }

    return;
}


void FEM_Interpolation(vd x, vd y, vd z, vd &coef, double lambda){

    int n = x.size(), m = z.size();

    vvd A(m, vd(m,0));
    vd b(m,0);

    for(int i = 0, j = 0; i < m - 1; ++i)
    {
        double h = z[i + 1] - z[i];
        A[i][i] += lambda/h;
        A[i][i+1] -= lambda/h;
        A[i+1][i] -= lambda/h;
        A[i+1][i+1] += lambda/h;
        for(; j < n && x[j] < z[i + 1]; ++j){

            double N = 1 - (x[j] - z[i])/h;
            double M = (x[j] - z[i])/h;

            A[i][i] += N*N;
            A[i][i+1] += N*M;
            A[i+1][i] += N*M;
            A[i+1][i+1] += M*M;
            b[i] += y[j]*N;
            b[i+1] += y[j]*M;
        }
    }

    coef.assign(m, 0);

    vvd H(m, vd(m, 0));

    Cholesky_TriD(A, H);
    Triangular_Inferior(H, b, coef);
    Triangular_Superior(H, coef, coef);

}

double mcmc_integration_2d(double a, double b, double (*f)(double)){

    int n = 1e6;
    double h = (b - a)/n;
    double maximo = 0, minimo = 1e18;

    for(double x = a; x < b; x += h){
        maximo = max(maximo, f(x)), minimo = min(minimo, f(x));
    }

    maximo = 1.3*maximo, minimo = 1.3*minimo;

    double sum = 0;
    int accp = 0, accn = 0, totp = 0, totn = 0;

    srand(time(NULL));
    mt19937 gen(rand());
    mt19937 genx(rand());
    uniform_real_distribution<> dist_y(minimo, maximo), dist_x(a, b);

    int global = 0;
    while(global < 100000){

        ++global;
        for(int i = 0; i < 5000; ++i){
            double rng = dist_y(gen), x_i = dist_x(genx);
            if(f(x_i) >= 0 && f(x_i) > rng)++accp;
            if(f(x_i) < 0 && f(x_i) < rng)++accn;
            if(f(x_i) >= 0 )++totp;
            if(f(x_i) < 0 )++totn;
        }

        double tmp1 = totp ? double(accp)/totp : 0, tmp2 = totn ? double(accn)/totn : 0;
        double tmp = (tmp1 - tmp2)*(b - a)*(maximo - minimo);
        if(fabs(tmp - sum) < 0.0000001)return tmp;
        sum = tmp;

    }

    return sum;
}

double mcmc_integration_3d(double a, double b, double c, double d, double (*f)(double, double)){

    int n = 1e4;
    double h1 = (b - a)/n, h2 = (d - c)/n;
    double maximo = 0, minimo = 1e18;

    for(double x = a; x < b; x += h1){
        for(double y = c; y < d; y += h2){
            maximo = max(maximo, f(x,y)), minimo = min(minimo, f(x,y));
        }
    }

    maximo = 1.5*maximo, minimo = 1.5*minimo;

    double sum = 0;
    int accp = 0, accn = 0, totp = 0, totn = 0;
    srand(time(NULL));
    mt19937 gen(rand());
    mt19937 genx(rand());
    mt19937 geny(rand());
    uniform_real_distribution<> dis(minimo, maximo), plx(a,b), ply(c,d);

    int global = 0;
    while(global < 100000){

        ++global;
        for(int i = 0; i < 5000; ++i){
            double rng = dis(gen), x_i = plx(genx), y_i = ply(geny);
            if(f(x_i,y_i) >= 0 && f(x_i,y_i) > rng)++accp;
            if(f(x_i,y_i) < 0 && f(x_i, y_i) < rng)++accn;
            if(f(x_i, y_i) >= 0 )++totp;
            if(f(x_i, y_i) < 0 )++totn;
        }

        double tmp1 = totp ? double(accp)/totp : 0, tmp2 = totn ? double(accn)/totn : 0;
        double tmp = (tmp1-tmp2)*(b-a)*(d-c)*(maximo - minimo);
        if(fabs(tmp-sum) < 0.000001)return tmp;
        sum = tmp;
    }

    return sum;

}

double Newton_Cotes(double a, double b, int n, double (*f)(double)){

    double h = (b - a)/n;
    vd x(n, 0);
    for(int i = 0; i < n; ++i){
        x[i] = a + i*h;
    }

    double I = 0;
    for(int i = 1; i < n - 1; i++){
        I += f(x[i - 1]) + 4*f(x[i]) + f(x[i + 1]);
    }

    return (h/6)*I;
}

double Extra_Richardson(double a, double b, int n, double (*f)(double)){

    int n1 = n, n2 = 2*n;

    return (16*Newton_Cotes(a, b, n2, f) - Newton_Cotes(a, b, n1, f))/15;
}


double Romberg(double a, double b, double eps, double (*f)(double)){

    double h;
    vvd R(1, vd(1, 0.5*(b - a)*(f(a) + f(b))));

    int k = R[0].size();
    for(int i = 0; i < k; i++){
        printf(" %lf ", R[0][i]);
    }
    cout << endl;

    int n = 1;
    while(true){
        h = (b - a)/pow(2, n);

        double acum = 0;
        for(int k = 1; k < pow(2, (n - 1)) + 1; ++k){
            acum += f(a + (2*k - 1)*h);
        }

        R.push_back(vd(n + 1, 0));
        R[n][0] = 0.5*R[n - 1][0] + h*acum;
        for (int m = 1; m < n + 1; ++m){
            R[n][m] = R[n][m - 1] + (R[n][m - 1] - R[n - 1][m - 1]) / (pow(4, m) - 1);
        }

        k = R[n].size();

        for(int i = 0; i < k; i++){
            printf(" %lf ", R[n][i]);
        }
        cout << endl;

        if (abs(R[n][n - 1] - R[n][n]) < eps){
            return R[n][n];
        }

        n++;
    }

}

double Cuadratura_Gaussiana(double a, double b, int nodos, double (*f)(double)){

    vd x_i(3, 0);
    vd w_i(3, 0);
    vd c(2, 0);
    vd x(nodos, 0);
    vd y(nodos, 0);
    vvd coef(4, vd(nodos - 1, 0));
    int st, ed;
    double add, z;
    double h = (b - a)/nodos;

    x_i[0] = 0;
    x_i[1] = 0.774597, x_i[2] = -0.774597;
    w_i[0] = 0.888889;
    w_i[1] = w_i[2] = 0.555556;

    c[0] = (b - a)/2;
    c[1] = (b + a)/2;

    for(int i = 0; i < nodos; ++i){
        x[i] = a + i*h;
        y[i] = f(x[i]);
    }

    cubic_Splines;

    double I = 0;
    for(int i = 0; i < 3; ++i){

        st = 1, ed = x.size() - 1;
        z = c[0]*x_i[i] + c[1];

        while(st <= ed){
            int md = (st + ed)/2;
            if(x[md] > z){
                ed = md - 1;
            }else st = md + 1;
        }

        add = coef[0][ed]*(z - x[ed])*(z - x[ed])*(z - x[ed]);
        add -= coef[1][ed]*(z - x[ed + 1])*(z - x[ed + 1])*(z - x[ed + 1]);
        add += coef[2][ed]*(z - x[ed]) - coef[3][ed]*(z - x[ed + 1]);

        I += w_i[i]*add;
    }

    return ((b - a)/2)*I;
}




///Diferenciacion Numerica///

double Diferencias_Finitas(double x , double (*f)(double), int n, double h){

    double df = 0, conv;

    for(int i = 0; i <= n; i++){
        conv = factorial(n)/(factorial(n - i)*factorial(i));
        if(i%2 == 0)df += conv*f(x + (n/2 - i)*h);
        else df -= conv*f(x + (n/2 - i)*h);
    }

    return df/pow(h, n);
}

double Lagrange_Derivative(vd x, vd y, double z, int m){

    poly L = Lagrange_Interpolation(x, y);
    vd coef = L.p;
    int n = coef.size();

    for(int k = 0; k < m; k++){
        for(int i = 0; i < n; i++){
            coef[i] = coef[i]*i;
        }
        coef.erase(coef.begin());
        n = n - 1;
    }

    double px = 0;
    for(int i = 0; i < n; i++){
        px += coef[i]*pow(z, i);
    }

    return px;
}

double Extra_Richardson_diff(double x, double (*f)(double), int n, double h){

    double N1 = Diferencias_Finitas(x, f, n, 2*h);
    double N2 = Diferencias_Finitas(x, f, n, 4*h);
    return N1 + (N1 - N2)/3;
}




    ///Proyecto Final///

void PCA(vvd X, vvd &Y, int n_comp, double &contribution){

    int n = X.size();
    int m = X[0].size();
    vvd S(m, vd(m, 0));
    double mean;

    for(int i = 0; i < m; i++){
        mean = 0;
        for(int j = 0; j < n; j++){
            mean += X[j][i];
        }
        mean = mean/n;
        X[i] = X[i] - mean;
    }

    for(int j = 0; j < m; j++){
        for(int k = 0; k < m; k++){
            for(int i = 0; i < n; i++){
                S[j][k] += X[i][j]*X[i][k];
            }
            S[j][k] = S[j][k]/(n - 1);
        }
    }

    vvd eigenvectors;
    vd eigenvalues;

    metodo_Potencia(S, eigenvectors, eigenvalues, 0.000001, n_comp, 100000);

    Y = X*Transpose(eigenvectors);

    double traza = 0;
    for(int i = 0; i < m; i++){
        traza += S[i][i];
    }

    for(int i = 0; i < n_comp; i++){
        eigenvalues[i] = eigenvalues[i]/traza;
    }

    contribution = sum(eigenvalues)*100;

    return;
}

double Softmax(vd x, int k, vvd beta){

    int n = x.size();
    int k_par = beta.size();

    double num = exp(beta[k]*x);
    double den = 0;
    for(int j = 0; j < k_par; j++){
        den += exp(beta[j]*x);
    }

    return num/den;
}

void Softmax_Estimation(vvd X, vd y, vvd &beta, double alpha, int max_it){

    vvd beta_temp;
    int m = X.size();
    int n = X[0].size();
    int k_par = beta.size();
    int n_par = beta[0].size();

    X = Transpose(X);
    X.insert(X.begin(), vd(m, 1));
    X = Transpose(X);

    vvd gradient(k_par, vd(n_par, 0));
    int ind;
    double num, den;

    for(int it = 0; it < max_it; it++){

        for(int g = 0; g < k_par; g++){
            for(int h = 0; h < n_par; h++){
                gradient[g][h] = 0;
                for(int i = 0; i < m; i++){
                    num = (beta[g]*X[i]);
                    den = 0;
                    for(int j = 0; j < k_par; j++){
                        den += exp(beta[j]*X[i]);
                    }
                    if(y[i] = g)ind = 1;
                    else ind = 0;

                    gradient[g][h] += (num/den - ind)*X[i][h];
                }
            }
        }

        cout << "\n" << endl;

        for(int g = 0; g < k_par; g++){
            for(int h = 0; h < n_par; h++){
                beta[g][h] -= alpha*gradient[g][h];
            }
        }

    }

    return;
}

void Softmax_Predictions(vvd new_x, vvd beta, vd &pred){

    int n = new_x.size();
    int m = new_x[0].size();
    int k_class = 0;

    pred.assign(n, 0);
    double pred_temp;
    int max_prob = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < k_class; j++){
            pred_temp = Softmax(new_x[i], j, beta);
            if(pred_temp > max_prob){
                max_prob = pred_temp;
                pred[i] = j;
            }
        }
    }

    return;
}

void One_Hot(vd y, vvd &Y, int classes){

    int n = y.size();
    Y.assign(n, vd(classes, 0));

    for(int i = 0; i < classes; i++){
        for(int j = 0; j < n; j++){
            if(y[j] == i){
                Y[i][j] = 1;
            }
        }
    }

    return;
}

double Accuracy(vd pred, vd real){

    int n = pred.size();

    double TP = 0, TN = 0;
    double FP = 0, FN = 0;

    for(int i = 0; i < n; i++){
        if(pred[i] == 1 && real[i] == 1) TP++;
        else if(pred[i] == 0 && real[i] == 0) TN++;
        else if(pred[i] == 0 && real[i] == 1) FN++;
        else if(pred[i] == 1 && real[i] == 0) FP++;
    }

    double acc = (TP + TN)/(TP + TN + FN + FP);

    return acc;
}

double Balanced_Accuracy(vd pred, vd real){

    int n = pred.size();

    double TP = 0, TN = 0;
    double FP = 0, FN = 0;

    for(int i = 0; i < n; i++){
        if(pred[i] == 1 && real[i] == 1) TP++;
        else if(pred[i] == 0 && real[i] == 0) TN++;
        else if(pred[i] == 0 && real[i] == 1) FN++;
        else if(pred[i] == 1 && real[i] == 0) FP++;
    }

    double Sensitivity = TP/(TP + FP);
    double Specificity = TN/(TN + FN);

    double b_acc = (Sensitivity + Specificity)/2;

    return b_acc;
}




