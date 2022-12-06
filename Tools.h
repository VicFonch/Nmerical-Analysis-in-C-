#include "typedef.h"

using namespace std;

void Read_Matrix(vvd &M, string filename);
void Read_Vector(vd &V, string filename);
void Write_Matrix(const vvd &M, string filename, bool flag_append);
void Write_Vector(const vd &V, string filename, bool flag_append);
void Write_Pairs(vd &V1, vd &V2, string filename);
void Print_Matrix(vvd M);
void Print_Vector(vd V);
void Scan_Consol_Vector(vd &V);
void Scan_Consol_Matrix(vvd &M);



void Read_Matrix(vvd &M, string filename){

    int n = M.size();
    FILE *in_M = fopen(filename.c_str(), "r");

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            fscanf(in_M, "%lf", &M[i][j]);
        }
    }

    fclose(in_M);

    return;
}

void Read_Vector(vd &V, string filename){

    int n = V.size();
    FILE *in_b = fopen(filename.c_str(), "r");

    for(int i = 0; i < n; i++){
        fscanf(in_b, "%lf", &V[i]);
    }

    fclose(in_b);

    return;
}

void Write_Matrix(const vvd &M, string filename){

    int n = M.size();
    int m = M[0].size();
    FILE *out = fopen(filename.c_str(), "w");

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            fprintf(out, " %lf ", M[i][j]);
        }
        fprintf(out, "\n");
    }

    fclose(out);

    return;
}

void Write_Pairs(vd &V1, vd &V2, string filename){

    int n = V1.size();
    freopen(filename.c_str(), "w", stdout);

    for(int i = 0; i < n; i++){
        cout << V1[i] << " " << V2[i] << '\n';
    }

    fclose(stdout);

    return;
}

void Write_Vector(const vd &V, string filename){

    int n = V.size();
    FILE *out = fopen(filename.c_str(), "w");

    for(int i = 0; i < n; i++){
        fprintf(out, " %lf\n", V[i]);
    }

    fclose(out);

    return;
}


void Print_Matrix(vvd M){

    for(int i = 0; i < M.size(); i++){
        for(int j = 0; j < M[0].size(); j++){
            printf(" %e ", M[i][j]);
        }
        printf("\n");
    }

    return;
}


void Print_Vector(vd V){

    int n = V.size();

    for(int i = 0; i < n; i++){
        printf(" %f ", V[i]);
    }

    return;
}

void Scan_Consol_Matrix(vvd &M){

    for(int i = 0; i < M.size(); i++){

        cout << "Escriba la fila " << i << "\n" << endl;

        for(int j = 0; j < M[0].size(); j++){
            cin >> M[i][j];
        }

    }

    return;
}


void Scan_Consol_Vector(vd &V){

    int n = V.size();

    printf("Introduzca los elementos del vector b\n");

    for(int i = 0; i < n; i++){
        cin >> V[i];
    }

    return;
}





