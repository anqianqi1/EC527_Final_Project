#include <iostream>
#include <stdio.h>
#include <vector>
#include<cmath>
using namespace std;

#define strength_Matrix_row 5
#define strength_Matrix_col 5

#define smooth_Matrix_row 4
#define smooth_Matrix_col 4
#define ITERATION_LIMIT 1000

typedef struct{
    int value;

    //if the vertex has been handled
    bool fresh;
    bool cpoint;
    bool fpoint;
} vertex;

class coarse_matrix{
public:
    int row;
    int col;

    //data for coarse
    vertex **v;//a matrix of v(vertex)


    //result
    vector<vector<int> > cpoint_coor;


    coarse_matrix(int ROW, int COL){
        //this->matrix = matrix;
        row = ROW;
        col = COL;

        //vertex v[row][col];

        v = new vertex*[row];
        for( int i = 0; i < row; i++) {
            vertex* j = new vertex[col];
            v[i] = j;
        }


//        if(!(v = (vertex**)malloc(row * sizeof(vertex*)))){
//            printf("error on malloc matrix");
//        }
//        for(int c = 0; c < col; c++){
//            if(!(v[c] = (vertex*)malloc(sizeof(vertex)))){
//                printf("error on malloc matirx");
//            }
//        }

        //initial arribute
        for(int r = 0; r < row; r++){
            for(int c = 0; c < col; c++){
                v[r][c].fresh = true;
                v[r][c].cpoint = false;
                v[r][c].fpoint = false;
                v[r][c].value = -100000;//-1 means value is invalid
            }
        }
        print_value();
    }

    void input_data(int input[][strength_Matrix_col]){
        for(int i = 0; i < row; i++){
            for(int j = 0; j < col; j++){
                //printf("!%d\n",data[i][j]);
                v[i][j].value = input[i][j];
                //printf("(%d,%d)#%d\n",i,j, v[i][j].value);
            }
        }
        print_value();
    }


    void input_RHS_data(int input[][smooth_Matrix_col]){
        for(int i = 0; i < row; i++){
            for(int j = 0; j < col; j++){
                //printf("!%d\n",data[i][j]);
                v[i][j].value = input[i][j];
                //printf("(%d,%d)#%d\n",i,j, v[i][j].value);
            }
        }
        print_value();
    }

    void print_value(){
        printf("\nMatrix Value:\n");
        for(int r = 0; r < row; r++){
            for(int c = 0; c < col; c++){
                printf("%d  ",v[r][c].value);
            }
            printf("\n");
        }
        printf("\n");
    }

    void print_stats(){
        printf("\nMatrix Cpoint:\n");
        for(int r = 0; r < row; r++){
            for(int c = 0; c < col; c++){
                if(v[r][c].cpoint == true){
                    printf("c ");
                }else if(v[r][c].fpoint == true){
                    printf("f ");
                }else if(v[r][c].fresh == true) {
                    printf("%d ",v[r][c].value);
                }else{
                    printf("e");
                }
            }
            printf("\n");
        }
        printf("\n");
    }

    void coarse(){
        //result data

        //function variables
        int max_x, max_y;
        int max_value = -10000000;
        bool done = 0;

        //start coarse
        while(!done){
            //printf("*****************************\n");

            //find max value of matrix(find cpoint)
            for(int i = 0; i < row; i++){
                for(int j = 0; j < col; j++){
                    if(v[i][j].value > max_value && v[i][j].fresh == true){
                        max_value = v[i][j].value;
                        max_x = i;
                        max_y = j;
                    }
                }
            }
            //store coordinate of cpoint
            vector<int> tmp = {max_x,max_y};
            cpoint_coor.push_back(tmp);

            //print cpoint coordinates
//            printf("cpoint coordinates:\n");
//            for(int i = 0; i < cpoint_coor.size(); i++){
//                printf("(%d,%d) ",cpoint_coor[i][0],cpoint_coor[i][1]);
//            }
//            printf("\n");

//            printf("max_value(%d,%d) = %d\n",max_x,max_y,max_value);

            //////////Whether update cpoint and fpoint value?
            //        for(int i = -1; i <= 1; i++){
            //            for(int j = -1; j <= 1; j++){
            //                if(i == 0 && j == 0){//self vertex which is cpoint
            //                    v[max_x+i][max_y+j].fresh = false;
            //                    v[max_x+i][max_y+j].cpoint = true;
            //                }else{//iterate around vertex which are fpoints
            //                    v[max_x+i][max_y+j].fresh = false;
            //                    v[max_x+i][max_y+j].fpoint = true;
            //                    /////////////
            //                    printf("Before adding value:");
            //                    print_stats();
            //                    /////////////
            //                    for(int x = -1; x <=1; x++){
            //                        for(int y = -1; y <=1; y++){
            //                            if((x == 0 && y == 0) || (max_x+i+x < 0) || (max_x+i+x >= col) || (max_y+j+y < 0) || (max_y+j+y >= row)){
            //                                continue;
            //                            }else if(v[max_x+i+x][max_y+j+y].fresh == true){
            //                                v[max_x+i+x][max_y+j+y].value += 1;
            //                            }else{
            //                                continue;
            //                            }
            //                        }
            //                    }
            //                }
            //            }
            //        }
            ////////whether update cpoint and fpoint value?

            //iterate around cpoint if they are "fresh", make them as fpoint
            vector<vector<int> > new_fpoint;
            for(int i = -1; i <= 1; i++){
                for(int j = -1; j <= 1; j++){
                    if(i == 0 && j == 0 && v[max_x+i][max_y+j].fresh == true){//self vertex which is cpoint
                        v[max_x+i][max_y+j].fresh = false;
                        v[max_x+i][max_y+j].cpoint = true;
                    }else if(v[max_x+i][max_y+j].fresh == true){//iterate around vertex which are fpoints
                        v[max_x+i][max_y+j].fresh = false;
                        v[max_x+i][max_y+j].fpoint = true;

                        vector<int>v = {max_x+i,max_y+j};
                        new_fpoint.push_back(v);

                        /////////////
                        //printf("Before adding value:");
                        //print_stats();
                        /////////////
                    }else{
//                        printf("!!!Error in setting fpoint.\n");
//                        printf("fpoint? %d\n",v[max_x+i][max_y+j].fpoint);
//                        printf("cpoint? %d\n",v[max_x+i][max_y+j].cpoint);
//                        printf("frsh? %d\n",v[max_x+i][max_y+j].fresh);
                        continue;
                    }
                }
            }

            //increase fresh vertex's value around new fpoint
            for(int i = 0; i < new_fpoint.size(); i++){
                int new_fpoint_x = new_fpoint[i][0];
                int new_fpoint_y = new_fpoint[i][1];
                for(int x = -1; x <= 1; x++){
                    for(int y = -1; y <= 1; y++){
                        if ((x == 0 && y == 0) || (new_fpoint_x + x < 0) || (new_fpoint_x + x >= col) ||
                            (new_fpoint_y + y < 0) || (new_fpoint_y + y >= row)) {
                            continue;
                        }else if(v[new_fpoint_x + x][new_fpoint_y + y].fresh == true){
                            v[new_fpoint_x + x][new_fpoint_y + y].value += 1;
                        }
                        else{
//                            printf("!!!Error at incrementing value of surrounding fresh point.\n");
//                            printf("fpoint? %d\n",v[new_fpoint_x + x][new_fpoint_y + y].fpoint);
//                            printf("cpoint? %d\n",v[new_fpoint_x + x][new_fpoint_y + y].cpoint);
//                            printf("fresh? %d\n",v[new_fpoint_x + x][new_fpoint_y + y].fresh);
                            continue;
                        }
                    }
                }
            }


            //error version
//            for(int i = -1; i <= 1; i++) {
//                for (int j = -1; j <= 1; j++) {
//                    for (int x = -1; x <= 1; x++) {
//                        for (int y = -1; y <= 1; y++) {
//                            if ((x == 0 && y == 0) || (max_x + i + x < 0) || (max_x + i + x >= col) ||
//                                (max_y + j + y < 0) || (max_y + j + y >= row)) {
//                                continue;
//                            } else if (v[max_x + i + x][max_y + j + y].fresh == true) {
//                                v[max_x + i + x][max_y + j + y].value += 1;
//                            } else {
//                                continue;
//                            }
//                        }
//                    }
//                }
//            }

            //check whether all the elements has been handled, then set "done" = 1
            for(int i = 0; i < row; i ++){
                for(int j = 0; j < col; j++){
                    if(v[i][j].fresh == true) {
                        done = 0;
                        break;
                    }else{
                        done = 1;
                    }
                }
            }
            //printf("Done? %d\n",done);
            //printf("After coarse:\n");
            //print_stats();

            //clear data
            max_value = -100000;
            new_fpoint.clear();
        }
    }
};

class smooth_matrix{
public:
    int row;
    int col;

    double **A;
    double *b;

    smooth_matrix(int ROW,int COL){
        row = ROW;
        col = COL;

        A =  new double*[row];
        for(int i = 0; i < row; i++){
            A[i] = new double[col];
        }

        b = new double[row];
    }

    void input_data(double A_matrix[][smooth_Matrix_col], double b_matrix[smooth_Matrix_row]){
        for(int i = 0; i < row; i++){
            for(int j = 0; j < col; j++){
                A[i][j] = A_matrix[i][j];
            }
        }
        for(int i = 0; i < row; i++){
            b[i] = b_matrix[i];
        }
    }

    void print_value(){
        printf("\nA Matrix Value:\n");
        for(int r = 0; r < row; r++){
            for(int c = 0; c < col; c++){
                printf("%f  ",A[r][c]);
            }
            printf("\n");
        }
        printf("\nb Vector Value:\n");
        for(int r = 0; r < row; r++){
            printf("%f  ",b[r]);
        }
        printf("\n");
    }

    void print_b(int B[]){
        printf("[");
        for(int r = 0; r < row; r++){
            if(r == row - 1) printf("%d",B[r]);
            else printf("%d  ",B[r]);
        }
        printf("]\n");
    }

    void print_b(double B[]){
        printf("[");
        for(int r = 0; r < row; r++){
            if(r == row - 1) printf("%f",B[r]);
            else printf("%f  ",B[r]);
        }
        printf("]\n");
    }

    double oneD_dot(int ROW, double x[], int start_index, int end_index){
        int result;
        for(int i = start_index; i < end_index; i++){
            result =+ A[ROW][i] * x[i];
        }
        return result;
    }

    double* twoD_dot(double x[],double result[]){
        //double result[row];
        for(int i = 0; i < row; i++){
            double tmp = 0;
            for(int j = 0; j < col; j++){
                tmp += A[i][j] * x[i];
            }
            result[i] = tmp;
        }
        return result;
    };


    bool allclose(double a[], double b[], double tolerance){
        for(int i = 0; i < row; i++){
            if( abs(a[i] - b[i]) > tolerance) return false;
            else continue;
        }
        return true;
    }

    void smooth(){
        double *x = new double[row];
        for(int i = 0; i < row; i++) x[i] = 0;
        for (int it_count = 0; it_count < ITERATION_LIMIT; it_count++){
            double *x_new = new double[row];
            for(int i = 0; i < row; i++) x_new[i] = 0;
            printf("Iteration %d: ",it_count);print_b(x);
            for (int i = 0; i < row; i++){
                float s1,s2;
                s1 = oneD_dot(i,x_new,0,i);
                s2 = oneD_dot(i,x,i+1,row);
                x_new[i] = (b[i] - s1 - s2) / A[i][i];
            }
            if (allclose(x,x_new,1e-8)){
                break;
            }
            x = x_new;
        }
        printf("Solution: ");print_b(x);
        double vague_b[row];
        twoD_dot(x,vague_b);
        double error[row];
        for(int i = 0; i < row; i++) error[i] = vague_b[i] - b[i];
        printf("\nError: ");print_b(error);
    }

    double* SPMV(double result[]){
        //double result[row];
        vector<int> row_start;
        vector<int> col_index;
        vector<double> value;

        for(int i = 0; i < row; ++i){
            row_start.push_back(col_index.size());
            for(int j = 0; j < col; ++j){
                if(A[i][j] != 0){
//                    printf("%f\n",A[i][j]);
                    col_index.push_back(j);
                    value.push_back((A[i][j]));
                }
            }
        }
        row_start.push_back(value.size());

        ////////print vector value for debug
//        printf("row_start:\n");
//        for(int i = 0; i < row_start.size(); ++i){
//            printf("%d ",row_start[i]);
//        }
//        printf("\n");
//
//        printf("col_index:\n");
//        for(int i = 0; i < col_index.size(); ++i){
//            printf("%d ",col_index[i]);
//        }
//        printf("\n");
//
//        printf("value:\n");
//        for(int i = 0; i < value.size(); ++i){
//            printf("%f ",value[i]);
//        }
//        printf("\n");
        ////////print vector value for debug

        for(int i = 0; i < row; ++i){
            for(int j = row_start[i]; j < row_start[i+1]; ++j){
                result[i] += value[j] * b[col_index[j]];
            }
//            printf("result = %f",result[i]);
        }

        printf("Result of SPMV:\n");
        print_b(result);

        return result;
    };

};


int main(){

    coarse_matrix coarse_matrix(strength_Matrix_row,strength_Matrix_col);


    int coarse_data[strength_Matrix_row][strength_Matrix_col]= {{3,5,5,5,3},{5,8,8,8,5},
                                       {5,8,8,8,5},{5,8,8,8,5},{3,5,5,5,3}};

    // int data[Matrix_row][Matrix_col]= {{3,3,3,3,3},{5,5,5,5,5},
    // {5,5,5,5,5},{5,5,5,5,5},{3,3,3,3,3}};

    coarse_matrix.input_data(coarse_data);
//    strength_matrix.input_data_inclass();


    // strength_matrix.v[4][4].value = 7;
    //printf(strength_matrix.v[0][0].fresh? "true\n" : "false\n");

    //printf("%d\n",strength_matrix.v[4][4].value);
    // printf("%d\n",data[1][2]);

    printf("Before Coarse:");
    coarse_matrix.print_value();
    coarse_matrix.coarse();
    coarse_matrix.print_stats();
    printf("\nAfter Coarse:");
    coarse_matrix.print_value();

    vector<vector<int>> cpoint = coarse_matrix.cpoint_coor;
    printf("Cpoints:\n");
    for(int i = 0; i < cpoint.size(); i++){
        printf("(%d,%d) ",cpoint[i][0],cpoint[i][1]);

    }

    printf("\n*******smmoth*********\n");
    double smooth_data[smooth_Matrix_row][smooth_Matrix_col]= {{10,-1,2,0},{-1,11,-1,3},
                                                            {2,-1,10,-1},{0,3,-1,8}};
    double RHS_data[smooth_Matrix_row] = {6,25,-11,15};
    class smooth_matrix smooth_matrix(smooth_Matrix_row,smooth_Matrix_col);
    smooth_matrix.input_data(smooth_data,RHS_data);
    smooth_matrix.print_value();

    smooth_matrix.smooth();


    printf("\n********SPMV********\n");
    double SPMV_A[4][4]= {{1,2,0,0},
                          {3,8,4,0},
                          {0,5,6,0},
                          {0,0,4,3}};
    double SPMV_v[4] = {1,1,1,1};

    class smooth_matrix spmv(4,4);
    spmv.input_data(SPMV_A,SPMV_v);
    double result[smooth_Matrix_row];
    spmv.SPMV(result);
}