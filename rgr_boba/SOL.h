#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <iostream>

using namespace std;

class SOL {
private:
    vector<vector<float>> X;
    vector<vector<float>> copy;

    int sizei, sizej;
    float epsilon = 0.0001;

    void Matrin();
    void print(vector<vector<float>> M);
    bool checkDiag(vector<vector<float>> Matrix);
    void matrixFormula();
    void simpleIt();
    void DiagDomin();
    void MethodGaus();
    void printSolution(vector<vector<float>> Matrix);
public:
    SOL() {
        Matrin();
        copy = X;
        cout << endl << "   Original matrix in Ax+B:" << endl << endl;
        print(X);
        cout << endl << "   Solution:" << endl << endl;
        MethodGaus();
        if (!checkDiag(X))DiagDomin();
        matrixFormula();
        cout << endl << "   Matrix for solving by iterations:" << endl << endl;
        print(X);
        cout << endl << "   Solution:" << endl << endl;
        simpleIt();
        system("pause");
    }
};
