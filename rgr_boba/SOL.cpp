#include "SOL.h"

using namespace std;

//Ввод матрицы
void SOL::Matrin() {
    ifstream file("file.txt");
    float si1, sj1;
    file >> si1 >> sj1;
    sizei = (int)si1;
    sizej = (int)sj1;
    for (int i = 0; i < sizei; i++) {
        vector<float> Avs;
        for (int j = 0; j < sizej; j++) {
            float  buff;
            file >> buff;
            Avs.push_back(buff);
        }
        X.push_back(Avs);
    }
};

void SOL::print(vector<vector<float>> M) {
    for (int i = 0; i < sizei; i++) {
        for (int j = 0; j < sizej; j++) {
            cout << setiosflags(ios::right) << setprecision(5);
            cout << setw(10) << M[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;
};

//Проверка на диагональное преобладание
bool SOL::checkDiag(vector<vector<float>> C) {
    vector<vector<float>> B = C;
    for (int i = 0; i < sizei; i++) {
        B[i][i] = 0;
        float buff = 0;
        for (int j = 0; j < sizei; j++) buff += abs(B[i][j]);
        if (abs(C[i][i]) <= buff)return false;
    }
    return true;
};

//Запись формул для простой итерации
void SOL::matrixFormula() {
    for (int i = 0; i < sizei; i++) {
        float buff = -X[i][i];
        for (int j = 0; j < sizej; j++) if (X[i][j] != 0)X[i][j] /= buff;
        X[i][sizej - 1] *= -1;
        X[i][i] = 0;
    }
};

//Метод простой итерации
void SOL::simpleIt() {
    vector<vector<float>> aprox;
    vector<float> vsp;
    int count = 1;
    for (int i = 0; i < sizei; i++) vsp.push_back(X[i][sizej - 1]);
    for (int i = 0; i < sizei; i++) vsp.push_back(0);
    for (int i = 0; i < sizei; i++) {
        float buff1 = 0;
        for (int j = 0; j < sizei; j++) buff1 += vsp[j] * copy[i][j];
        vsp.push_back(buff1);
    }
    aprox.push_back(vsp);

    //Итерации
    while (1)
    {
        vector<float> wsp, wsp2;
        float buff = 0;
        for (int i = 0; i < sizei; i++) {
            buff = 0;
            for (int j = 0; j < sizei; j++) buff += X[i][j] * aprox[count - 1][j];
            buff += X[i][sizej - 1];
            wsp.push_back(buff);
        }
        for (int i = 0; i < sizei; i++) {
            buff = abs(wsp[i] - aprox[count - 1][i]) / abs(wsp[i]);
            wsp.push_back(buff);
        }
        count++;
        for (int i = 0; i < sizei; i++) {
            float buff1 = 0;
            for (int j = 0; j < sizei; j++) buff1 += wsp[j] * copy[i][j];
            wsp.push_back(buff1);
        }
        aprox.push_back(wsp);
        bool flag = true;
        for (int i = 0; i < sizei; i++) {
            if (aprox[count - 1][i + sizei] > epsilon)flag = false;
        }
        if (flag == true)break;
    }
    printSolution(aprox);
};

//Поиск матрицы с диагональным преобладанием
void SOL::DiagDomin()
{
    vector<vector<float>> B = X, C = X;

    //Перебор матриц по строчкам
    for (int g = 0; g < sizej; g++) {
        for (int i = 1; i < sizei; i++) {
            float buff1, buff2;
            for (int j = 0; j < sizej; j++) {
                buff1 = B[i][j];
                buff2 = B[i - 1][j];
                B[i - 1][j] = buff1;
                B[i][j] = buff2;
            }
            if (checkDiag(B)) {
                X = B;
                return;
            }
        }
    }

    for (int g = 0; g < sizej; g++) {
        for (int j = 1; j < sizei; j++) {
            float buff1, buff2;
            for (int i = 0; i < sizei; i++) {
                buff1 = C[i][j];
                buff2 = C[i][j - 1];
                B[i][j - 1] = buff1;
                B[i][j] = buff2;
            }
            if (checkDiag(C)) {
                X = C;
                return;
            }
        }
    }
    //Вывод ошибки
    throw runtime_error("ERROR");
};

//Метод Гаусса
void SOL::MethodGaus() {
    vector<vector<float>> G = X;
    float tmp;
    float* xx = new float[sizej];
    short int i, j, k;
    for (i = 0; i < sizei; i++)
    {
        tmp = G[i][i];
        for (j = sizei; j >= i; j--)
            G[i][j] /= tmp;
        for (j = i + 1; j < sizei; j++)
        {
            tmp = G[j][i];
            for (k = sizei; k >= i; k--)
                G[j][k] -= tmp * G[i][k];
        }
    }
    xx[sizei - 1] = G[sizei - 1][sizei];
    for (i = sizei - 2; i >= 0; i--)
    {
        xx[i] = G[i][sizei];
        for (j = i + 1; j < sizei; j++) xx[i] -= G[i][j] * xx[j];
    }
    for (int g = 0; g < sizei; g++) {
        cout << setw(10) << setprecision(5) << "x" << g + 1 << ": " << xx[g] << endl;
    }
};

//Вывод ответа
void SOL::printSolution(vector<vector<float>> Matrix) {
    cout << setiosflags(ios::right) << setw(13) << " Iteration | ";
    for (int i = 0; i < sizei; i++) {
        cout << setiosflags(ios::right) << setw(9) << "x" << i + 1 << " | ";
    }
    for (int i = 0; i < sizei; i++) {
        cout << setiosflags(ios::right) << setw(9) << "del" << i + 1 << " | ";
    }
    for (int i = 0; i < sizei; i++) {
        cout << setiosflags(ios::right) << setw(9) << "r" << i + 1 << " | ";
    }
    cout << endl;
    for (int i = 0; i < 4 * sizei; i++) cout << setw(10) << "-----------";
    if (Matrix.size() > 10) {
        for (int i = Matrix.size() - 10; i < Matrix.size(); i++) {
            cout << endl;
            cout << setiosflags(ios::right) << setprecision(5) << setw(10) << i << " | ";
            for (int j = 0; j < 2 * sizei; j++) {
                cout << setiosflags(ios::right) << setprecision(5) << setw(10) << Matrix[i][j] << " | ";
            }
            cout << endl;
            for (int g = 0; g < 4 * sizei; g++) cout << setw(10) << "----------";
            cout << endl;
        }
    }
    else {
        for (int i = 0; i < Matrix.size(); i++) {
            cout << endl;
            cout << setiosflags(ios::right) << setprecision(5) << setw(10) << i << " | ";
            for (int j = 0; j < 3 * sizei; j++) {
                cout << setiosflags(ios::right) << setprecision(5) << setw(10) << Matrix[i][j] << " | ";
            }
            cout << endl;
            for (int g = 0; g < 4 * sizei; g++) cout << setw(10) << "-----------";
        }
    }
};
