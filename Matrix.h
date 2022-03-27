#pragma once

#ifdef MATRIX_EXPORTS
#define MATRIX_API __declspec(dllexport)
#else
#define MATRIX_API __declspec(dllimport)
#endif


#include <vector>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <initializer_list>




const double EPS = 1E-8;


class MATRIX_API IndexError {
private:
    std::string ind_error;

public:
    IndexError(std::string error)
        : ind_error(error)
    {
    }

    const char* getError() { return IndexError::ind_error.c_str(); }

};

class MATRIX_API MatrixError {
private:
    std::string mtr_error;

public:
    MatrixError(std::string error)
        : mtr_error(error)
    {
    }

    const char* getError() { return MatrixError::mtr_error.c_str(); }

};




class MATRIX_API FileError {
private:
    std::string fl_error;

public:
    FileError(std::string error)
        : fl_error(error)
    {
    }

    const char* getError() { return FileError::fl_error.c_str(); }

};








class MATRIX_API Matrix {
    friend class PCA;
public:
    int strk; // количество строк
    int stlb; // количество столбцов
    std::vector<std::vector<double>> values;

public:

    // конструкторы
    Matrix(int a, int b)
    {
        // a - количество строк, b - количество столбцов
        strk = a;
        stlb = b;
        for (int i = 0; i < a; i++) {
            std::vector<double> tmp;
            values.push_back(tmp);
            for (int j = 0; j < b; j++) {
                values[i].push_back(0);
            }
        }
    };

    Matrix()
    {
        strk = 0;
        stlb = 0;

    };


    Matrix(std::vector<std::vector<double>> vals) {


        for (int i = 0; i < vals.size(); i++) {
            if (vals[0].size() != vals[i].size()) {
                throw IndexError("Matrix.cell(i, j): wrong vector sizes in Matrix");
            }
        }
        strk = vals.size();
        stlb = vals[0].size();
        values = vals;

    };

    Matrix(std::initializer_list <double> a) {
        if (sqrt(a.size()) == (int)sqrt(a.size())) {
            strk = sqrt(a.size());
            stlb = sqrt(a.size());

            const double* a_ptr = a.begin();

            for (int i = 0; i < stlb; i++) {
                std::vector<double> tmp;
                values.push_back(tmp);
                for (int j = 0; j < strk; j++) {
                    values[i].push_back(*a_ptr);
                    a_ptr++;
                }
            }

        }
        else {
            throw MatrixError("Need 1,4,9,16,25,etc. cound of numbers.");
        }
    }

    // чтоб изменять значение в ячейке
    double& cell(int i, int j) { // строка, столбец
        if (i >= strk || j >= stlb || i < 0 || j < 0) {
            throw IndexError("Matrix.cell(i, j): no such indexes in Matrix");
        }

        return values[i][j];
    };

    // чтоб получить значение в ячейке с запретом на изменение
    double print(int i, int j) const { // строка, столбец
        if (i >= strk || j >= stlb || i < 0 || j < 0) {
            throw IndexError("Matrix.cell(i, j): no such indexes in Matrix");
            return -1;
        }

        return values[i][j];
    };



    // сумма матриц
    Matrix operator+(const Matrix& tmp) {
        if (this->stlb != tmp.stlb || this->strk != tmp.strk) {
            throw MatrixError("Matrix + Matrix: wrong sizes");
        }

        std::vector<std::vector<double>> vectmp;

        for (int i = 0; i < this->strk; i++) {
            std::vector<double> temp;
            vectmp.push_back(temp);
            for (int j = 0; j < this->stlb; j++) {
                vectmp[i].push_back(tmp.values[i][j] + this->values[i][j]);
            }
        }

        return Matrix(vectmp);

    }

    // вычитание матриц
    Matrix operator-(const Matrix& tmp) {
        if (this->stlb != tmp.stlb || this->strk != tmp.strk) {
            throw MatrixError("Matrix - Matrix: wrong sizes");
        }
        std::vector<std::vector<double>> vectmp;

        for (int i = 0; i < this->strk; i++) {
            std::vector<double> temp;
            vectmp.push_back(temp);
            for (int j = 0; j < this->stlb; j++) {
                vectmp[i].push_back(tmp.values[i][j] - this->values[i][j]);
            }
        }
        return Matrix(vectmp);

    }

    // умножение матриц
    Matrix operator*(const Matrix& tmp) const {
        if (this->stlb != tmp.strk) {
            throw MatrixError("Matrix * Matrix: wrong sizes");
        }
        std::vector<std::vector<double>> vectmp;


        double tmpsum = 0;

        for (int i = 0; i < this->strk; i++) {
            std::vector<double> temp;
            vectmp.push_back(temp);
            for (int j = 0; j < tmp.stlb; j++) {

                for (int k = 0; k < this->stlb; k++) {

                    tmpsum += this->print(i, k) * tmp.print(k, j);
                    //std::cout << tmpsum << "t";
                }

                /*if (fabs(tmpsum) < EPS) {
                    tmpsum = 0;

                }*/

                vectmp[i].push_back(tmpsum);
                tmpsum = 0;
            }



        }




        return Matrix(vectmp);
    }

    // умножение на число
    Matrix operator*(double num) {
        std::vector<std::vector<double>> tmpvec;

        for (int i = 0; i < this->strk; i++) {
            std::vector<double> temp;
            tmpvec.push_back(temp);
            for (int j = 0; j < this->stlb; j++) {
                tmpvec[i].push_back(this->values[i][j] * num);
            }


        }

        return Matrix(tmpvec);
    }

    // деление матрицы на число
    Matrix operator/(double num) {
        double tmp;

        if (num == 0) {
            throw IndexError(" can't divide to zero ");
        }
        std::vector<std::vector<double>> tmpvec;

        for (int i = 0; i < this->strk; i++) {
            std::vector<double> temp;
            tmpvec.push_back(temp);
            for (int j = 0; j < this->stlb; j++) {
                //if (fabs(this->values[i][j] / num) < EPS) {
                //    
                //    tmp = 0;
                //}
                //else {
                tmp = this->values[i][j] / num;
                //}
                tmpvec[i].push_back(tmp);

            }


        }
        return Matrix(tmpvec);
    }





    // константные методы чтоб получить число строк/столбцов, изменять нельзя
    int stroks() const {
        return this->strk;
    }

    int stlbs() const {
        return this->stlb;
    }



    // проверка будет ли матрица квадратной
    bool issquare() const {
        if (this->stroks() == this->stlbs()) {
            return true;
        }
        else {
            return false;
        }

    };




    // след матрицы
    double trace() const {
        if (this->issquare()) {
            double sum = 0;
            for (int i = 0; i < this->stroks(); i++) {
                sum += this->print(i, i);
            }
            return sum;
        }
        else {
            throw MatrixError("Matrix must be square");
        }
    }



    //определитель матрицы методом Гаусса

    double det() const {
        if (this->issquare() == false) {
            throw MatrixError("Matrix must be square");
        }
        else {
            Matrix temp(this->values);
            int n = temp.stroks();
            double d = 1;
            for (int i = 0; i < n; ++i) {
                int k = i;
                for (int j = i + 1; j < n; ++j)
                    if (abs(temp.print(j, i)) > abs(temp.print(k, i))) // поиск максимального элемента в столбце
                        k = j;
                if (abs(temp.print(k, i)) < EPS) { // если максимальный элемент в столбце очень мал, то определитель 0
                    d = 0;
                    break;
                }
                for (int j = 0; j < n; j++) { //  переставляем максимальную строку
                    // с максимальным по модулю элементом наверх
                    std::swap(temp.cell(i, j), temp.cell(k, j));
                }

                if (i != k) // меняем знак каждый раз когда переставляем строки
                    d = -d; // если максимальный элемент оказался в нужном месте, ничего не делаем
                d *= temp.cell(i, i); // умножаем det на элемент на главной диагонали
                for (int j = i + 1; j < n; ++j)
                    temp.cell(i, j) /= temp.cell(i, i); // делим элементы в строке на максимальный
                for (int j = 0; j < n; ++j)
                    if (j != i && abs(temp.cell(j, i)) > EPS)
                        for (int k = i + 1; k < n; ++k) {
                            temp.cell(j, k) -= temp.cell(i, k) * temp.cell(j, i);
                        }
            }
            return d;
        }
    }



    double fbnorm() const {
        double norm = 0;
        for (int i = 0; i < this->stlbs(); i++) {
            for (int j = 0; j < this->stroks(); j++) {
                norm += this->values[i][j] * this->values[i][j];
            }
        }
        return sqrt(norm);
    }










    Matrix T() {
        int strk1 = this->stlb;
        int stlb1 = this->strk;
        std::vector<std::vector<double>> tmp;
        for (int j = 0; j < this->stlb; j++) {
            std::vector<double> temp;
            tmp.push_back(temp);
            for (int i = 0; i < this->strk; i++) {
                tmp[j].push_back(this->cell(i, j));
            }
        }
        return Matrix(tmp);

    }

    int rank() {

        int i, j, k, l;
        double r;
        std::vector<std::vector<double>> a = this->values;


        i = 0; j = 0;
        while (i < this->strk && j < this->stlb) {
            // Инвариант: минор матрицы в столбцах 0..j-1
            //   уже приведен к ступенчатому виду, и строка
            //   с индексом i-1 содержит ненулевой эл-т
            //   в столбце с номером, меньшим чем j

            // Ищем максимальный элемент в j-м столбце,
            // начиная с i-й строки
            r = 0.0;
            for (k = i; k < this->strk; ++k) {
                if (fabs(a[k][j]) > r) {
                    l = k;      // Запомним номер строки
                    r = fabs(a[k][j]); // и макс. эл-т
                }
            }
            if (r <= EPS) {
                // Все элементы j-го столбца по абсолютной
                // величине не превосходят eps.
                // Обнулим столбец, начиная с i-й строки
                for (k = i; k < this->strk; ++k) {
                    a[k][j] = 0.0;
                }
                ++j;      // Увеличим индекс столбца
                continue; // Переходим к следующей итерации
            }

            if (l != i) {
                // Меняем местами i-ю и l-ю строки
                for (k = j; k < this->stlb; ++k) {
                    r = a[i][k];
                    a[i][k] = a[l][k];
                    a[l][k] = (-r); // Меняем знак строки
                }
            }

            // Утверждение: fabs(a[i*n + k]) > eps

            // Обнуляем j-й столбец, начиная со строки i+1,
            // применяя элем. преобразования второго рода
            for (k = i + 1; k < this->strk; ++k) {
                r = (-a[k][j] / a[i][j]);

                // К k-й строке прибавляем i-ю, умноженную на r
                a[k][j] = 0.0;
                for (l = j + 1; l < this->stlb; ++l) {
                    a[k][l] += r * a[i][l];
                }
            }

            ++i; ++j;   // Переходим к следующему минору
        }

        return i; // Возвращаем число ненулевых строк
    }



    // обратная матрица
    Matrix inv()
    {

        if (this->issquare() == false) {
            throw MatrixError("must be square matrix");
        }

        if (this->det() == 0) {
            throw MatrixError("Matrix is degenerate (Vurozhdena)");
        }



        std::vector<std::vector<double>> mtr = this->values; // копируем исходную матрицу
        std::vector<std::vector<double>> tmp; // составим тут единичную матрицу

        int n = this->values.size();
        double num;
        for (int i = 0; i < n; i++) {
            std::vector<double> temp;
            tmp.push_back(temp);
            for (int j = 0; j < n; j++)
            {
                tmp[i].push_back(0.0);

                if (i == j)
                    tmp[i][j] = 1.0;
            }
        }// единичную матрицу составили




        // методом гаусса приводим исходную к единичной матрице, при этом повторяем действия по строком с единичной
        for (int k = 0; k < n; k++)
        {
            num = mtr[k][k];

            for (int j = 0; j < n; j++)
            {
                mtr[k][j] /= num;
                tmp[k][j] /= num;
            }

            for (int i = k + 1; i < n; i++)
            {
                num = mtr[i][k];

                for (int j = 0; j < n; j++)
                {
                    mtr[i][j] -= mtr[k][j] * num;
                    tmp[i][j] -= tmp[k][j] * num;

                }
            }
        }

        for (int k = n - 1; k > 0; k--)
        {
            for (int i = k - 1; i >= 0; i--)
            {
                num = mtr[i][k];

                for (int j = 0; j < n; j++)
                {
                    mtr[i][j] -= mtr[k][j] * num;
                    tmp[i][j] -= tmp[k][j] * num;



                }
            }
        }




        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (fabs(tmp[i][j]) < EPS) {
                    tmp[i][j] = 0;
                }
            }
        }


        return Matrix(tmp);
    }


    int write(std::string filename)
    {
        std::ofstream out(filename, std::ios::binary);    //Ставим режим "бинарный файл"
        if (!out) {
            throw FileError("Error while openning file to write");
        }
        out.write((char*)&this->strk, sizeof(int));       //Записываем в файл колво столбцов
        out.write((char*)&this->stlb, sizeof(int));       //Записываем в файл колво строк
        for (int i = 0; i < this->strk; i++) {
            for (int j = 0; j < this->stlb; j++) {
                out.write((char*)&this->values[i][j], sizeof(double));
            }

        }
        out.close();
        return 0;
    }




    int read(std::string filename)
    {
        std::ifstream in(filename, std::ios::binary);            //Ставим режим "бинарный файл"
        if (!in) {
            throw FileError("Error while openning file to read");
        }

        for (int i = 0; i < this->strk; i++) {
            this->values[i].clear();
        }
        this->values.clear(); // будет ли утечка памяти если массив двойной?




        in.read((char*)&this->strk, sizeof(int));        //
        if (!in) throw FileError("wrong file");
        in.read((char*)&this->stlb, sizeof(int));        //
        if (!in) throw FileError("wrong file");



        for (int i = 0; i < this->strk; i++) {
            std::vector<double> temp;
            this->values.push_back(temp);
            for (int j = 0; j < this->stlb; j++) {
                double tmp;
                in.read((char*)&tmp, sizeof(double));
                if (!in) throw FileError("not enough elements or wrong file");
                // если элемент не был прочитан или прочитан неправильно или элементы кончились
                // выдаём ошибку

                this->values[i].push_back(tmp);
            }

        }
        in.close();
        return 0;


    }












    friend Matrix operator*(int num, const Matrix& tmp);
    friend Matrix adamar(const Matrix& tmp1, const Matrix& tmp2);
    friend std::ostream& operator<< (std::ostream& out, const Matrix& mtr);
    friend std::istream& operator>> (std::istream& in, Matrix& mtr); // доделать


};



Matrix operator*(int num, const Matrix& tmp) {
    std::vector<std::vector<double>> tmpvec;

    for (int i = 0; i < tmp.strk; i++) {
        std::vector<double> temp;
        tmpvec.push_back(temp);
        for (int j = 0; j < tmp.stlb; j++) {
            tmpvec[i].push_back(tmp.values[i][j] * num);
        }
    }


    return Matrix(tmpvec);
}


Matrix adamar(const Matrix& tmp1, const Matrix& tmp2) {
    if (tmp1.stlb != tmp2.stlb || tmp1.strk != tmp2.strk) {
        throw MatrixError("Matrix adamar: wrong sizes");
    }
    std::vector<std::vector<double>> vectmp;

    for (int i = 0; i < tmp1.strk; i++) {
        std::vector<double> temp;
        vectmp.push_back(temp);
        for (int j = 0; j < tmp1.stlb; j++) {
            vectmp[i].push_back(tmp1.values[i][j] * tmp2.values[i][j]);
        }
    }

    return Matrix(vectmp);
}


std::ostream& operator<< (std::ostream& out, const Matrix& mtr)
{
    int size = 0;



    size = 7;


    out.setf(std::ios::left);
    out << std::setprecision(5);



    for (int j = 0; j < mtr.strk; j++) {
        for (int i = 0; i < mtr.stlb; i++) {
            out << std::setfill(' ') << std::setw(size);
            out << mtr.print(j, i) << " ";
        }
        out << "\n";
    }


    return out;

}


std::istream& operator>> (std::istream& in, Matrix& mtr) {

    std::string line;
    for (int i = 0; i < mtr.strk; i++) {
        mtr.values[i].clear();
    }
    mtr.values.clear();

    std::vector<std::vector<double>> vec;

    int strk = 0;
    int stlb = 0;
    int tmpstlb = 0;
    std::string tmp;


    while (getline(in, line)) {
        std::replace(line.begin(), line.end(), ',', '.'); // меняет запятые на точки
        std::stringstream ss;
        ss << line;
        std::vector<double> temp;                          // чтоб формат double с . и , читались
        vec.push_back(temp);


        while (!ss.eof()) {
            ss >> tmp;
            if (strk == 0) { // если заполняется первая строка 
                stlb += 1;
                tmpstlb += 1;
            }
            else {
                tmpstlb += 1;
            }
            vec[strk].push_back(std::stof(tmp));
        }
        if (tmpstlb != stlb) {
            throw FileError{ "file have wrong sizes to be matrix" };
        }

        tmpstlb = 0;

        strk += 1;

    }



    /*
    std::string line;
    std::string tmp;
    while (getline(inf, line)) {
        std::replace(line.begin(), line.end(), ',', '.');
        std::stringstream ss;
        ss << line;
        while (!ss.eof()) {
            ss >> tmp;
            std::cout << std::stof(tmp);
        }

        //std::cout << line;
        tmp = "";
    }
    */


    mtr.stlb = stlb;
    mtr.strk = strk;
    mtr.values = vec;



    return in;


}






class MATRIX_API Vector : public Matrix {
public:
    friend class PCA;
    Vector() {
        stlb = 0;
        strk = 0;
    };


    Vector(std::vector<double> vec) {
        strk = vec.size();
        stlb = 1;

        for (int i = 0; i < vec.size(); i++) {
            std::vector<double> tmp;
            values.push_back(tmp);
            values[i].push_back(vec[i]);
        }

    }


    Vector(Matrix A, int j) { // матрица и номер столбца

        for (int i = 0; i < A.stroks(); i++) {
            std::vector<double> tmp;
            values.push_back(tmp);
            values[i].push_back(A.cell(i, j));
        }
        strk = A.stroks();
        stlb = 1;


    }





    double operator*(const Vector& tmp) {
        if (this->strk != tmp.strk) {
            throw MatrixError("vector * vector: wrong sizes");
        }

        double sum = 0;

        for (int i = 0; i < this->strk; i++) {
            sum += this->values[i][0] * tmp.values[i][0];
        }

        return sum;

    }


    Matrix operator*(const Matrix& tmp) {

        return Matrix(this->values) * tmp;
    }




    double evklidnorm() {
        return sqrt(*this * *this);
    }



    double maxnorm() {
        double max = 0;
        for (int i = 0; i < this->strk; i++) {
            for (int j = 0; j < this->stlb; j++) {
                if (abs(this->values[i][j]) > max) {
                    max = abs(this->values[i][j]);
                }
            }
        }
        return max;
    }



    friend double angle(Vector& v1, Vector& v2);


};

double angle(Vector& v1, Vector& v2) {

    return  acos((v1 * v2) / (v1.evklidnorm() * v2.evklidnorm()));


}




// написать подклассы для матриц указанных видов -  
//  верхняя и нижняя треугольные матрицы, симметричная матрица.


class MATRIX_API eye : public Matrix {
public:
    eye(int a) {


        strk = a;
        stlb = a;

        for (int i = 0; i < a; i++) {
            std::vector<double> tmp;
            values.push_back(tmp);
            for (int j = 0; j < a; j++) {
                if (i == j) {
                    values[i].push_back(1);
                }
                else {
                    values[i].push_back(0);
                }
            }
        }
    };

};


class MATRIX_API diagonal : public Matrix {
public:
    diagonal(std::vector<double> vec) {

        strk = vec.size();
        stlb = vec.size();

        for (int i = 0; i < strk; i++) {
            std::vector<double> tmp;
            values.push_back(tmp);
            for (int j = 0; j < strk; j++) {
                if (i == j) {
                    values[i].push_back(vec[i]);
                }
                else {
                    values[i].push_back(0);
                }
            }
        }
    };

    diagonal(Matrix mtr) {
        if (mtr.issquare() != true) {
            throw MatrixError("diagonal: matrix must be square");
        }
        for (int i = 0; i < mtr.stroks(); i++) {
            std::vector<double> tmp;
            values.push_back(tmp);
            for (int j = 0; j < mtr.stroks(); j++) {
                if (i == j) {
                    values[i].push_back(mtr.cell(i, j));
                }
                else {
                    values[i].push_back(0);
                }
            }
        }
        strk = mtr.stroks();
        stlb = strk;
    }

};


class MATRIX_API upper : public Matrix {
public:
    upper(int c) {
        strk = c;
        stlb = c;
        for (int i = 0; i < c; i++) {
            std::vector<double> tmp;
            values.push_back(tmp);
            for (int j = 0; j < c; j++) {
                if (j >= i) {
                    double a = (double)(rand() % 10000) / 1000;
                    values[i].push_back(a);
                }
                else {
                    values[i].push_back(0);
                }
            }
        }
    };

    double& cell(int i, int j) {
        if (i > j) {
            throw MatrixError("This value is constant");
        }
        else {
            return Matrix::cell(i, j);
        }

    }


};


class MATRIX_API lower : public Matrix {
public:
    lower(int c) {
        strk = c;
        stlb = c;
        for (int i = 0; i < c; i++) {
            std::vector<double> tmp;
            values.push_back(tmp);
            for (int j = 0; j < c; j++) {
                if (i >= j) {
                    double a = (double)(rand() % 10000) / 1000;
                    values[i].push_back(a);
                }
                else {
                    values[i].push_back(0);
                }
            }
        }
    };

    double& cell(int i, int j) {
        if (j > i) {
            throw MatrixError("This value is constant");
        }
        else {
            return Matrix::cell(i, j);
        }

    }

};


class MATRIX_API pair { // чисто для симметричной матрицы существует
    // чтобы в симметричной матрице менялись оба элемента всегда.
    // если это не нужно можно преобразовать симметричную к обычной матрице уже после генерации
public:
    double* a;
    double* b;
    pair(double* tmp1, double* tmp2) {
        a = tmp1;
        b = tmp2;
    }


    void operator= (const double num) {
        *a = num;
        *b = num;
    }

};


class MATRIX_API simetric : public Matrix {
public:
    simetric(int c) {
        strk = c;
        stlb = c;
        for (int i = 0; i < c; i++) {
            std::vector<double> tmp;
            values.push_back(tmp);
            for (int j = 0; j < c; j++) {
                if (j >= i) {
                    double a = (double)(rand() % 10000) / 1000;
                    values[i].push_back(a);
                }
                else {
                    values[i].push_back(this->print(j, i)); // всегда сработает так как матрица заполняется построчно.
                    // то есть рандомится вся первая строка, вторая строка вся рандомится кроме одного числа которое берётся из первой строки
                    // в третьей строки уже 2 числа берутся из предыдущих строк а остальные рандомятся и так далее, ошибки никогда не возникнет.
                }
            }
        }
    };




    pair cell(int i, int j) {
        return pair(&Matrix::cell(i, j), &Matrix::cell(j, i));
    }



};













/*Дополнить библиотеку классов, разработанную в рамках ЛР1 и ЛР2. Реализовать операции над векторами и матрицами:

Угол между векторами  (done)   angle(v1, v2)

Ранг Матрицы (с помощью алгоритма Гаусса)   (done) rank()

Обратная матрица (если существует)   

Транспонирование матрицы (done)  Matrix.T()


*/













class MATRIX_API PCA{
protected:
public:
    Matrix X;
    Matrix T;
    Matrix P;
    Matrix E;
    int components = 0;
    friend class Matrix;
    friend class Vector;
    PCA(Matrix A, int a) {
        X = Matrix(A.values);
        E = X;
        components = a;
    }


    double m(int j) {
        double sum = 0;

        for (int i = 0; i < X.stroks(); i++) {
            sum += X.cell(i, j);
        }

        return sum / X.stroks();
    }



    void center() {
        for (int i = 0; i < E.stlbs(); i++) {
            
            double sum = m(i);
            
            for (int j = 0; j < E.stroks(); j++) {
                E.cell(j ,i) = E.print(j, i) - sum;
            }
         
        
        }
        

    }
   

    double s(int j) {
        double sum = 0;
        double mtmp = m(j);
 
        for (int i = 0; i < E.stroks(); i++) {
            sum += pow((E.cell(i, j) - mtmp), 2);
        }

        sum = sum / (E.stroks() - 1);

        return sqrt(sum);

    }

    
    void scaling() {
        
        for (int j = 0; j < X.stlbs(); j++) {
            int sum = 0;
            double stmp = s(j);

            for (int i = 0; i < E.stroks(); i++) {
                E.cell(i, j) = E.print(i, j) / stmp;
            }



        }





    }




    void autoshkal() {
        for (int j = 0; j < X.stlbs(); j++) {
            double mtmp = m(j);
            double stmp = s(j);
            for (int i = 0; i < X.stroks(); i++) {
                E.cell(i, j) = (X.print(i, j) - mtmp) / stmp;
                
            }
            //std::cout << stmp << " / ";
        }


    }


    void main_func() {
        double eps = 1E-8;
        Vector d;
        Vector p;
        

        for (int h = 0; h < components; h++) {
            Vector t = Vector(E, h); // из матрицы E берём столбец с индесом h

            do {
                p = Vector(((t.T() * E) / (t * t)).T(), 0);
                //std::cout << p << std::endl;
                p = Vector((p / p.evklidnorm()), 0) ;
                
                Vector t_old = t;

                t = Vector((E * p) / (p * p), 0);

                d = Vector(t_old - t, 0);
                //std::cout << d.evklidnorm() << std::endl;
                //std::cout << t_old - t;
            } while (d.evklidnorm() > eps);
            //std::cout << E;
            E = E - (Matrix)t * p.T();
            

            T.values.push_back(t.T().values[0]);
            T.stlb = t.T().values[0].size();
            T.strk++;
            P.values.push_back(p.T().values[0]);
            P.stlb = p.T().values[0].size();
            P.strk++;


        }

        //std::cout << T;
        T = T.T();
        P = P.T();

    }



    void print2() {
        std::cout << X;
    }



    void print() {
        std::cout << E; 
    }
    
    void check() {
        std::cout << "\n\n\n";
        std::cout << T;
        std::cout << "\n\n\n";
        std::cout << P;
        std::cout << "\n\n\n";
        std::cout << T * P.T() << std::endl;
        
    }


    Matrix razmahi() {
        std::vector<std::vector<double>> tmp;
        std::vector<double> razmahs;
        
        for (int i = 0; i < T.values.size(); i++) {
            razmahs.push_back((Vector(T.T(), i).T() * (T.T() * T).inv() * Vector(T.T(), i)).values[0][0]);
            //std::cout << Vector(T.T(), i).T() * (T.T() * T).inv() * Vector(T.T(), i);
        }
        tmp.push_back(razmahs);
        //std::cout << Matrix(tmp);
        return Matrix(tmp);
    }



    Matrix dispersion() { // сумма квадратов по строкам матрицы E
        std::vector<std::vector<double>> tmp;
        std::vector<double> dspr;
        for (int i = 0; i < E.values.size(); i++) {
            double sum = 0;
            for (int j = 0; j < E.values[0].size(); j++) {
                sum += pow(E.values[i][j], 2) ;
            }
            dspr.push_back(sum); // сохарняем сумму остатокв
            //std::cout << Vector(T.T(), i).T() * (T.T() * T).inv() * Vector(T.T(), i);
        }
        tmp.push_back(dspr);
        //std::cout << Matrix(tmp);
        return Matrix(tmp);
    }


    
    double TRV() {
        Matrix tmp = PCA::dispersion();
        double v0 = 0;
        for (int i = 0; i < tmp.values[0].size(); i++) {
            v0 += tmp.values[0][i];
        }

        v0 /= tmp.values[0].size();

        return (v0 / components);

    }


    double ERV() {
        Matrix tmp = PCA::dispersion();
        double v0 = 0;
        for (int i = 0; i < tmp.values[0].size(); i++) {
            v0 += tmp.values[0][i];
        }
        
        double sum = 0;
        for (int i = 0; i < X.values.size(); i++) {
            for (int j = 0; j < X.values[0].size(); j++) {
                sum += pow(X.values[i][j], 2);
                

            }
        }
        return (1 - v0 / sum);


    }











};










