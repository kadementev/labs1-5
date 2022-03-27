#include <vector>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <math.h>

const double EPS = 1E-9;


class IndexError {
private:
    std::string ind_error;

public:
    IndexError(std::string error)
        : ind_error(error)
    {
    }

    const char* getError() { return IndexError::ind_error.c_str(); }

};




class Matrix {
    protected:
        int strk; // количество строк
        int stlb; // количество столбцов
        std::vector<std::vector<float>> values;

    public:

        // конструкторы
        Matrix(int a, int b)
        {
            // a - количество строк, b - количество столбцов
            strk = a;
            stlb = b;
            for (int i = 0; i < a; i++) {
                std::vector<float> tmp;
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


        Matrix(std::vector<std::vector<float>> vals) {


            for (int i = 0; i < vals.size(); i++) {
                if (vals[0].size() != vals[i].size()) {
                    throw IndexError("Matrix.cell(i, j): wrong vector sizes in Matrix");
                }
            }
            strk = vals.size();
            stlb = vals[0].size();
            values = vals;

        };



        // чтоб изменять значение в ячейке
        float& cell(int i, int j) { // строка, столбец
            if (i >= strk || j >= stlb || i < 0 || j < 0) {
                throw IndexError("Matrix.cell(i, j): no such indexes in Matrix");
            }

            return values[i][j] ;
        };
    
        // чтоб получить значение в ячейке с запретом на изменение
        float print(int i, int j) const { // строка, столбец
            if (i >= strk || j >= stlb || i < 0 || j < 0) {
                throw IndexError("Matrix.cell(i, j): no such indexes in Matrix");
                return -1;
            }

            return values[i][j];
        };



        // сумма матриц
        Matrix operator+(const Matrix& tmp) {
            if (this->stlb != tmp.stlb || this->strk != tmp.strk) {
                throw IndexError("Matrix + Matrix: wrong sizes");
            }

            std::vector<std::vector<float>> vectmp;

            for (int i = 0; i < this->strk; i++) {
                std::vector<float> temp;
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
                throw IndexError("Matrix - Matrix: wrong sizes");
            }
            std::vector<std::vector<float>> vectmp;

            for (int i = 0; i < this->strk; i++) {
                std::vector<float> temp;
                vectmp.push_back(temp);
                for (int j = 0; j < this->stlb; j++) {
                    vectmp[i].push_back(tmp.values[i][j] - this->values[i][j]);
                }
            }
            return Matrix(vectmp);

        }

        // умножение матриц
        Matrix operator*(const Matrix& tmp) const{
            if (this->stlb != tmp.strk) {
                throw IndexError("Matrix * Matrix: wrong sizes");
            }
            std::vector<std::vector<float>> vectmp;

            
            float tmpsum = 0;

            for (int i = 0; i < this->strk; i++) {
                std::vector<float> temp;
                vectmp.push_back(temp);
                for (int j = 0; j < tmp.stlb; j++) {

                    for (int k = 0; k < this->stlb; k++) {
                        tmpsum += this->print(i, k) * tmp.print(k, j);
                    }



                    vectmp[i].push_back(tmpsum);
                    tmpsum = 0;
                }



            }




            return Matrix(vectmp);
        }

        // умножение на число
        Matrix operator*(int num) {
            std::vector<std::vector<float>> tmpvec;

            for (int i = 0; i < this->strk; i++) {
                std::vector<float> temp;
                tmpvec.push_back(temp);
                for (int j = 0; j < this->stlb; j++) {
                    tmpvec[i].push_back(this->values[i][j] * num);
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
        float trace() const {
            if (this->issquare()) {
                float sum = 0;
                for (int i = 0; i < this->stroks(); i++) {
                    sum += this->print(i, i);
                }
                return sum;
            }
            else {
                throw IndexError("Matrix must be squared");
            }
        }



        //определитель матрицы методом Гаусса

        float det() const {
            if (this->issquare() == false) {
                throw IndexError("Matrix must be squared");
            }
            else {
                Matrix temp(this->values);
                int n = temp.stroks();
                float d = 1;
                for (int i = 0; i < n; ++i) {
                    int k = i;
                    for (int j = i + 1; j < n; ++j)
                        if (abs(temp.print(j,i)) > abs(temp.print(k,i))) // поиск максимального элемента в столбце
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
                        temp.cell(i,j) /= temp.cell(i, i); // делим элементы в строке на максимальный
                    for (int j = 0; j < n; ++j)
                        if (j != i && abs(temp.cell(j, i)) > EPS)
                            for (int k = i + 1; k < n; ++k) {
                                temp.cell(j, k) -= temp.cell(i, k) * temp.cell(j, i);
                            }
                }
                return d;
            }
        }



        float fbnorm() const {
            float norm = 0;
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
            std::vector<std::vector<float>> tmp;
            for (int j = 0; j < this->stlb; j++) {
                std::vector<float> temp;
                tmp.push_back(temp);
                for (int i = 0; i < this->strk; i++) {
                    tmp[j].push_back(this->cell(i, j));
                }
            }
            return Matrix(tmp);

        }





        friend Matrix operator*(int num, const Matrix& tmp);
        friend Matrix adamar(const Matrix& tmp1, const Matrix& tmp2);
        friend std::ostream& operator<< (std::ostream& out, const Matrix& mtr);

};



Matrix operator*(int num, const Matrix& tmp) {
    std::vector<std::vector<float>> tmpvec;

    for (int i = 0; i < tmp.strk; i++) {
        std::vector<float> temp;
        tmpvec.push_back(temp);
        for (int j = 0; j < tmp.stlb; j++) {
            tmpvec[i].push_back(tmp.values[i][j] * num);
        }
    }


    return Matrix(tmpvec);
}


Matrix adamar(const Matrix& tmp1, const Matrix& tmp2) {
    if (tmp1.stlb != tmp2.stlb || tmp1.strk != tmp2.strk) {
        throw IndexError("Matrix adamar: wrong sizes");
    }
    std::vector<std::vector<float>> vectmp;

    for (int i = 0; i < tmp1.strk; i++) {
        std::vector<float> temp;
        vectmp.push_back(temp);
        for (int j = 0; j < tmp1.stlb; j++) {
            vectmp[i].push_back(tmp1.values[i][j] * tmp2.values[i][j]);
        }
    }

    return Matrix(vectmp);
}




std::ostream& operator<< (std::ostream& out, const Matrix& mtr)
    {
    float max = mtr.values[0][0];
    float min = mtr.values[0][0];
    int size = 0;
    
    

    size = 6;


    out.setf(std::ios::left);
    out << std::setprecision(6);
    if (mtr.strk == 1) {
        out << "( ";
        for (int i = 0; i < mtr.stlb; i++) {
            out << std::setfill(' ') << std::setw(size);
            out << mtr.print(0, i) << " ";
        }
        out << ")\n";


    }

    else {

        out << "/ ";
        for (int i = 0; i < mtr.stlb; i++) {
            out << std::setfill(' ') << std::setw(size);
            out << mtr.print(0, i) << " ";
        }
        out << "\\ \n";

        for (int j = 1; j < mtr.strk - 1; j++) {
            out << "| ";
            for (int i = 0; i < mtr.stlb; i++) {
                out << std::setfill(' ') << std::setw(size);
                out << mtr.print(j, i) << " ";
            }
            out << "|\n";
        }

        out << "\\ ";
        for (int i = 0; i < mtr.stlb; i++) {
            out << std::setfill(' ') << std::setw(size);
            out << mtr.print(mtr.strk - 1, i) << " ";
        }
        out << "/ \n";

    }

    return out;

    }




class Vector : public Matrix {
public:

    Vector() {
        stlb = 0;
        strk = 0;
    };

    Vector(std::vector<float> vec) {
        strk = vec.size();
        stlb = 1;

        for (int i = 0; i < vec.size(); i++) {
            std::vector<float> tmp;
            values.push_back(tmp);
            values[i].push_back(vec[i]);
        }

    }

    

    
    float operator*(const Vector& tmp) {
        if (this->strk != tmp.strk) {
            throw IndexError("vector * vector: wrong sizes");
        }
        
        float sum = 0;

        for (int i = 0; i < this->strk; i++) {
            sum += this->values[i][0] * tmp.values[i][0];
        }

        return sum;

    }


    Matrix operator*(const Matrix& tmp) {
        
        return Matrix(this->values) * tmp;
    }




    float evklidnorm() {
        return sqrt(*this * *this);
    }



    float maxnorm() {
        float max = 0;
        for (int i = 0; i < this->strk; i++) {
            for (int j = 0; j < this->stlb; j++) {
                if (abs(this->values[i][j]) > max) {
                    max = abs(this->values[i][j]);
                }
            }
        }
        return max;
    }



    friend float angle(Vector& v1, Vector& v2);


};

float angle(Vector& v1, Vector& v2) {

    return  acos((v1 * v2) / (v1.evklidnorm() * v2.evklidnorm()));


}




// написать подклассы для матриц указанных видов -  
//  верхняя и нижняя треугольные матрицы, симметричная матрица.


class eye : public Matrix {
public:
    eye(int a) {


        strk = a;
        stlb = a;

        for (int i = 0; i < a; i++) {
            std::vector<float> tmp;
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


class diagonal : public Matrix {
public:
    diagonal(std::vector<float> vec) {

        strk = vec.size();
        stlb = vec.size();

        for (int i = 0; i < strk; i++) {
            std::vector<float> tmp;
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
            throw IndexError("diagonal: matrix must be square");
        }
        for (int i = 0; i < mtr.stroks(); i++) {
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


class upper : public Matrix {
public:
    upper(int c) {
        strk = c;
        stlb = c;
        for (int i = 0; i < c; i++) {
            for (int j = 0; j < c; j++) {
                if (j >= i) {
                    float a = (float)(rand() % 10000) / 1000;
                    values[i].push_back(a);
                }
                else {
                    values[i].push_back(0);
                }
            }
        }
    };

    float& cell(int i, int j) {
        if (i > j) {
            throw IndexError("This value is constant");
        }
        else {
            return Matrix::cell(i, j);
        }
        
    }


};

class lower : public Matrix {
public:
    lower(int c) {
        strk = c;
        stlb = c;
        for (int i = 0; i < c; i++) {
            for (int j = 0; j < c; j++) {
                if (i >= j) {
                    float a = (float)(rand() % 10000) / 1000;
                    values[i].push_back(a);
                }
                else {
                    values[i].push_back(0);
                }
            }
        }
    };

    float& cell(int i, int j) {
        if (j > i) {
            throw IndexError("This value is constant");
        }
        else {
            return Matrix::cell(i, j);
        }

    }

};





class pair { // чисто для симметричной матрицы существует
    // чтобы в симметричной матрице менялись оба элемента всегда.
    // если это не нужно можно преобразовать симметричную к обычной матрице уже после генерации
public:
    float* a;
    float* b;
    pair(float* tmp1, float* tmp2) {
        a = tmp1;
        b = tmp2;
    }
    

    void operator= (const float num) {
        *a = num;
        *b = num;
    }

};



class simetric : public Matrix {
public:
    simetric(int c) {
        strk = c;
        stlb = c;
        for (int i = 0; i < c; i++) {
            std::vector<float> tmp;
            values.push_back(tmp);
            for (int j = 0; j < c; j++) {
                if (j >= i) {
                    float a = (float)(rand() % 10000) / 1000;
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







/*Реализовать операции над векторами и матрицами :
* 
След матрицы (done) trace

Определитель матрицы(методом Гаусса) (done) det

Скалярное произведение векторов (done) vector * vector

Норма вектора(евклидова норма, максимальная норма) (done) evkildnorm, maxnorm

Норма матрицы(норма Фробениуса)  (done)
*/












/*Дополнить библиотеку классов, разработанную в рамках ЛР1 и ЛР2. Реализовать операции над векторами и матрицами:

Угол между векторами  (done)   angle(v1, v2)

Ранг Матрицы (с помощью алгоритма Гаусса)   

Обратная матрица (если существует)   

Транспонирование матрицы (done)  Matrix.T()


*/






