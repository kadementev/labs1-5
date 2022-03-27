#include <vector>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <iostream>


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


    public:
        std::vector<float> values;
        // конструкторы
        Matrix(int a, int b)
        {
            // a - количество строк, b - количество столбцов
            strk = a;
            stlb = b;
            for (int i = 0; i < a * b; i++) {
                values.push_back(0);
            }
        };

        Matrix()
        {
            strk = 0;
            stlb = 0;

        };


        Matrix(int a, int b, std::vector<float> vals) {

            if (a * b != vals.size()) {
                throw IndexError("Matrix.cell(i, j): no such indexes in Matrix");
            }

            strk = a;
            stlb = b;
            values = vals;

        };



        // чтоб изменить
        float& cell(int i, int j) { // строка, столбец
            if (i >= strk || j >= stlb || i < 0 || j < 0) {
                throw IndexError("Matrix.cell(i, j): no such indexes in Matrix");
            }

            return values[i * stlb + j];
        };
    
        // чтоб написать
        float print(int i, int j) const { // строка, столбец
            if (i >= strk || j >= stlb || i < 0 || j < 0) {
                throw IndexError("Matrix.cell(i, j): no such indexes in Matrix");
                return -1;
            }

            return values[i * stlb + j];
        };




        Matrix operator+(const Matrix& tmp) {
            if (this->stlb != tmp.stlb || this->strk != tmp.strk) {
                throw IndexError("Matrix + Matrix: wrong sizes");
            }

            std::vector<float> vectmp;

            for (int i = 0; i < tmp.values.size(); i++) {
                vectmp.push_back(tmp.values[i] + this->values[i]);
            }

            return Matrix(this->strk, this->stlb, vectmp);

        }


        Matrix operator-(const Matrix& tmp) {
            if (this->stlb != tmp.stlb || this->strk != tmp.strk) {
                throw IndexError("Matrix - Matrix: wrong sizes");
            }
            std::vector<float> vectmp;

            for (int i = 0; i < tmp.values.size(); i++) {
                vectmp.push_back(this->values[i] - tmp.values[i]);
            }
            return Matrix(this->strk, this->stlb, vectmp);

        }


        Matrix operator*(const Matrix& tmp) const{
            if (this->stlb != tmp.strk) {
                throw IndexError("Matrix * Matrix: wrong sizes");
            }
            std::vector<float> vectmp;

            
            float tmpsum = 0;

            for (int i = 0; i < this->strk; i++) {
                for (int j = 0; j < tmp.stlb; j++) {

                    for (int k = 0; k < this->stlb; k++) {
                        tmpsum += this->print(i, k) * tmp.print(k, j);
                    }



                    vectmp.push_back(tmpsum);
                    tmpsum = 0;
                }



            }




            return Matrix(this->strk, tmp.stlb, vectmp);
        }


        Matrix operator*(int num) {
            std::vector<float> tmpvec;

            for (int i = 0; i < this->values.size(); i++) {
                tmpvec.push_back(this->values[i] * num);
            }

            return Matrix(this->strk, this->stlb, tmpvec);
        }

        bool issquare() {
            if (this->strk == this->stlb) {
                return true;
            }
            else {
                return false;
            }
        
        };

        int stroks() const {
            return this->strk;
        }

        
        friend Matrix operator*(int num, const Matrix& tmp);
        friend Matrix adamar(const Matrix& tmp1, const Matrix& tmp2);
        friend std::ostream& operator<< (std::ostream& out, const Matrix& mtr);
};



Matrix operator*(int num, const Matrix& tmp) {
    std::vector<float> tmpvec;

    for (int i = 0; i < tmp.values.size(); i++) {
        tmpvec.push_back(tmp.values[i] * num);
    }


    return Matrix(tmp.strk, tmp.stlb, tmpvec);
}


Matrix adamar(const Matrix& tmp1, const Matrix& tmp2) {
    if (tmp1.stlb != tmp2.stlb || tmp1.strk != tmp2.strk) {
        throw IndexError("Matrix adamar: wrong sizes");
    }
    std::vector<float> vectmp;

    for (int i = 0; i < tmp1.values.size(); i++) {
        vectmp.push_back(tmp1.values[i] * tmp2.values[i]);
    }

    return Matrix(tmp1.strk, tmp2.stlb, vectmp);
}




std::ostream& operator<< (std::ostream& out, const Matrix& mtr)
    {
    float max = mtr.values[0];
    float min = mtr.values[0];
    int size = 0;
    
    /*
    for (int i = 0; i < mtr.values.size(); i++) {
        if (mtr.values[i] > max) {
            max = mtr.values[i];
        }
        if (mtr.values[i] < min) {
            min = mtr.values[i];
        }
    };


    
    if (std::to_string(max).size() > std::to_string(min).size()) {
        size = std::to_string(max).size();

        }
    else {
        size = std::to_string(min).size();
    }
    */

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
        out << "\\\ \n";

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

    Vector(int a, std::vector<float> vec) {
        
        values = vec;
        strk = a;
        stlb = 1;
    }

    

    
    int operator*(const Vector& tmp) {
        if (this->strk != tmp.strk) {
            throw IndexError("Matrix adamar: wrong sizes");
        }
        
        int sum = 0;

        for (int i = 0; i < this->strk; i++) {
            sum += this->values[i] * tmp.values[i];
        }

        return sum;

    }


    Matrix operator*(const Matrix& tmp) {
        
        return Matrix(this->stlb, this->strk, this->values) * tmp;
    }


};






// написать подклассы для матриц указанных видов -  
//  верхняя и нижняя треугольные матрицы, симметричная матрица.


class eye : public Matrix {
public:
    eye(int a) {


        strk = a;
        stlb = a;

        for (int i = 0; i < a; i++) {
            for (int j = 0; j < a; j++) {
                if (i == j) {
                    values.push_back(1);
                }
                else {
                    values.push_back(0);
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
            for (int j = 0; j < strk; j++) {
                if (i == j) {
                    values.push_back(vec[i]);
                }
                else {
                    values.push_back(0);
                }
            }
        }
    };

    diagonal(Matrix mtr) {
        if (mtr.issquare() != true) {
            throw IndexError("Matrix: wrong size to make diagonal");
        }
        for (int i = 0; i < mtr.stroks(); i++) {
            for (int j = 0; j < mtr.stroks(); j++) {
                if (i == j) {
                    values.push_back(mtr.cell(i, j));
                }
                else {
                    values.push_back(0);
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
                    values.push_back(a);
                }
                else {
                    values.push_back(0);
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
                    values.push_back(a);
                }
                else {
                    values.push_back(0);
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
    

    void operator= (const int num) {
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
            for (int j = 0; j < c; j++) {
                if (j >= i) {
                    float a = (float)(rand() % 10000) / 1000;
                    values.push_back(a);
                }
                else {
                    values.push_back(this->print(j, i)); // всегда сработает так как матрица заполняется построчно.
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
След матрицы

Определитель матрицы(методом Гаусса)

Скалярное произведение векторов (уже готово)

Норма вектора(евклидова норма, максимальная норма)

Норма матрицы(норма Фробениуса)
*/










/*Дополнить библиотеку классов, разработанную в рамках ЛР1 и ЛР2. Реализовать операции над векторами и матрицами:

Угол между векторами

Ранг Матрицы (с помощью алгоритма Гаусса)

Обратная матрица (если существует)

Транспонирование матрицы*/






