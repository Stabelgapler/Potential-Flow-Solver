#pragma once

#include <stdlib.h>
#include <vector>

#include <SFML\Graphics.hpp>

class Matrix
{   
    public:
    unsigned int rows, columns;
    std::vector<double> elems;

    Matrix(unsigned int nrows, unsigned int ncols);

    unsigned int ind(unsigned int row, unsigned int col) const;

    void set_elem(unsigned int row, unsigned int col, double value);
    double get_elem(unsigned int row, unsigned int col) const;

    double get_min_elem() const;
    double get_max_elem() const;

    void transpose();
    Matrix get_transpose() const;
    void assemble(Matrix loc, std::vector<unsigned int> IDs); //works only for square matrices, check that or implement non square!
    void remove_dof(unsigned int dof);
    double get_frobenius() const;
    void overwrite(const Matrix& nmatrix);
    double scalar_prod(const Matrix& rhs) const;

    void print() const;
    void print_MATLAB() const;

    Matrix operator*(const Matrix& rhs) const; //Should operate on object --> performance
    Matrix operator*(const double& rhs) const;
    Matrix operator-(const Matrix& rhs) const;
    Matrix operator+(const Matrix& rhs) const;
    Matrix& operator=(const Matrix& rhs);

    static Matrix solve_LGS_Jacobi(const Matrix A, const Matrix b, Matrix* x0 = nullptr, double omega = 1.0, double epsilon = 1E-5, unsigned int max_count = 7500);
    static Matrix solve_LGS_GS(const Matrix A, const Matrix b, Matrix* x0 = nullptr, double omega = 1.5, double epsilon = 1E-5, unsigned int max_count = 7500);
    static Matrix solve_LGS_Grad(const Matrix A, const Matrix b, double epsilon = 1E-5, unsigned int max_count = 7500);
};

class Vector_Field
{
    protected:
    Matrix x_component, y_component;

    public:
    static std::vector<Vector_Field *> Vector_Field_List;

    double offset_x, offset_y, size_x, size_y, num_x, num_y;

    Vector_Field(unsigned int nnum_x, unsigned int nnum_y, double noffset_x, double noffset_y, double nsize_x, double nsize_y);

    void set_entry(unsigned int x_num, unsigned int y_num, unsigned int comp, double value);
    double get_entry(unsigned int x_num, unsigned int y_num, unsigned int comp) const;

    void add_entry(unsigned int x_num, unsigned int y_num, unsigned int comp, double value);
    
    double get_magnitude(unsigned int x_num, unsigned int y_num);
    double get_angle(unsigned int x_num, unsigned int y_num);

    std::vector<double> get_entry_pos(unsigned int num_x, unsigned int num_y);

    void draw_field(sf::RenderWindow& window, double scale, double gamma);
};

class Scalar_Field
{
    protected:
    Matrix field;

    public:
    static std::vector<Scalar_Field *> Scalar_Field_List;

    double offset_x, offset_y, size_x, size_y, num_x, num_y;

    Scalar_Field(unsigned int nnum_x, unsigned int nnum_y, double noffset_x, double noffset_y, double nsize_x, double nsize_y);

    void set_entry(unsigned int x_num, unsigned int y_num, double value);
    double get_entry(unsigned int x_num, unsigned int y_num) const;

    void add_entry(unsigned int x_num, unsigned int y_num, double value);

    double get_min() const;
    double get_max() const;

    std::vector<double> get_entry_pos(unsigned int num_x, unsigned int num_y);

    void draw_field(sf::RenderWindow& window, double scale, double gamma);
    void draw_field_2(sf::RenderWindow& window, double gamma);
};
