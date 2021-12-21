#pragma once

#include <stdlib.h>
#include <vector>

#include <SFML\Graphics.hpp>

//Simple static vector
template <typename T, unsigned int N>
struct stat_vector
{
    T elems[N];
    unsigned int size = N;

    //Indexing operator overloaded for quicker access to elements
    T operator[](const unsigned int& idx) const
    {
        return this.elems[idx];
    }

    //Overloading for write access
    T& operator[](const unsigned int& idx)
    {   
        return this->elems[idx];
    }

};

struct vec2d
{
    double x, y;

    double magnitude() const
    {
        return sqrt(this->x * this->x + this->y * this->y);
    }

    double magnitude_sqrd() const
    {
        return (this->x * this->x + this->y * this->y);
    }

    double angle() const
    {
        return atan2(this->y, this->x);
    }

    vec2d& operator+(const vec2d& rhs)
    {
        this->x += rhs.x;
        this->y += rhs.y;
        return *this;
    }

    vec2d add(const vec2d& rhs) const
    {   
        vec2d vec;
        vec.x = this->x + rhs.x;
        vec.y = this->y + rhs.y;
        return vec;
    }

    vec2d& operator-(const vec2d& rhs)
    {
        this->x -= rhs.x;
        this->y -= rhs.y;
        return *this;
    }

    vec2d subtract(const vec2d& rhs) const
    {   
        vec2d vec;
        vec.x = this->x - rhs.x;
        vec.y = this->y - rhs.y;
        return vec;
    }

    vec2d& operator*(const double scalar)
    {
        this->x *= scalar;
        this->y *= scalar;
        return *this;
    }

    vec2d multiply(const double scalar) const
    {   
        vec2d vec;
        vec.x = this->x * scalar;
        vec.y = this->y * scalar;
        return vec;
    }

    vec2d& operator/(const double scalar)
    {
        this->x /= scalar;
        this->y /= scalar;
        return *this;
    }

    vec2d divide(const double scalar) const
    {   
        vec2d vec;
        vec.x = this->x / scalar;
        vec.y = this->y / scalar;
        return vec;
    }

    vec2d& operator=(const vec2d& rhs)
    {
        this->x = rhs.x;
        this->y = rhs.y;
        return *this;
    }
};

struct vec3d
{
    double x, y, z;
};

class Matrix
{   
    public:
    unsigned int rows, columns;
    std::vector<double> elems;

    Matrix(unsigned int nrows, unsigned int ncols);
    //Matrix elements are stored in a row-major order
    unsigned int ind(unsigned int row, unsigned int col) const;

    void set_elem(unsigned int row, unsigned int col, double value);
    double get_elem(unsigned int row, unsigned int col) const;

    void transpose();
    Matrix get_transpose() const;
    void assemble(Matrix loc, std::vector<unsigned int> IDs); //works only for square matrices, check that or implement non square!
    void remove_dof(unsigned int dof);
    double get_frobenius() const;

    double get_max() const;
    double get_min() const;
    
    void overwrite(const Matrix& nmatrix);
    double scalar_prod(const Matrix& rhs) const;

    void print() const;
    void print_MATLAB() const;

    //Operate on object --> Fast
    Matrix& operator*(const Matrix& rhs);
    Matrix& operator*(const double& rhs);
    Matrix& operator-(const Matrix& rhs);
    Matrix& operator+(const Matrix& rhs);
    
    Matrix& operator=(const Matrix& rhs);

    //Operate on copy --> leave object const, are however slower
    Matrix multiply(const Matrix& rhs) const;
    Matrix multiply(const double& rhs) const;
    Matrix add(const Matrix& rhs) const;
    Matrix subtract(const Matrix& rhs) const;

    //SIMD operations
    void multiply(const Matrix& lhs, const Matrix& rhs, Matrix& res) const;

    static Matrix solve_LSE_Jacobi(const Matrix A, const Matrix b, Matrix* x0 = nullptr, double omega = 1.0, double epsilon = 1E-5, unsigned int max_count = 7500);
    //Usually faster than the Jacobi method
    static Matrix solve_LSE_GS(const Matrix A, const Matrix b, Matrix* x0 = nullptr, double omega = 1.5, double epsilon = 1E-5, unsigned int max_count = 7500);
    //WIP not correctly implemented
    static Matrix solve_LSE_Grad(const Matrix A, const Matrix b, double epsilon = 1E-5, unsigned int max_count = 7500);
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
