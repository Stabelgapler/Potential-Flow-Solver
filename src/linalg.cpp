#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <stdexcept>

#include "immintrin.h" //SIMD AVX

#include "linalg.hpp"
#include "utility.hpp"

#include <SFML\Graphics.hpp>


Matrix::Matrix(unsigned int nrows, unsigned int ncolumns)
{
    this->rows = nrows;
    this->columns = ncolumns;

    for(unsigned int u = 0; u < rows * columns; ++u)
    {
        elems.push_back(0.0);
    }

    return;
}

unsigned int Matrix::ind(unsigned int row, unsigned int column) const
 {  
    if(row > rows || column > columns)
    {
        throw std::invalid_argument("Indice exceeds Matrix size!");
    }

    return (columns * (row - 1)) + (column - 1);
 }

void Matrix::set_elem(unsigned int row, unsigned int column, double value)
 {
    elems[ind(row, column)] = value;
    return;
 }

double Matrix::get_elem(unsigned int row, unsigned int column) const
{
   return elems[ind(row, column)];
}

void Matrix::transpose()
{
    std::vector<double> temp = elems;

    for(unsigned int u = 1; u <= rows; ++u)
    {
        for(unsigned int v = 1; v <= columns; ++v)
        {   
            if(u == v)
            {
                continue;
            }
            elems[ind(v,u)] = temp[ind(u,v)];
        } 
    }

    unsigned int row_temp = rows;
    rows = columns;
    columns = row_temp;

    return;
}

Matrix Matrix::get_transpose() const
{  
    Matrix T(this->columns, this->rows);

    for(unsigned int u = 1; u <= this->rows; ++u)
    {
        for(unsigned int v = 1; v <= this->columns; ++v)
        {   
            T.elems[ind(v,u)] = this->elems[ind(u,v)];
        } 
    }

    return T;
}

void Matrix::assemble(Matrix loc, std::vector<unsigned int> IDs)
{   
    for(unsigned int u = 0; u < IDs.size(); u++)
    {
        for(unsigned int v = 0; v < IDs.size(); v++)
        {   
            this->elems[this->ind(IDs[u] + 1, IDs[v] + 1)] += loc.elems[loc.ind(u+1,v+1)];
        }
    }

    return;
}

void Matrix::remove_dof(unsigned int dof)
{
    elems.erase(elems.begin() + (columns * (dof - 1)), elems.begin() + (columns * dof));
    rows = rows - 1;

    if(columns == 1){return;}

    for(unsigned int u = 0; u < rows; ++u)
    {   
        elems.erase(elems.begin() + ((columns * u) + (dof - 1)) - u);
    }

    columns = columns - 1;
    return;
}

double Matrix::get_frobenius() const
{   
    double norm = 0;
    for(unsigned int u = 0; u < elems.size(); ++u)
    {
        norm += elems[u] * elems[u];
    }

    return sqrt(norm);
}

double Matrix::get_max() const
{
    double max = elems[0];
    for(unsigned int u = 1; u < elems.size(); ++u)
    {
        if(elems[u] > max){max = elems[u];}
    }
    return max;
}

double Matrix::get_min() const
{
    double min = elems[0];
    for(unsigned int u = 1; u < elems.size(); ++u)
    {
        if(elems[u] < min){min = elems[u];}
    }
    return min;
}

double Matrix::scalar_prod(const Matrix& rhs) const
{
    if(this->elems.size() != rhs.elems.size())
    {
        throw std::invalid_argument("Cant compute scalar product, Dimensions do not agree!");
    }
    double prod = 0;

    for(unsigned int u = 0; u < this->elems.size(); ++u)
    {
        prod += this->elems[u] * rhs.elems[u];
    }

    return prod;
}

void Matrix::overwrite(const Matrix& nmatrix)
{
    this->rows = nmatrix.rows;
    this->columns = nmatrix.columns;
    this->elems = nmatrix.elems;

    return;
}

void Matrix::print() const
{   
    std::cout << std::endl;
    for(unsigned int u = 1; u <= rows; ++u)
    {   
        std::cout << "|";
        for(unsigned int v = 1; v <= columns; ++v)
        {   
            std::cout << elems[ind(u,v)];
            if(v < columns){std::cout << " \t";}
        }
        std::cout << "|" << std::endl;
    }
    return;
}

void Matrix::print_MATLAB() const
{   
    std::cout << std::endl;
    std::cout << "[";
    for(unsigned int u = 1; u <= rows; ++u)
    {   
        for(unsigned int v = 1; v <= columns; ++v)
        {   
            std::cout << elems[ind(u,v)];
            if(v < columns){std::cout << " ";}
        }
        std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;
    return;
}

Matrix& Matrix::operator*(const Matrix& rhs)
{   
    if(this->columns != rhs.rows)
    {
        throw std::invalid_argument("Matrix dimensions do not agree!");
    }

    Matrix product(this->rows, rhs.columns);

    for (int row = 1; row <= this->rows; row++)
    {
        for (int col = 1; col <= rhs.columns; col++)
        {
            for (int inner = 1; inner <= this->columns; inner++)
            {
                product.elems[product.ind(row, col)] += this->elems[this->ind(row, inner)] * rhs.elems[rhs.ind(inner, col)];
            }
        }
    }
    this->overwrite(product);
    return *this;
}

Matrix& Matrix::operator*(const double& rhs)
{
    for (unsigned int u = 0; u < this->elems.size(); ++u)
    {
        this->elems[u] = rhs * this->elems[u];
    }
    
    return *this;
}

Matrix& Matrix::operator+(const Matrix& rhs)
{
    if(this->rows != rhs.rows || this->columns != rhs.columns)
    {
        throw std::invalid_argument("Matrix dimensions do not agree!");
    }

    for(unsigned int u = 0; u < this->elems.size(); ++u)
    {
        this->elems[u] = this->elems[u] + rhs.elems[u];
    }

    return *this;
}

Matrix& Matrix::operator-(const Matrix& rhs)
{
    if(this->rows != rhs.rows || this->columns != rhs.columns)
    {
        throw std::invalid_argument("Matrix dimensions do not agree!");
    }

    for(unsigned int u = 0; u < this->elems.size(); ++u)
    {
        this->elems[u] = this->elems[u] - rhs.elems[u];
    }

    return *this;
}

Matrix& Matrix::operator=(const Matrix& rhs)
{   
    if(this->rows != rhs.rows || this->columns != rhs.columns)
    {
        throw std::invalid_argument("Matrix dimensions do not agree!");
    }

    this->elems = rhs.elems;
    
    return *this;
}

Matrix Matrix::multiply(const Matrix& rhs) const
{   
    if(this->columns != rhs.rows)
    {
        throw std::invalid_argument("Matrix dimensions do not agree!");
    }

    Matrix product(this->rows, rhs.columns);

    for (int row = 1; row <= this->rows; row++)
    {
        for (int col = 1; col <= rhs.columns; col++)
        {
            for (int inner = 1; inner <= this->columns; inner++)
            {
                product.elems[product.ind(row, col)] += this->elems[this->ind(row, inner)] * rhs.elems[rhs.ind(inner, col)];
            }
        }
    }
    
    return product;
}

Matrix Matrix::multiply(const double& rhs) const
{
    Matrix product(this->rows, this->columns);

    for (unsigned int u = 0; u < this->elems.size(); ++u)
    {
        product.elems[u] = rhs * this->elems[u];
    }
    
    return product;
}

Matrix Matrix::add(const Matrix& rhs) const
{
    if(this->rows != rhs.rows || this->columns != rhs.columns)
    {
        throw std::invalid_argument("Matrix dimensions do not agree!");
    }

    Matrix sum(this->rows, this->columns);

    for(unsigned int u = 0; u < this->elems.size(); ++u)
    {
        sum.elems[u] = this->elems[u] + rhs.elems[u];
    }

    return sum;
}

Matrix Matrix::subtract(const Matrix& rhs) const
{
    if(this->rows != rhs.rows || this->columns != rhs.columns)
    {
        throw std::invalid_argument("Matrix dimensions do not agree!");
    }

    Matrix difference(this->rows, this->columns);

    for(unsigned int u = 0; u < this->elems.size(); ++u)
    {
        difference.elems[u] = this->elems[u] - rhs.elems[u];
    }

    return difference;
}

//WIP
void Matrix::multiply(const Matrix& lhs, const Matrix& rhs, Matrix& res) const
{   
    //Pad data to be visible by 4
    //Load data from unaligned memory
    __m256d vec_lhs = _mm256_loadu_pd(&lhs.elems[0]); //n x m, row major order
    __m256d vec_rhs = _mm256_loadu_pd(&rhs.elems[0]); //m x k
    __m256d vec_res = _mm256_setzero_pd(); //n x k
    
    //_mm256_fmadd_pd(vec_lhs, vec_rhs, vec_res);
    vec_res = _mm256_add_pd(vec_lhs, vec_rhs);

    double* result = (double*)&vec_res;
    printf("%lf %lf %lf %lf\n", result[0], result[1], result[2], result[3]);
}

Matrix Matrix::solve_LSE_Jacobi(const Matrix A, const Matrix b, Matrix* x0, double omega, double epsilon, unsigned int max_count)
{
    if(A.rows != b.rows)
    {
        throw std::invalid_argument("Cant solve LSE, Dimensions do not agree!");
    }
    for(unsigned int u = 1; u <= A.rows; ++u)
    {
        if(A.elems[A.ind(u,u)] == 0)
        {
            throw std::invalid_argument("Cant solve LSE, Cant invert Diagonal!");
        }
    }

    bool flag = true;
    unsigned int count = 0;

    Matrix x(A.rows, 1);
    Matrix r(A.rows, 1);

    if(x0 != nullptr){x = *x0;}

    while(r.get_frobenius() > epsilon || flag)
    {   
        flag = false;

        r = b.subtract(A.multiply(x));
        for(unsigned int u = 0; u < A.rows; ++u)
        {
            x.elems[u] = x.elems[u] + omega * (r.elems[u] / A.elems[A.ind(u+1,u+1)]); 
        }

        count ++;
        if(count > max_count)
        {
            std::cout << "LSE Solution didnt converge, residual norm: " << r.get_frobenius() << std::endl;
            return x;
        }
    }

    std::cout << "Solved LSE in " << count << " steps." << std::endl;
    return x;
}

Matrix Matrix::solve_LSE_GS(const Matrix A, const Matrix b, Matrix* x0, double omega, double epsilon, unsigned int max_count)
{
    if(A.rows != b.rows)
    {
        throw std::invalid_argument("Cant solve LSE, Dimensions do not agree!");
    }
    for(unsigned int u = 1; u <= A.rows; ++u)
    {
        if(A.elems[A.ind(u,u)] == 0)
        {
            throw std::invalid_argument("Cant solve LSE, Cant invert Diagonal!");
        }
    }

    bool flag = true;
    unsigned int count = 0;

    double sum;

    Matrix x(A.rows, 1);
    Matrix r(A.rows, 1);

    if(x0 != nullptr){x = *x0;}
    
    while(r.get_frobenius() > epsilon || flag)
    {   
        flag = false;

        r = b.subtract(A.multiply(x));
        for(unsigned int u = 0; u < A.rows; ++u)
        {
            sum = 0;
        
            for(unsigned int v = 0; v < A.rows; ++v)
            {   
                if(u == v){continue;}
                sum += A.elems[A.ind(u+1,v+1)] * x.elems[v];
            }

            x.elems[u] = (1-omega) * x.elems[u] + (omega / A.elems[A.ind(u+1,u+1)]) * (b.elems[u] - sum); 
        }

        if(++count > max_count)
        {
            std::cout << "LSE Solution didnt converge, residual norm: " << r.get_frobenius() << std::endl;
            return x;
        }
    }

    std::cout << "Solved LSE in " << count << " steps." << std::endl;
    return x;
}

Matrix Matrix::solve_LSE_Grad(const Matrix A, const Matrix b, double epsilon, unsigned int max_count)
{
    if(A.columns != b.rows)
    {
        throw std::invalid_argument("Cant solve LSE, Dimensions do not agree!");
    }
    for(unsigned int u = 1; u <= A.rows; ++u)
    {   
        for(unsigned int v = 1; v <= A.columns; ++v)
        {
            if(A.elems[A.ind(u,v)] != A.elems[A.ind(v,u)])
            {
                throw std::invalid_argument("Cant solve LSE with Gradient-Method, Matrix not symmetric!");
            }
        }
    }

    unsigned int count = 0;
    double step_length = 0;

    Matrix x(A.rows, 1);
    Matrix r(A.rows, 1);

    r = b;

    do
    {
        count ++;
        
        step_length = (r.scalar_prod(r)) / (r.scalar_prod(A.multiply(r)));
        x = (x + (r * step_length));
        r = (b.subtract(A.multiply(x)));

        if(count > max_count)
        {
            std::cout << "LSE Solution didnt converge, residual norm: " << r.get_frobenius() << std::endl;
            return x;
        }

    } while (r.get_frobenius() > epsilon);

    return x;

}


std::vector<Vector_Field *> Vector_Field::Vector_Field_List;

Vector_Field::Vector_Field(unsigned int nnum_x, unsigned int nnum_y, double noffset_x, double noffset_y, double nsize_x, double nsize_y): x_component(Matrix(nnum_x, nnum_y)), y_component(Matrix(nnum_x, nnum_y)) 
{
    this->num_x = nnum_x;
    this->num_y = nnum_y;

    this->offset_x = noffset_x;
    this->offset_y = noffset_y;

    this->size_x = nsize_x;
    this->size_y = nsize_y;

    Vector_Field::Vector_Field_List.push_back(this);
}

void Vector_Field::set_entry(unsigned int x_num, unsigned int y_num, unsigned int comp, double value)
{
    if(comp == 1)
    {
        x_component.set_elem(x_num, y_num, value);
        return;
    }
    y_component.set_elem(x_num, y_num, value);
}
    
double Vector_Field::get_entry(unsigned int x_num, unsigned int y_num, unsigned int comp) const
{
    if(comp == 1)
    {
        return x_component.get_elem(x_num, y_num);
    }
    return y_component.get_elem(x_num, y_num);
}

void Vector_Field::add_entry(unsigned int x_num, unsigned int y_num, unsigned int comp, double value)
{
    if(comp == 1)
    {
        x_component.set_elem(x_num, y_num, x_component.get_elem(x_num, y_num) + value);
        return;
    }
    y_component.set_elem(x_num, y_num, y_component.get_elem(x_num, y_num) + value);
}

double Vector_Field::get_magnitude(unsigned int x_num, unsigned int y_num)
{   
    double magn = 0;
    double x_comp = this->x_component.get_elem(x_num, y_num);
    double y_comp = this->y_component.get_elem(x_num, y_num);

    magn = sqrt(x_comp * x_comp + y_comp * y_comp);

    return magn;
}

double Vector_Field::get_angle(unsigned int x_num, unsigned int y_num)
{
    double x_comp = this->x_component.get_elem(x_num, y_num);
    double y_comp = this->y_component.get_elem(x_num, y_num);

    return atan2(-y_comp, x_comp) * 180 / M_PI;
}

vec2d Vector_Field::get_entry_pos(unsigned int nnum_x, unsigned int nnum_y) const
{   
    vec2d position;

    position.x = (((nnum_x-1) / (double)(this->num_x-1)) * this->size_x) - (this->size_x / 2.0) + this->offset_x;
    position.y = (((nnum_y-1) / (double)(this->num_y-1)) * this->size_y) - (this->size_y / 2.0) + this->offset_y;

    return position;
}

void Vector_Field::update_position(const vec2d& size, const vec2d& offset)
{
    this->size_x = size.x;
    this->size_y = size.y;

    this->offset_x = offset.x;
    this->offset_y = offset.y;
}

void Vector_Field::draw_field(sf::RenderWindow& window, double scale, double gamma) 
{   
    double circle_radius = 1.5;
    sf::CircleShape temp_c(circle_radius);
    temp_c.setOrigin(circle_radius, circle_radius);
    temp_c.setFillColor(sf::Color::Red);

    double line_width = 0.75;
    sf::RectangleShape temp_l(sf::Vector2f(scale, line_width));
    temp_l.setOrigin(0, line_width/2);
    temp_l.setFillColor(sf::Color::White);

    vec2d coord, w_pxl;
    
    for(unsigned int u = 1; u <= this->num_x; ++u)
    {   
        for(unsigned int v = 1; v <= this->num_y; ++v)
        {   
            coord = get_entry_pos(u,v);
            w_pxl = Mapping::coord_to_pxl(coord);
            
            temp_l.setSize(sf::Vector2f(scale * pow(this->get_magnitude(u,v), gamma), line_width));
            temp_l.setRotation(this->get_angle(u,v));

            temp_c.setPosition(w_pxl.x, w_pxl.y);
            temp_l.setPosition(w_pxl.x, w_pxl.y);
            
            
            window.draw(temp_l);
            window.draw(temp_c);
        }
    }
}


std::vector<Scalar_Field *> Scalar_Field::Scalar_Field_List;

Scalar_Field::Scalar_Field(unsigned int nnum_x, unsigned int nnum_y, double noffset_x, double noffset_y, double nsize_x, double nsize_y): field(Matrix(nnum_x, nnum_y)) 
{
    this->num_x = nnum_x;
    this->num_y = nnum_y;

    this->offset_x = noffset_x;
    this->offset_y = noffset_y;

    this->size_x = nsize_x;
    this->size_y = nsize_y;

    Scalar_Field::Scalar_Field_List.push_back(this);
}

void Scalar_Field::set_entry(unsigned int x_num, unsigned int y_num, double value)
{
    field.set_elem(x_num, y_num, value);
}
    
double Scalar_Field::get_entry(unsigned int x_num, unsigned int y_num) const
{
    return field.get_elem(x_num, y_num);
}

void Scalar_Field::add_entry(unsigned int x_num, unsigned int y_num, double value)
{
    field.set_elem(x_num, y_num, field.get_elem(x_num, y_num) + value);
}

double Scalar_Field::get_min() const
{
    return this->field.get_min();
}

double Scalar_Field::get_max() const
{
    return this->field.get_max();
}

vec2d Scalar_Field::get_entry_pos(unsigned int nnum_x, unsigned int nnum_y) const
{   
    vec2d position;
    position.x = (((nnum_x-1) / (double)(this->num_x-1)) * this->size_x) - (this->size_x / 2.0) + this->offset_x;
    position.y = (((nnum_y-1) / (double)(this->num_y-1)) * this->size_y) - (this->size_y / 2.0) + this->offset_y;

    return position;
}

void Scalar_Field::update_position(const vec2d& size, const vec2d& offset)
{
    this->size_x = size.x;
    this->size_y = size.y;

    this->offset_x = offset.x;
    this->offset_y = offset.y;
}

void Scalar_Field::draw_field(sf::RenderWindow& window, double scale, double gamma) 
{   
    const double max_radius = 15;
    double value = 1.0;
    sf::CircleShape temp_c(value);

    vec2d coord, w_pxl;
    
    for(unsigned int u = 1; u <= this->num_x; ++u)
    {   
        for(unsigned int v = 1; v <= this->num_y; ++v)
        {  
            value = this->get_entry(u,v);

            if(value >= 0){temp_c.setFillColor(sf::Color::Red);}
            else
            {
                temp_c.setFillColor(sf::Color::Yellow);
                value = -value;
            }
            
            value = scale * pow(value, gamma);

            if(value > max_radius){continue;}

            coord = get_entry_pos(u,v);
            w_pxl = Mapping::coord_to_pxl(coord);
            
            temp_c.setRadius(value);
            temp_c.setOrigin(value, value);
            temp_c.setPosition(w_pxl.x, w_pxl.y);
            
            window.draw(temp_c);
        }
    }
}

void Scalar_Field::draw_field_2(sf::RenderWindow& window, double gamma) 
{   
    double min_val, max_val, cmin_val, cmax_val;
    
    min_val = this->get_min();
    max_val = this->get_max();

    cmin_val = Mapping::gamma_corr(min_val, gamma);
    cmax_val = Mapping::gamma_corr(max_val, gamma);

    double value = 1.0;
    vec2d coord, w_pxl;

    double cell_width, cell_height;
    cell_width = (this->size_x/(this->num_x - 1)) / Settings::mtr_per_pxl_x;
    cell_height = (this->size_y/(this->num_y - 1)) / Settings::mtr_per_pxl_y;

    sf::RectangleShape temp_l(sf::Vector2f(cell_width, cell_height));
    temp_l.setOrigin(cell_width/2, cell_height/2);

    for(unsigned int u = 1; u <= this->num_x; ++u)
    {   
        for(unsigned int v = 1; v <= this->num_y; ++v)
        {  
            value = this->get_entry(u,v);            
            value = Mapping::gamma_corr(value, gamma);
            
            temp_l.setFillColor(Mapping::colormap_rgb(min_val, max_val, value, 255));

            coord = get_entry_pos(u,v);
            w_pxl = Mapping::coord_to_pxl(coord);
            temp_l.setPosition(w_pxl.x, w_pxl.y);
            
            window.draw(temp_l);
        }
    }

    Mapping::draw_colorbar(window, 720, 250, 200, 6, 50, min_val, max_val, gamma, 0, 0);
}
