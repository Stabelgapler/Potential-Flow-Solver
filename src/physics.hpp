#pragma once

#include <stdlib.h>
#include <vector>

#include "linalg.hpp"


class Source
{
    public:
    static std::vector<Source *> Source_List;

    double x_pos, y_pos, intensity;

    virtual double* calc_velocity(double pos[]) = 0;
    virtual void superpose_velocity(Vector_Field& field) = 0;

    static void remove_sources(unsigned int ind_1, unsigned int ind_2 = 0);

    static void draw_sources(sf::RenderWindow& window);
};

class Uniform: public Source
{
    public:
    double x_vel, y_vel;

    Uniform(double nx_vel, double ny_vel);

    double calc_angle();

    double* calc_velocity(double pos[]);
    void superpose_velocity(Vector_Field& field);

    void change_flow();
};

class Point: public Source
{
    public:
    Point(double nx_pos, double ny_pos, double nintensity);

    double* calc_velocity(double pos[]);
    void superpose_velocity(Vector_Field& field);
};

class Vortex: public Source
{
    public:
    Vortex(double nx_pos, double ny_pos, double nintensity);

    double* calc_velocity(double pos[]);
    void superpose_velocity(Vector_Field& field);
};

class Doublet: public Source
{
    public:
    Doublet(double nx_pos, double ny_pos, double nintensity);

    double* calc_velocity(double pos[]);
    void superpose_velocity(Vector_Field& field);
};

class Physics
{
    public:
    static double* get_velocity(double x_pos, double y_pos);

    static std::vector<double> integrate_streamline(double x_start, double y_start, double x_end, double step);
    static void draw_streamline(sf::RenderWindow& window, std::vector<double> pos_vec);

    static void calc_pressure_field(Scalar_Field& pfield, double rho);

    static void calc_velocity_field(Vector_Field& vfield);
};

class Body
{   
    protected:
    std::vector<vec2d> vertices;
    double offset_x, offset_y, scale_x, scale_y;

    public:
    static std::vector<Body *> Body_List;

    Matrix panel_sol;

    Body(std::vector<vec2d> nvertices);

    void set_offset(double noffset_x, double noffset_y);
    void set_scale(double nscale_x, double nscale_y);

    void add_vertex(double x_pos, double y_pos);
    void get_scaled_vertices(std::vector<vec2d>& scaled_vertices) const;

    void calc_source_panel();
    void calc_vortex_panel();
    
    void draw_body(sf::RenderWindow& window);
};