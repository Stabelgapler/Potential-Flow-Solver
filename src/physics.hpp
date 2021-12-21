#pragma once

#include <stdlib.h>
#include <vector>

#include "linalg.hpp"


class Source
{
    public:
    static std::vector<Source *> Source_List;

    vec2d position;
    double intensity;

    virtual void calc_velocity(const vec2d& coord) const = 0;
    virtual void superpose_velocity(Vector_Field& field) const;

    static void remove_sources(unsigned int ind_1, unsigned int ind_2 = 0);

    static void draw_sources(sf::RenderWindow& window);
};

class Uniform: public Source
{
    public:
    double x_vel, y_vel;

    Uniform(double nx_vel, double ny_vel);

    double calc_angle();

    void calc_velocity(const vec2d& coord) const;
    void superpose_velocity(Vector_Field& field) const;

    void change_flow();
};

class Point: public Source
{
    public:
    Point(double nx_pos, double ny_pos, double nintensity);

    void calc_velocity(const vec2d& coord) const;
};

class Vortex: public Source
{
    public:
    Vortex(double nx_pos, double ny_pos, double nintensity);

    void calc_velocity(const vec2d& coord) const;
};

class Doublet: public Source
{
    public:
    Doublet(double nx_pos, double ny_pos, double nintensity);

    void calc_velocity(const vec2d& coord) const;
};

class Physics
{
    public:
    static void get_velocity(const vec2d& coord);

    static std::vector<vec2d> integrate_streamline(double x_start, double y_start, double x_end, double step);
    static void draw_streamline(sf::RenderWindow& window, std::vector<vec2d> pos_vec);

    static void calc_pressure_field(Scalar_Field& pfield, double rho = 1.0);

    static void calc_velocity_field(Vector_Field& vfield);
};

class Body
{   
    protected:
    std::vector<vec2d> vertices; //Vertices, not scaled
    std::vector<vec2d> panel_vertices; //Scaled vertices
    std::vector<vec2d> panel_mid_points;
    std::vector<double> panel_lengths;
    std::vector<double> panel_angles;

    double offset_x, offset_y, scale_x, scale_y;

    public:
    static std::vector<Body *> Body_List;
    double net_source_strength;

    //Panel Methods:
    Matrix LSE; //Linear System of equations
    Matrix SOL; //Last solution vector, can be used as start for next iteration

    Body(std::vector<vec2d> nvertices);

    void set_offset(double noffset_x, double noffset_y);
    void set_scale(double nscale_x, double nscale_y);

    void add_vertex(double x_pos, double y_pos);
    void interpolate_vertices(); //WIP

    void get_scaled_vertices(std::vector<vec2d>& scaled_vertices) const;

    //Panel Methods:
    void setup_source_panel(); //Determines LSE
    void calc_source_panel(); //Solves LSE

    void calc_vortex_panel();
    
    void draw_body(sf::RenderWindow& window);
};