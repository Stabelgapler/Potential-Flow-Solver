#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "physics.hpp"
#include "linalg.hpp"
#include "utility.hpp"


std::vector<Source *> Source::Source_List;

void Source::remove_sources(unsigned int ind_1, unsigned int ind_2)
{   
    if(ind_2 == 0)
    {
        ind_2 = Source::Source_List.size();
    }
    if(ind_1 <= ind_2)
    {
        for(unsigned int u = ind_1; u < ind_2; ++u)
        {
            delete Source::Source_List[u];
        }
        Source::Source_List.erase(Source::Source_List.begin() + ind_1, Source::Source_List.begin() + ind_2);
    }
}

void Source::draw_sources(sf::RenderWindow& window)
{
    double circle_radius = 2;
    sf::CircleShape temp_c(circle_radius);
    temp_c.setOrigin(circle_radius, circle_radius);
    temp_c.setFillColor(sf::Color::Magenta);

    double x_pos, y_pos;

    for(unsigned int u=0; u < Source::Source_List.size(); ++u)
    {
        x_pos =Source::Source_List[u]->x_pos;
        y_pos =Source::Source_List[u]->y_pos;

        if(!isnan(x_pos) && !isnan(y_pos))
        {
            temp_c.setPosition(x_pos, y_pos);
            window.draw(temp_c);
        }
    }

    return;
}


Uniform::Uniform(double nx_vel, double ny_vel)
{
    Source::Source_List.push_back(this);

    this->x_pos = NAN;
    this->y_pos = NAN;

    this->intensity = NAN;

    this->x_vel = nx_vel;
    this->y_vel = ny_vel; 
}

double Uniform::calc_angle()
{
    return atan2(this->y_vel, this->x_vel);
}

void Uniform::calc_velocity(double pos[])
{
    PreAllocated::calc_velocity.x = this->x_vel;
    PreAllocated::calc_velocity.y = this->y_vel;
}

void Uniform::superpose_velocity(Vector_Field& field)
{
    for(unsigned int u = 1; u <= field.num_x; ++u)
    {   
        for(unsigned int v = 1; v <= field.num_y; ++v)
        {   
            field.add_entry(u,v,1,this->x_vel);
            field.add_entry(u,v,2,this->y_vel);
        }
    }

    return;
}

void Uniform::change_flow() //Changes Uniform flow y-velocity every frame
{
    static double de_vel = Settings::uniform_flow_change_step; //Velocity change increment

    this->y_vel += de_vel;

    if(abs(this->y_vel) >= Settings::uniform_flow_max_y_val){de_vel = -de_vel;}
}


Point::Point(double nx_pos, double ny_pos, double nintensity)
{
    Source::Source_List.push_back(this);

    this->x_pos = nx_pos;
    this->y_pos = ny_pos;

    this->intensity = nintensity;
}

void Point::calc_velocity(double pos[])
{
    double x_dist, y_dist;

    x_dist = pos[0] - this->x_pos;
    y_dist = pos[1] - this->y_pos;

    PreAllocated::calc_velocity.x = this->intensity * x_dist / (2 * M_PI * (x_dist * x_dist + y_dist * y_dist));
    PreAllocated::calc_velocity.y = -this->intensity * y_dist / (2 * M_PI * (x_dist * x_dist + y_dist * y_dist));
}

void Point::superpose_velocity(Vector_Field& field)
{   
    std::vector<double> pos;

    for(unsigned int u = 1; u <= field.num_x; ++u)
    {   
        for(unsigned int v = 1; v <= field.num_y; ++v)
        {   
            pos = field.get_entry_pos(u,v);
            
            calc_velocity(&pos[0]);

            field.add_entry(u,v,1,PreAllocated::calc_velocity.x);
            field.add_entry(u,v,2,PreAllocated::calc_velocity.y);
        }
    }

    return;
}


Vortex::Vortex(double nx_pos, double ny_pos, double nintensity)
{
    Source::Source_List.push_back(this);

    this->x_pos = nx_pos;
    this->y_pos = ny_pos;

    this->intensity = nintensity;
}

void Vortex::calc_velocity(double pos[])
{
    double x_dist, y_dist;

    x_dist = pos[0] - this->x_pos;
    y_dist = pos[1] - this->y_pos;

    PreAllocated::calc_velocity.x = -this->intensity * y_dist / (2 * M_PI * (x_dist * x_dist + y_dist * y_dist));
    PreAllocated::calc_velocity.y = -this->intensity * x_dist / (2 * M_PI * (x_dist * x_dist + y_dist * y_dist));
}

void Vortex::superpose_velocity(Vector_Field& field)
{   
    std::vector<double> pos;

    for(unsigned int u = 1; u <= field.num_x; ++u)
    {   
        for(unsigned int v = 1; v <= field.num_y; ++v)
        {   
            pos = field.get_entry_pos(u,v);

            calc_velocity(&pos[0]);

            field.add_entry(u,v,1,PreAllocated::calc_velocity.x);
            field.add_entry(u,v,2,PreAllocated::calc_velocity.y);
        }
    }

    return;
}


Doublet::Doublet(double nx_pos, double ny_pos, double nintensity)
{
    Source::Source_List.push_back(this);

    this->x_pos = nx_pos;
    this->y_pos = ny_pos;

    this->intensity = nintensity;
}

void Doublet::calc_velocity(double pos[])
{
    double x_dist, y_dist;

    x_dist = pos[0] - this->x_pos;
    y_dist = pos[1] - this->y_pos;

    PreAllocated::calc_velocity.x = this->intensity * (y_dist * y_dist - x_dist * x_dist) / (2 * M_PI * (x_dist * x_dist + y_dist * y_dist) * (x_dist * x_dist + y_dist * y_dist));
    PreAllocated::calc_velocity.y = 2 * this->intensity * x_dist * y_dist / (2 * M_PI * (x_dist * x_dist + y_dist * y_dist) * (x_dist * x_dist + y_dist * y_dist));
}

void Doublet::superpose_velocity(Vector_Field& field)
{   
    std::vector<double> pos;

    for(unsigned int u = 1; u <= field.num_x; ++u)
    {   
        for(unsigned int v = 1; v <= field.num_y; ++v)
        {   
            pos = field.get_entry_pos(u,v);
            
            calc_velocity(&pos[0]);

            field.add_entry(u,v,1,PreAllocated::calc_velocity.x);
            field.add_entry(u,v,2,PreAllocated::calc_velocity.y);
        }
    }

    return;
}


void Physics::get_velocity(double x_pos, double y_pos)
{
    PreAllocated::get_velocity.x = 0;
    PreAllocated::get_velocity.y = 0;

    double pos[2] = {x_pos, y_pos};

    for(unsigned int u=0; u < Source::Source_List.size(); ++u)
    {
        Source::Source_List[u]->calc_velocity(pos);
        PreAllocated::get_velocity.x += PreAllocated::calc_velocity.x;
        PreAllocated::get_velocity.y += PreAllocated::calc_velocity.y;
    }
}

std::vector<vec2d> Physics::integrate_streamline(double x_start, double y_start, double x_end, double step)
{
    std::vector<vec2d> pos_vec;

    vec2d pos;
    double scale;
    double vel_mag;

    const double eps = 1E-1;
    const unsigned long int max_its = 1E3;
    unsigned long int its = 0;

    pos.x = x_start;
    pos.y = y_start;

    pos_vec.push_back(pos);

    while(pos.x < x_end && pos.x >= x_start && its < max_its)
    {
        Physics::get_velocity(pos.x, pos.y);

        vel_mag = PreAllocated::get_velocity.magnitude();
        scale = step / vel_mag;

        if(vel_mag <= eps) //Stagnation point
        {
            return pos_vec;
        }

        pos.x += scale * PreAllocated::get_velocity.x;
        pos.y -= scale * PreAllocated::get_velocity.y;

        pos_vec.push_back(pos);

        ++its;
    }

    return pos_vec;
}

void Physics::draw_streamline(sf::RenderWindow& window, std::vector<vec2d> pos_vec)
{
    double x_dist, y_dist, dist, angle;
    
    double line_width = 1.5;
    sf::RectangleShape temp_l(sf::Vector2f(1, line_width));
    temp_l.setOrigin(0, line_width/2);
    temp_l.setFillColor(sf::Color::Blue);
    
    for(unsigned int u = 0; u < pos_vec.size()-1; ++u)
    {   
        x_dist = pos_vec[u+1].x - pos_vec[u].x;
        y_dist = pos_vec[u+1].y - pos_vec[u].y;

        dist = sqrt(x_dist * x_dist + y_dist * y_dist);
        angle = atan2(y_dist, x_dist) * 180 / M_PI;

        temp_l.setSize(sf::Vector2f(dist, line_width));
        temp_l.setRotation(angle);
        temp_l.setPosition(pos_vec[u].x, pos_vec[u].y);
        
        window.draw(temp_l);
    }

    return;
}

void Physics::calc_pressure_field(Scalar_Field& pfield, double rho)
{
    double v_x = dynamic_cast<Uniform*>(Source::Source_List[0])->x_vel;
    double v_y = dynamic_cast<Uniform*>(Source::Source_List[0])->y_vel;
    const double v_inf = v_x * v_x + v_y * v_y;

    double v_mag = 0;

    double p = 0;

    for(unsigned int u = 1; u <= pfield.num_x; ++u)
    {   
        for(unsigned int v = 1; v <= pfield.num_y; ++v)
        {   
            std::vector<double> pos = pfield.get_entry_pos(u,v);

            Physics::get_velocity(pos[0], pos[1]);
            v_mag = PreAllocated::get_velocity.magnitude_sqrd();
            
            p = 0.5 * rho * (v_inf - v_mag);

            pfield.set_entry(u,v,p);
        }
    } 

    return;
}

void Physics::calc_velocity_field(Vector_Field& vfield)
{   
    for(unsigned int u = 1; u <= vfield.num_x; ++u)
    {
        for(unsigned int v = 1; v <= vfield.num_y; ++v)
        {
            vfield.set_entry(u,v,1,0);
            vfield.set_entry(u,v,2,0);
        }
    }

    for(unsigned int u=0; u < Source::Source_List.size(); ++u)
    {
        Source::Source_List[u]->superpose_velocity(vfield);
    }
}


std::vector<Body *> Body::Body_List;

Body::Body(std::vector<vec2d> nvertices): SOL(nvertices.size(), 1), LSE(nvertices.size(), nvertices.size())
{   
    Body::Body_List.push_back(this);

    this->vertices = nvertices;

    this->offset_x = 0;
    this->offset_y = 0;

    this->scale_x = 1;
    this->scale_y = 1;
}

void Body::set_offset(double noffset_x, double noffset_y)
{
    this->offset_x = noffset_x;
    this->offset_y = noffset_y;
}

void Body::set_scale(double nscale_x, double nscale_y)
{
    this->scale_x = nscale_x;
    this->scale_y = nscale_y;
}
    
void Body::add_vertex(double x_pos, double y_pos)
{   
    vec2d vertex;
    vertex.x = x_pos;
    vertex.y = y_pos;

    this->vertices.push_back(vertex);
}

//WIP weird nan error probably in calc source panels, some entries in the linear system are nan
void Body::interpolate_vertices()
{   
    std::vector<vec2d> nvertices;
    unsigned int end_idx = this->vertices.size()-1;

    for(unsigned int u = 0; u < end_idx; ++u)
    {   
        nvertices.push_back(this->vertices[u]);
        nvertices.push_back(this->vertices[u].add(this->vertices[u+1]) / 2.0);
    }
    nvertices.push_back(this->vertices[end_idx]);
    nvertices.push_back(this->vertices[end_idx].add(this->vertices[0]) / 2.0);

    this->vertices = nvertices;

    Matrix nsol(nvertices.size(),1);
    this->SOL.overwrite(nsol); 
}

void Body::get_scaled_vertices(std::vector<vec2d>& scaled_vertices) const
{   
    vec2d temp_vertex;
    for(unsigned int u = 0; u < this->vertices.size(); ++u)
    {
        temp_vertex.x = this->vertices[u].x * this->scale_x;
        temp_vertex.y = this->vertices[u].y * this->scale_y;
        scaled_vertices.push_back(temp_vertex);
    }
    return;
}

//Sets up the Linear system of equation for the classic source panel method
void Body::setup_source_panel()
{   
    //Setup panels
    get_scaled_vertices(this->panel_vertices);
    unsigned int end_idx = this->panel_vertices.size()-1;

    for(unsigned int u = 0; u < end_idx; ++u)
    {   
        this->panel_mid_points.push_back(this->panel_vertices[u].add(this->panel_vertices[u+1]) / 2.0);
        this->panel_lengths.push_back(this->panel_vertices[u+1].subtract(this->panel_vertices[u]).magnitude());
        this->panel_angles.push_back(atan2((this->panel_vertices[u+1].y - this->panel_vertices[u].y), (this->panel_vertices[u+1].x - this->panel_vertices[u].x)));
    }

    this->panel_mid_points.push_back(this->panel_vertices[end_idx].add(this->panel_vertices[0]) / 2.0);
    this->panel_lengths.push_back(sqrt((this->panel_vertices[0].x - this->panel_vertices[end_idx].x) * (this->panel_vertices[0].x - this->panel_vertices[end_idx].x) + (this->panel_vertices[0].y - this->panel_vertices[end_idx].y) * (this->panel_vertices[0].y - this->panel_vertices[end_idx].y)));
    this->panel_angles.push_back(atan2((this->panel_vertices[0].y - this->panel_vertices[end_idx].y), (this->panel_vertices[0].x - this->panel_vertices[end_idx].x)));


    //Assemble linear system of equations
    double A, B, C, D, E, S, INTEGRAL;

    for(unsigned int u = 0; u < this->LSE.rows; ++u)
    {
        for(unsigned int v = 0; v < this->LSE.columns; ++v)
        {
            if(u == v){this->LSE.set_elem(u+1,v+1,M_PI);}
            else{
                A = -(this->panel_mid_points[u].x - this->panel_vertices[v].x) * cos(this->panel_angles[v]) - (this->panel_mid_points[u].y - this->panel_vertices[v].y) * sin(this->panel_angles[v]);
                B = (this->panel_mid_points[u].x - this->panel_vertices[v].x) * (this->panel_mid_points[u].x - this->panel_vertices[v].x) + (this->panel_mid_points[u].y - this->panel_vertices[v].y) * (this->panel_mid_points[u].y - this->panel_vertices[v].y);
                C = sin(this->panel_angles[u] - this->panel_angles[v]);
                D = -(this->panel_mid_points[u].x - this->panel_vertices[v].x) * sin(this->panel_angles[u]) + (this->panel_mid_points[u].y - this->panel_vertices[v].y) * cos(this->panel_angles[u]);
                E = sqrt(B - A*A);
                S = this->panel_lengths[v];
                INTEGRAL = (C/2) * log((S*S + 2*A*S + B) / B) + ((D - A*C) / E) * (atan((S + A) / E) - atan(A / E));
                
                this->LSE.set_elem(u+1,v+1,INTEGRAL);
            }
        }
    }

    //Invert solution if needed
    if(Settings::invert_solution)
    {
        this->LSE = this->LSE * -1;
    }

    return;
}

//Calculates the source panel method, i.e. source distribution on body surface to make its surface a streamline
void Body::calc_source_panel()
{   
    //Assemble right hand side of linear system 
    Matrix RHS(this->LSE.rows, 1);
    double Vx = dynamic_cast<Uniform*>(Source::Source_List[0])->x_vel;
    double Vy = dynamic_cast<Uniform*>(Source::Source_List[0])->y_vel;
    double V = sqrt(Vx * Vx + Vy * Vy);
    double Vangle = dynamic_cast<Uniform*>(Source::Source_List[0])->calc_angle();

    for(unsigned int u = 0; u < this->LSE.rows; ++u)
    {
        RHS.set_elem(u+1,1, 2 * M_PI * V * sin(this->panel_angles[u] - Vangle));
    }

    //Solve using Gauss Seidel iteration
    this->SOL = Matrix::solve_LSE_GS(this->LSE, RHS, &this->SOL, 1.8);

    //Distribute sources and calculate net source strength
    this->net_source_strength = 0;

    for(unsigned int u = 0; u < this->SOL.rows; ++u)
    {   
        this->net_source_strength += this->SOL.get_elem(u+1,1) * this->panel_lengths[u];
        Point* temp = new Point(this->panel_mid_points[u].x + this->offset_x, -this->panel_mid_points[u].y + this->offset_y, this->SOL.get_elem(u+1,1));
    }

    //std::cout << "Net-Panel-Source-Strength: " << this->net_source_strength << "\n";

    return;
}

/*
//Work in Progress, not yet working
void Body::calc_vortex_panel()
{   
    std::vector<double> s_points; //Body vertices scaled
    std::vector<double> c_point_x; //Points centered at panels, such that the c_points lie always between the s_points and vice versa
    std::vector<double> c_point_y;
    std::vector<double> length; //Length of edges formed by c_points
    std::vector<double> angle; //Angle of edges formed by c_points

    for(unsigned int u = 0; u < this->vertices.size(); ++u)
    {
        s_points.push_back(this->vertices[u].x * this->scale_x);
        s_points.push_back(this->vertices[u].y * this->scale_x);
    }

    for(unsigned int u = 0; u < s_points.size()-3; u = u+2)
    {
        c_point_x.push_back((s_points[u] + s_points[u+2]) / 2.0);
        c_point_y.push_back((s_points[u+1] + s_points[u+3]) / 2.0);
        length.push_back(sqrt((s_points[u+2] - s_points[u]) * (s_points[u+2] - s_points[u]) + (s_points[u+3] - s_points[u+1]) * (s_points[u+3] - s_points[u+1])));
        angle.push_back(atan2((s_points[u+3] - s_points[u+1]), (s_points[u+2] - s_points[u])));
    }

    c_point_x.push_back((s_points[s_points.size()-2] + s_points[0]) / 2.0);
    c_point_y.push_back((s_points[s_points.size()-1] + s_points[1]) / 2.0);
    length.push_back(sqrt((s_points[0] - s_points[s_points.size()-2]) * (s_points[0] - s_points[s_points.size()-2]) + (s_points[1] - s_points[s_points.size()-1]) * (s_points[1] - s_points[s_points.size()-1])));
    angle.push_back(atan2((s_points[1] - s_points[s_points.size()-1]), (s_points[0] - s_points[s_points.size()-2])));


    Matrix LSE(c_point_x.size(), c_point_x.size());
    Matrix RHS(c_point_x.size(), 1);
    double A, B, C, D, E, S, INTEGRAL;
    double Vx = dynamic_cast<Uniform*>(Source::Source_List[0])->x_vel; //Uniform velocities
    double Vy = dynamic_cast<Uniform*>(Source::Source_List[0])->y_vel;
    double V = sqrt(Vx * Vx + Vy * Vy);
    double Vangle = dynamic_cast<Uniform*>(Source::Source_List[0])->calc_angle();

    for(unsigned int u = 0; u < LSE.rows; ++u)
    {
        for(unsigned int v = 0; v < LSE.columns; ++v)
        {
            if(u == v){LSE.set_elem(u+1,v+1,M_PI);}
            else{
                A = -(c_point_x[u] - s_points[v*2]) * cos(angle[v]) - (c_point_y[u] - s_points[v*2 + 1]) * sin(angle[v]);
                B = (c_point_x[u] - s_points[v*2]) * (c_point_x[u] - s_points[v*2]) + (c_point_y[u] - s_points[v*2 + 1]) * (c_point_y[u] - s_points[v*2 + 1]);
                C = sin(angle[u] - angle[v]);
                D = -(c_point_x[u] - s_points[v*2]) * sin(angle[u]) + (c_point_y[u] - s_points[v*2 + 1]) * cos(angle[u]);
                E = sqrt(B - A*A);
                S = length[v];
                INTEGRAL = (C/2) * log((S*S + 2*A*S + B) / B) + ((D - A*C) / E) * (atan((S + A) / E) - atan(A / E));

                LSE.set_elem(u+1,v+1,INTEGRAL);
            }
        }
        RHS.set_elem(u+1,1, -2 * M_PI * V * sin(angle[u] - Vangle));
    }

    this->panel_sol = Matrix::solve_LSE_Jacobi(LSE, RHS, &this->panel_sol, 1);

    double net_source = 0;

    for(unsigned int u = 0; u < this->panel_sol.rows; ++u)
    {   
        net_source += this->panel_sol.get_elem(u+1,1);
        Point* temp = new Point(c_point_x[u] + this->offset_x, -c_point_y[u] + this->offset_y, this->panel_sol.get_elem(u+1,1));
    }

    std::cout << "Net-Vorticity: " << net_source << "\n";

    return;
}
*/

void Body::draw_body(sf::RenderWindow& window)
{   
    sf::ConvexShape polygon;
    polygon.setPointCount(this->vertices.size());
    polygon.setOutlineColor(sf::Color::Black);
    polygon.setOutlineThickness(1);
    
    for(unsigned int u = 0; u < this->vertices.size(); ++u)
    {   
        polygon.setPoint(u, sf::Vector2f(this->vertices[u].x * this->scale_x + this->offset_x, this->vertices[u].y * -this->scale_y + this->offset_y));
    }
    window.draw(polygon);
}