#include <stdlib.h>
#include <SFML\Graphics.hpp>
#include <SFML\Window.hpp>
#include <chrono>
#include <thread>
#include <iostream>
#include <math.h>
#include <vector>
#include <typeinfo>

#include "linalg.hpp"
#include "utility.hpp"
#include "physics.hpp"

//ToDo
//replace gamma correction
//window coordinate transform --> scrolling, integrate streamline/draw streamline (-y_dist)angle?
//vector2d insertion --> streamlines, panel method, velocity cacluclation --> add vertex operators +,-,*,/ and magnitude, angle
//split vec2d implementation to header and source file
//Body vertex interpolation fix, NAN occurs also for some panel geometries
//Linear Algebra Revision --> efficient Matrix multiplication (compare to EIGEN library)
//Implement algorithm that integrates influence of arbitrary panel distribution on other panels --> setup LSE such that total flow through body is minimized (not only controll poitns)
//SIMD operations, Loop unrolling
//Estimate source distribution based on panel angle, length and free flow --> accelerate convergence
//Load body from function description
//Gather draw / display / print function in class container UserIO
//draw field 1 colorbar

int main()
{   
    Settings::initialize("../settings.dat");
    
    //Timing variables for constant frame rate
    unsigned long int ms_start, ms_end;
    int ms_delay;
    const unsigned int ms_durr = round(1000.0 / Settings::frame_rate);

    //Initializes the pressure and velocity fields
    Scalar_Field pressure_field(Settings::pressure_field_samples_x, Settings::pressure_field_samples_y,0,0,600,400); 
    Vector_Field velocity_field(Settings::velocity_field_samples_x, Settings::velocity_field_samples_y,0,0,600,400);
    Uniform uniform(Settings::uniform_flow_x, Settings::uniform_flow_y);

    Input_Reader input_reader("../source_input.txt");
    input_reader.get_Input();

    if(Settings::use_body)
    {
        Input_Reader body_reader(Settings::body_file_path);
        body_reader.load_body_from_memory(); //Load body contour from point-file
        
        /* Load Body from specified function WIP
        std::vector<double> param;
        double steps = 200;
        for(unsigned int u = 0; u < steps; ++u)
        {   
            param.push_back(u * (2 * M_PI / steps));
        }
        vec2d (*funct_ptr)(double){ &Functions::ellipse_coord };
        body_reader.load_body_from_function(param, funct_ptr);
        */

        //Body::Body_List[0]->interpolate_vertices(); //WIP
        
        Body::Body_List[0]->set_offset(Settings::body_offset_x, Settings::body_offset_y); //Center Body in Window    
        Body::Body_List[0]->set_scale(Settings::body_scale_x, Settings::body_scale_y); //Scale Body
        Body::Body_List[0]->setup_source_panel(); //Setup LSE for source panel method
    }

    sf::RenderWindow window(sf::VideoMode(Settings::window_size_x, Settings::window_size_y), "Potential-Flow", sf::Style::Default, Settings::graphic_settings);
    

    while (window.isOpen() && Settings::max_frame-- != 0) //Set frame counter to negative number for infinite frames
    {
        { //Accounts partially for the delay of the program runtime
            using namespace std::chrono;
            ms_start = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        }

        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed) window.close();
            else if(event.type == sf::Event::MouseWheelMoved) UserIO::update_zoom(event.mouseWheel.delta);
        }

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) Settings::coord_offset_x -= 1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) Settings::coord_offset_x += 1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) Settings::coord_offset_y -= 1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) Settings::coord_offset_y += 1;

        //Update velocity and pressure field size depending on zoom and offset:
        vec2d tmp1 = vec2d(Settings::display_size_x, Settings::display_size_y);
        vec2d tmp2 = vec2d(Settings::coord_offset_x, Settings::coord_offset_y);

        tmp1.x = tmp1.x * Settings::mtr_per_pxl_x;
        tmp1.y = tmp1.y * Settings::mtr_per_pxl_y;

        (*Scalar_Field::Scalar_Field_List[0]).update_position(tmp1, tmp2);
        (*Vector_Field::Vector_Field_List[0]).update_position(tmp1, tmp2);
        
        //Clear Window for next frame
        window.clear();

        if(Settings::use_body){Body::Body_List[0]->calc_source_panel();} //Solve source panel distribution for body

        Physics::calc_velocity_field(*Vector_Field::Vector_Field_List[0]);
        Physics::calc_pressure_field(*Scalar_Field::Scalar_Field_List[0]);  //this is a relative pressure

        if(Settings::pressure_display_style == 0) (*Scalar_Field::Scalar_Field_List[0]).draw_field(window, 10, 0.5); //Alternative way to visualize Pressure field
        else (*Scalar_Field::Scalar_Field_List[0]).draw_field_2(window, 1);
        
        (*Vector_Field::Vector_Field_List[0]).draw_field(window, Settings::velocity_field_vector_scale, 1);

        Source::draw_sources(window);
        if(Settings::use_body){Body::Body_List[0]->draw_body(window);}

        //Calculates and draws stream lines
        vec2d tmp; tmp.x = 0;
        double y_coord;
        for(unsigned int u = 0; u < Settings::number_streamlines; ++u)
        {   
            tmp.y = ((Settings::display_size_y * u) / (Settings::number_streamlines - 1)) + Settings::display_offset_y;
            y_coord = Mapping::pxl_to_coord(tmp).y;

            Physics::draw_streamline(window, Physics::integrate_streamline(y_coord, Settings::streamline_step_size));
        }

        UserIO::print_GUI(window);

        //Change uniform inflow every frame
        dynamic_cast<Uniform*>(Source::Source_List[0])->change_flow(); //Change y-component of uniform flow

        Source::remove_sources(1); //Remove panel-sources for next frame recomputation

        window.display();

        { //Delay between frames accounts partially for the delay of the program runtime
            using namespace std::chrono;
            ms_end = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            ms_delay = ms_durr - (ms_end - ms_start);
            if(ms_delay < 0){ms_delay = 0;}
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(ms_delay));
    }

    return 0;
} 