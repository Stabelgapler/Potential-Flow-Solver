#include <stdlib.h>
#include <SFML\Graphics.hpp>
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
//Avoid dynamic memory allocation: safe for individual functions preallocated memory
//replace gamma correction
//window coordinate transform
//vector2d insertion --> streamlines, panel method, velocity cacluclation --> add vertex operators +,-,*,/ and magnitude, angle
//split vec2d implementation to header and source file
//Body vertex interpolation fix
//Linear Algebra Revision --> efficient Matrix multiplication (compare to EIGEN library)
//Implement algorithm that integrates influence of arbitrary panel distribution on other panels --> setup LSE such that total flow through body is minimized (not only controll poitns)
//SIMD operations, Loop unrolling
//Estimate source distribution based on panel angle, length and free flow --> accelerate convergence

int main()
{   
    Settings::initialize("../settings.dat");

    //Timing variables
    unsigned long int ms_start, ms_end;
    int ms_delay;
    const unsigned int ms_durr = round(1000.0 / Settings::frame_rate);

    //Initializes the pressure and velocity fields
    Scalar_Field pressure_field(90,60,400,250,600,400); 
    Vector_Field velocity_field(20,15,400,250,600,400);
    Uniform uniform(Settings::uniform_flow_x, Settings::uniform_flow_y);

    Input_Reader input_reader("../source_input.txt");
    input_reader.get_Input();

    if(Settings::use_body)
    {
        Input_Reader body_reader(Settings::body_file_path);
        body_reader.load_body(); //Load body contour from point-file
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
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();

        if(Settings::use_body){Body::Body_List[0]->calc_source_panel();} //Solve source panel distribution for body

        Physics::calc_velocity_field(*Vector_Field::Vector_Field_List[0]);
        Physics::calc_pressure_field(*Scalar_Field::Scalar_Field_List[0]);  //this is a relative pressure

        //(*Scalar_Field::Scalar_Field_List[0]).draw_field(window, 10, 0.5); //Alternative way to visualize Pressure field
        (*Scalar_Field::Scalar_Field_List[0]).draw_field_2(window, 1);
        (*Vector_Field::Vector_Field_List[0]).draw_field(window, Settings::velocity_field_vector_scale, 1);

        Source::draw_sources(window);
        if(Settings::use_body){Body::Body_List[0]->draw_body(window);}

        for(unsigned int u = 0; u < 15; ++u) //Calculates and draws stream lines
        {
            Physics::draw_streamline(window, Physics::integrate_streamline(100, 50 + u*28.5714, 700, 5));
        }

        //Print angle of attack to screen
        Mapping::draw_angle_of_attack(window, 100, 30);

        dynamic_cast<Uniform*>(Source::Source_List[0])->change_flow(); //Change y-component of uniform flow

        Source::remove_sources(1); //Remove panel-sources for next frame recomputation

        window.display();

        { //Accounts partially for the delay of the program runtime
            using namespace std::chrono;
            ms_end = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
            ms_delay = ms_durr - (ms_end - ms_start);
            if(ms_delay < 0){ms_delay = 0;}
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(ms_delay));
    }

    return 0;
} 