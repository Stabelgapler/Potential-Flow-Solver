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
//performance
//replace gamma correction
//window coordinate transform


void change_flow(Uniform* flow) //Changes Uniform flow y-velocity every frame
{
    static double de_vel = Settings::uniform_flow_change_step; //Velocity change increment

    flow->y_vel += de_vel;

    if(abs(flow->y_vel) >= Settings::uniform_flow_max_y_val){de_vel = -de_vel;}
}

int main()
{   
    Settings::initialize("../settings.dat");

    //Timing variables
    unsigned long int ms_start, ms_end;
    int ms_delay;
    const unsigned int ms_durr = round(1000.0 / Settings::frame_rate);

    sf::ContextSettings graphic_settings;
    graphic_settings.antialiasingLevel = 8; //Graphics Setting for "smooth" edges

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
        Body::Body_List[0]->set_offset(Settings::body_offset_x, Settings::body_offset_y); //Center Body in Window    
        Body::Body_List[0]->set_scale(Settings::body_scale_x, Settings::body_scale_y); //Scale Body
    }

    sf::RenderWindow window(sf::VideoMode(Settings::window_size_x, Settings::window_size_y), "Potential-Flow", sf::Style::Default, graphic_settings);

    while (window.isOpen())
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

        Calculations::calc_velocity_field(*Vector_Field::Vector_Field_List[0]);
        Calculations::calc_pressure_field(*Scalar_Field::Scalar_Field_List[0], 1.0);

        //(*Scalar_Field::Scalar_Field_List[0]).draw_field(window, 10, 0.5); //Alternative way to visualize Pressure field
        (*Scalar_Field::Scalar_Field_List[0]).draw_field_2(window, 1);
        (*Vector_Field::Vector_Field_List[0]).draw_field(window, 8, 1);

        Source::draw_sources(window);
        if(Settings::use_body){Body::Body_List[0]->draw_body(window);}

        for(unsigned int u = 0; u < 21; ++u) //Calculates stream lines
        {
            Calculations::draw_streamline(window, Calculations::integrate_streamline(100, -50 + u*30, 700, 5));
        }

        change_flow(dynamic_cast<Uniform*>(Source::Source_List[0])); //Change uniform inflow

        Source::remove_sources(1); //Remove panel-source-solution for next frame

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