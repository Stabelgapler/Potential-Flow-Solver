#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <math.h>

#include "utility.hpp"
#include "linalg.hpp"
#include "physics.hpp"



Input_Reader::Input_Reader(const std::string input_file_path_n): input_file_path(input_file_path_n), input_file(std::ifstream(input_file_path_n))
{
    return;
}

void Input_Reader::get_Input()
{
    if(input_file.is_open())
    {   
        std::string line;
        while (std::getline(input_file, line) )
        {   
            if(line.size() == 0 || line[0] == '#'){continue;} //Comments start with #
            else interpret_line(line);
        }
        input_file.close();
    }
    else throw std::invalid_argument("Could not open Input File: " + input_file_path);
    
    return;
}

void Input_Reader::interpret_line(std::string line)
{
    std::size_t index_s = line.find("Vector-Field:");
    if(index_s != std::string::npos)
    {   
        index_s += 13;
        std::vector<double> arg_vec;
        std::string arg_str;

        std::stringstream line_stream(line.substr(index_s));

        while (std::getline(line_stream, arg_str, ';'))
        {   
            arg_vec.push_back(atoi(arg_str.c_str()));
        }

        if(arg_vec.size() == 6)
        {
            new Vector_Field((unsigned int) arg_vec[0], (unsigned int) arg_vec[1], arg_vec[2], arg_vec[3], arg_vec[4], arg_vec[5]);
            //Input_Reader::print(arg_vec);
            return;
        }
        else throw std::invalid_argument("Could not interpret Input Command: " + line);
    }

    index_s = line.find("Uniform-Source:");
    if(index_s != std::string::npos)
    {   
        index_s += 15;
        std::vector<double> arg_vec;
        std::string arg_str;

        std::stringstream line_stream(line.substr(index_s));

        while (std::getline(line_stream, arg_str, ';'))
        {   
            arg_vec.push_back(atof(arg_str.c_str()));
        }

        if(arg_vec.size() == 2)
        {   
            Uniform* temp = new Uniform(arg_vec[0], arg_vec[1]);
            //temp->superpose_velocity(*Vector_Field::Vector_Field_List[(unsigned int) arg_vec[0]]);
            return;
        }
        else throw std::invalid_argument("Could not interpret Input Command: " + line);
    }

    index_s = line.find("Point-Source:");
    if(index_s != std::string::npos)
    {   
        index_s += 13;
        std::vector<double> arg_vec;
        std::string arg_str;

        std::stringstream line_stream(line.substr(index_s));

        while (std::getline(line_stream, arg_str, ';'))
        {   
            arg_vec.push_back(atof(arg_str.c_str()));
        }

        if(arg_vec.size() == 3)
        {   
            Point *temp = new Point(arg_vec[0], arg_vec[1], arg_vec[2]);
            //temp->superpose_velocity(*Vector_Field::Vector_Field_List[(unsigned int) arg_vec[0]]);
            return;
        }
        else throw std::invalid_argument("Could not interpret Input Command: " + line);
    }

    index_s = line.find("Vortex-Source:");
    if(index_s != std::string::npos)
    {   
        index_s += 14;
        std::vector<double> arg_vec;
        std::string arg_str;

        std::stringstream line_stream(line.substr(index_s));

        while (std::getline(line_stream, arg_str, ';'))
        {   
            arg_vec.push_back(atof(arg_str.c_str()));
        }

        if(arg_vec.size() == 3)
        {
            Vortex *temp = new Vortex(arg_vec[0], arg_vec[1], arg_vec[2]);
            //temp->superpose_velocity(*Vector_Field::Vector_Field_List[(unsigned int) arg_vec[0]]);
            return;
        }
        else throw std::invalid_argument("Could not interpret Input Command: " + line);
    }

    index_s = line.find("Doublet-Source:");
    if(index_s != std::string::npos)
    {   
        index_s += 16;
        std::vector<double> arg_vec;
        std::string arg_str;

        std::stringstream line_stream(line.substr(index_s));

        while (std::getline(line_stream, arg_str, ';'))
        {   
            arg_vec.push_back(atof(arg_str.c_str()));
        }

        if(arg_vec.size() == 3)
        {
            Doublet *temp = new Doublet(arg_vec[0], arg_vec[1], arg_vec[2]);
            //temp->superpose_velocity(*Vector_Field::Vector_Field_List[(unsigned int) arg_vec[0]]);
            return;
        }
        else throw std::invalid_argument("Could not interpret Input Command: " + line);
    }

    return;
}

void Input_Reader::print(std::vector<double> vec)
{
    for(unsigned int u = 0; u < vec.size(); ++u)
    {
        std::cout << "Argument " << u << ": " << vec[u] << "\n";
    }
    return;
}

void Input_Reader::load_body()
{   
    std::vector<double> vertices;
    unsigned int down_sample = Settings::body_points_down_sample;

    if(input_file.is_open())
    {   
        std::string line;
        std::string number_str;
        size_t char_ind_1 = 0;
        size_t char_ind_2 = 0;
        unsigned int line_ctr = down_sample - 1;

        while (std::getline(input_file, line))
        {   
            if(line.size() < 3 || line[0] == '#'){continue;} //Comments start with #
            else
            {   
                if(++line_ctr == down_sample)
                {
                    char_ind_1 = 0;
                    for ( ; char_ind_1 < line.length(); char_ind_1++ ){ if(isdigit(line[char_ind_1]) || line[char_ind_1] == '-') break; }
                    char_ind_2 = char_ind_1 + 1;
                    for ( ; char_ind_2 < line.length(); char_ind_2++ ){ if(isspace(line[char_ind_2]) || line[char_ind_2] == ';') break; }
                    
                    number_str = line.substr(char_ind_1, char_ind_2 - char_ind_1);
                    vertices.push_back(atof(number_str.c_str()));

                    char_ind_1 = char_ind_2 + 1;
                    for ( ; char_ind_1 < line.length(); char_ind_1++ ){ if(isdigit(line[char_ind_1]) || line[char_ind_1] == '-') break; }
                    char_ind_2 = char_ind_1 + 1;
                    for ( ; char_ind_2 < line.length(); char_ind_2++ ){ if(isspace(line[char_ind_2]) || line[char_ind_2] == ';') break; }

                    number_str = line.substr(char_ind_1, char_ind_2 - char_ind_1);
                    vertices.push_back(atof(number_str.c_str()));

                    line_ctr = 0;
                }
            }
        }
        input_file.close();
        Body* temp = new Body(vertices);
    }
    else throw std::invalid_argument("Could not open Input File: " + input_file_path);
    
    return;
}

void Input_Reader::get_string(std::string& str_ref, std::string key_str)
{   
    if(input_file.is_open())
    {   
        std::string line;
        while (std::getline(input_file, line) )
        {   
            if(line.size() == 0 || line[0] == '#'){continue;} //Comments start with #
            else if(line.find(key_str) != std::string::npos)
            {
                std::size_t index_s = line.find(key_str) + key_str.length() + 1;
                if(index_s >= line.length()){throw std::invalid_argument("Could not find input-argument for: " + key_str);}

                str_ref = line.substr(index_s);
                return;
            }
        }
        throw std::invalid_argument("Could not find input: " + key_str);
    }
    else throw std::invalid_argument("Could not open Input File: " + input_file_path);
}

void Input_Reader::get_double(double* double_ref, std::string key_str)
{
    if(input_file.is_open())
    {   
        std::string line;
        while (std::getline(input_file, line) )
        {   
            if(line.size() == 0 || line[0] == '#'){continue;} //Comments start with #
            else if(line.find(key_str) != std::string::npos)
            {
                std::size_t index_s = line.find(key_str) + key_str.length() + 1;
                if(index_s >= line.length()){throw std::invalid_argument("Could not find input-argument for: " + key_str);}
                
                *double_ref = atof(line.substr(index_s).c_str());
                return;
            }
        }
        throw std::invalid_argument("Could not find input: " + key_str);
    }
    else throw std::invalid_argument("Could not open Input File: " + input_file_path);
}

void Input_Reader::get_int(int* int_ref, std::string key_str)
{
    if(input_file.is_open())
    {   
        std::string line;
        while (std::getline(input_file, line) )
        {   
            if(line.size() == 0 || line[0] == '#'){continue;} //Comments start with #
            else if(line.find(key_str) != std::string::npos)
            {
                std::size_t index_s = line.find(key_str) + key_str.length() + 1;
                if(index_s >= line.length()){throw std::invalid_argument("Could not find input-argument for: " + key_str);}

                *int_ref = atoi(line.substr(index_s).c_str());
                return;
            }
        }
        throw std::invalid_argument("Could not find input: " + key_str);
    }
    else throw std::invalid_argument("Could not open Input File: " + input_file_path);
}

void Input_Reader::close_file()
{
    input_file.close();
    return;
}



double Mapping::linear(double from_low, double from_high, double to_low, double to_high, double value)
{
    return ((value - from_low) * (to_high - to_low) / (from_high - from_low)) + to_low;
}

double Mapping::linear_bound(double from_low, double from_high, double to_low, double to_high, double value)
{
    if(value <= from_low){return to_low;}
    if(value >= from_high){return to_high;}
    return linear(from_low, from_high, to_low, to_high, value);
}

sf::Color Mapping::colormap_rgb(const double v_min, const double v_max, const double value, const double alpha)
{
    const double v_avg = (v_max + v_min) / 2.0;

    unsigned int r_val, g_val, b_val;

    b_val = round(linear_bound(v_min, v_avg, 255, 0, value));

    if(value <= v_avg)
    {
        g_val = round(linear_bound(v_min, v_avg, 0, 255, value));
    }
    else
    {
        g_val = round(linear_bound(v_avg, v_max, 255, 0, value));
    }
    
    r_val = round(linear_bound(v_avg, v_max, 0, 255, value));

    sf::Color clr(r_val, g_val, b_val, alpha);

    return clr;
}

double Mapping::gamma_corr(double value, double gamma)
{   
    if(gamma == 1){return value;}
    if(value > 0){return pow(value, gamma);}
    else{return -pow(-value, gamma);}
}

void Mapping::draw_colorbar(sf::RenderWindow& window, double x_pos, double y_pos, double height, double width, double steps, double min_val, double max_val, double gamma, double disp_min, double disp_max)
{   
    double value = 0;
    double cell_height;
    cell_height = height/steps;

    sf::RectangleShape temp_l(sf::Vector2f(width, cell_height));
    temp_l.setOrigin(width/2, cell_height/2);

    for(unsigned int u = 0; u < steps; ++u)
    {   
        value = Mapping::linear(0, steps-1, min_val, max_val, u);            
        value = Mapping::gamma_corr(value, gamma);
        
        temp_l.setFillColor(Mapping::colormap_rgb(min_val, max_val, value, 255));
        temp_l.setPosition(x_pos, Mapping::linear(0, steps-1, y_pos + height/2, y_pos - height/2, u));
        
        window.draw(temp_l);
    }

    temp_l.setSize(sf::Vector2f(width * 2, 1));
    temp_l.setOrigin(width, 0.5);
    temp_l.setFillColor(sf::Color::White);

    temp_l.setPosition(x_pos, y_pos + height/2);
    window.draw(temp_l);

    temp_l.setPosition(x_pos, y_pos);
    window.draw(temp_l);

    temp_l.setPosition(x_pos, y_pos - height/2);
    window.draw(temp_l);

    sf::Text text;
    text.setFont(Settings::font);
    text.setCharacterSize(12);
    text.setFillColor(sf::Color::White);

    text.setString(std::to_string(min_val));
    text.setOrigin(sf::Vector2f(0, text.getLocalBounds().height * 0.75));
    text.setPosition(x_pos + 10, y_pos + height/2);
    window.draw(text);

    text.setString(std::to_string((min_val + max_val) / 2.0));
    text.setPosition(x_pos + 10, y_pos);
    window.draw(text);

    text.setString(std::to_string(max_val));
    text.setPosition(x_pos + 10, y_pos - height/2);
    window.draw(text);

}

double Mapping::coord_to_pxl()
{
    return 0;
}


//Variable initialization with default parameters:
sf::Font Settings::font = sf::Font();

double Settings::frame_rate = 30;
int Settings::window_size_x = 800;
int Settings::window_size_y = 500;

int Settings::use_body = 0;
std::string Settings::body_file_path = "";
int Settings::body_points_down_sample = 1;

double Settings::body_offset_x = 300;
double Settings::body_offset_y = 250;
double Settings::body_scale_x = 250;
double Settings::body_scale_y = 250;

double Settings::uniform_flow_x = 1;
double Settings::uniform_flow_y = 0;

double Settings::uniform_flow_max_y_val = uniform_flow_x / 12;
double Settings::uniform_flow_change_step = uniform_flow_max_y_val / 20;

//Loads settings from configuration file
void Settings::initialize(std::string file_path)
{
    Settings::font.loadFromFile("../src/Oswald-VariableFont_wght.ttf"); //Load text font

    Input_Reader settings_reader(file_path);

    settings_reader.get_double(&Settings::frame_rate, "frame-rate");
    settings_reader.get_int(&Settings::window_size_x, "window_size_x");
    settings_reader.get_int(&Settings::window_size_y, "window_size_y");

    settings_reader.get_int(&Settings::use_body, "use-body");
    settings_reader.get_string(Settings::body_file_path, "body_file_path");
    settings_reader.get_int(&Settings::body_points_down_sample, "body_points_down_sample");

    settings_reader.get_double(&Settings::body_offset_x, "body_offset_x");
    settings_reader.get_double(&Settings::body_offset_y, "body_offset_y");
    settings_reader.get_double(&Settings::body_scale_x, "body_scale_x");
    settings_reader.get_double(&Settings::body_scale_y, "body_scale_y");

    settings_reader.get_double(&Settings::uniform_flow_x, "uniform_flow_x");
    settings_reader.get_double(&Settings::uniform_flow_y, "uniform_flow_y");

    settings_reader.get_double(&Settings::uniform_flow_max_y_val, "uniform_flow_max_y_val");
    settings_reader.get_double(&Settings::uniform_flow_change_step, "uniform_flow_change_step");

    settings_reader.close_file();
}



//Debugging tools/variables
unsigned int Debug::ctr = 0;

void Debug::print_pos()
{
    std::cout << "Made it to counter: " << Debug::ctr << std::endl;
    Debug::ctr++;
}