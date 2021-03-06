#pragma once

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include <SFML\Graphics.hpp>

class Input_Reader
{
    protected:
    const std::string input_file_path;
    std::ifstream input_file;

    public:
    Input_Reader(const std::string input_file_path_n);
    
    void get_Input();
    void interpret_line(const std::string line);

    void load_body();

    void get_string(std::string& str_ref, std::string key_str);
    void get_double(double* double_ref, std::string key_str);
    void get_int(int* int_ref, std::string key_str);

    void close_file();

    static void print(std::vector<double> vec);
};

class Mapping
{
    public:
    static double linear(double from_low, double from_high, double to_low, double to_high, double value);
    static double linear_bound(double from_low, double from_high, double to_low, double to_high, double value);
    static sf::Color colormap_rgb(const double v_min, const double v_max, const double value, const double alpha = 255);
    static double gamma_corr(double value, double gamma);
    static void draw_colorbar(sf::RenderWindow& window, double x_pos, double y_pos, double height, double width, double steps, double min_val, double max_val, double gamma, double disp_min, double disp_max);
    static void draw_angle_of_attack(sf::RenderWindow& window, double x_pos, double y_pos, unsigned int size = 15);
    static double coord_to_pxl();
};

class Settings
{
    public:
    static sf::Font font; //Textfont
    static sf::ContextSettings graphic_settings; //Graphic settings for SFML Library

    static double frame_rate;
    static int window_size_x;
    static int window_size_y;
    
    static int use_body;
    static int invert_solution;
    static std::string body_file_path;
    static int body_points_down_sample; //Enables down-sampling, i.e. use only every nth point in body file

    static double body_offset_x;
    static double body_offset_y;
    static double body_scale_x;
    static double body_scale_y;

    static double uniform_flow_x;
    static double uniform_flow_y;

    static double uniform_flow_max_y_val;
    static double uniform_flow_change_step;


    static void initialize(std::string file_path);
};

class Debug
{
    public:
    static unsigned int ctr;

    static void print_pos();
};