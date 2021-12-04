cd "src"
g++ -c main.cpp linalg.cpp utility.cpp physics.cpp -I "..\SFML-2.5.1\include"
g++ main.o linalg.o utility.o physics.o -o ../bin/main.exe -L "..\SFML-2.5.1\lib" -lsfml-graphics -lsfml-window -lsfml-system
pause