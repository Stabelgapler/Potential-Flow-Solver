cd "src"
g++ -c -O3 -mavx -pg main.cpp linalg.cpp utility.cpp physics.cpp -I "..\SFML-2.5.1\include"
g++ -O3 -mavx -pg main.o linalg.o utility.o physics.o -o ../bin/main.exe -L "..\SFML-2.5.1\lib" -lsfml-graphics -lsfml-window -lsfml-system
pause

cd ..\bin
main.exe

gprof main.exe gmon.out > profile.txt
pause

profile.txt