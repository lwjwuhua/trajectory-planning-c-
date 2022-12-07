test:gpscovert.o Polynomial_insert.o Spline.o trajectoryslice_main.o
	g++ gpscovert.o Polynomial_insert.o Spline.o trajectoryslice_main.o -o test
gpscovert.o:gpscovert.cpp gpscovert.h
	g++ -Wall -O -g -c gpscovert.cpp -o gpscovert.o
Polynomial_insert.o:Polynomial_insert.cpp Polynomial_insert.h
	g++ -Wall -O -g -c Polynomial_insert.cpp -o Polynomial_insert.o
Spline.o:Spline.cpp Spline.h
	g++ -Wall -O -g -c Spline.cpp -o Spline.o

trajectoryslice_main.o:trajectoryslice_main.cpp 
	g++ -Wall -O -g -c trajectoryslice_main.cpp -o trajectoryslice_main.o
.PHONY:clean
clean:
	rm *.o test
