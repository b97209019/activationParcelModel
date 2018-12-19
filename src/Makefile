CPP	= g++
FLAGS	:= -lm -fopenmp -fPIC -lgomp
FLAGS	+= -O3
#FLAGS	+= -g
FLAGS	+= -D PARALLEL
#FLAGS	+= -D __AS__
#FLAGS	+= -D __TCK__
FLAGS	+= -D $(FOO)
FLAGS	+= -I/usr/lib/gcc/x86_64-redhat-linux/5.3.1/include
OBJ	=	diffusion_growth_structure.o	\
		diffusion_growth_io_method.o	\
		diffusion_growth_calculation.o	\
		main.o			\

HD	=	micro_physics_const.h	\
		diffusion_growth_header_file.h


TOBJ	=	test.o

vpath %.h inc
vpath %.o obj
vpath %.cpp src
%.o: %.cpp 
	$(CPP) $(FLAGS) -c $< -o $@

main: $(OBJ) $(HD)
	$(CPP) $(FLAGS) $(filter %.o, $^) -o $@

diffusion_growth.so: $(OBJ) $(HD)
	$(CPP) $(FLAGS) $(filter %.o, $^) -o $@

test: test.o diffusion_growth_structure.o diffusion_growth_calculation.o
	$(CPP) $(FLAGS) $^ -o $@


useful_function.o : useful_function.h diffusion_growth_header_file.h
main.o: diffusion_growth_header_file.h

clean:
	rm -rf $(OBJ)
