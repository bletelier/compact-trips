SDSL_DIR = $$HOME
FLAGS = -std=c++11 -O3 -DNDEBUG -march=native
INCLUDES = -I$(SDSL_DIR)/include/
LIB = $(SDSL_DIR)/lib/libsdsl.a $(SDSL_DIR)/lib/libdivsufsort.a $(SDSL_DIR)/lib/libdivsufsort64.a 

all:
	g++ $(FLAGS) $(INCLUDES) create.cpp -o ../../build/nok3-v1_create $(LIB)
	g++ $(FLAGS) $(INCLUDES) ask_quantity.cpp -o ../../build/nok3-v1_ask_quantity $(LIB)
	g++ $(FLAGS) $(INCLUDES) ask_trips.cpp -o ../../build/nok3-v1_ask_trips $(LIB)
	g++ $(FLAGS) $(INCLUDES) ask_linestop.cpp -o ../../build/nok3-v1_ask_linestop $(LIB)
