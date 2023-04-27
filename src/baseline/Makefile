SDSL_DIR = $$HOME
INCLUDES = -I$(SDSL_DIR)/include/
LIB = $(SDSL_DIR)/lib/libsdsl.a $(SDSL_DIR)/lib/libdivsufsort.a $(SDSL_DIR)/lib/libdivsufsort64.a 
all:
	g++ -std=c++11 $(INCLUDES) create.cpp -o ../../build/baseline_create $(LIB)
	g++ -std=c++11 $(INCLUDES) ask_quantity.cpp -o ../../build/baseline_ask_quantity $(LIB)
	g++ -std=c++11 $(INCLUDES) ask_trips.cpp -o ../../build/baseline_ask_trips $(LIB)
	g++ -std=c++11 $(INCLUDES) ask_linestop.cpp -o ../../build/baseline_ask_linestop $(LIB)
