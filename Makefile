# Compiler
CC = g++

# Default file
FILE = Data/equilatero_norm.txt

# work
#FILE = Data/example_3.txt # OK!
#FILE = Data/Lee_1.txt # OK! but a bit cheating bcs this one should just stay order-1
#FILE = Data/equilatero_norm.txt # OK, I obtain only regions of order 2!
    
# doesn't work (ordered based on which I want to fix first)
#FILE = Data/four_points.txt # Not ok, I obtain both regions of order 1 and 2, I only need regions of order 2!
#FILE = Data/hardest_one.txt
#FILE = Data/example_1.txt
#FILE = Data/example_2.txt

#FILE = Data/Lee_fig2.txt
#FILE = Data/Lee_fig2_2.txt
#FILE = Data/Lee_2.txt
#FILE = Data/degenerate_1.txt # I have no idea how it should even be fixed



# Generic flags that are always valid 
CFLAGS = -std=c++17 -Isrc/MyGAL/include
LIBS =  -lsfml-graphics -lsfml-window -lsfml-system

# The flags 
DEBUG_FLAGS = -g -fsanitize=address -fno-omit-frame-pointer
SRC = src/Main.cpp src/MinKnapsack.cpp src/Point2D.cpp
OUT = myProgram

# Targets
all: $(OUT)

$(OUT): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(OUT) $(LIBS)

debug_build: CFLAGS += $(DEBUG_FLAGS)
debug_build: $(OUT)

run: $(OUT)
	./$(OUT) --file $(FILE) --visualize 1

valgrind: debug_build
	valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all ./$(OUT) --file $(FILE) --visualize 0

debug: debug_build
	gdb --args ./myProgram --file $(FILE) --visualize 0
	
clean:
	rm -f $(OUT)