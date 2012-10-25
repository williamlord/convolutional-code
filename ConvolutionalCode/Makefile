# specify the target file below
EXEC := conv_test 
# specify any compile options, such as -g -Wall ...
CXXFLAGS := #-g #-Wall 
# any libs to be linked
LIBS := 

# You shouldn't need to change anything below this point
ITPPFLAGS := `itpp-config --cflags`
ITPPLIBS :=`itpp-config --libs`
#ITPPFLAGS := `it-config-3.8.1 --flags`
#ITPPLIBS :=`it-config-3.8.1 --libs`
CXXFLAGS += $(ITPPFLAGS)

CXX := g++
#CXX := g++-4.0
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.d)
CPPFLAGS += -MMD

.Phony: all deps objs clean rebuild

all: $(EXEC)

deps: $(DEPS)

objs: $(OBJS)
	
clean: 
	@rm -f *.o
	@rm -f *.d

rebuild: clean all
	
$(DEPS):
	@rm -f $(@:.d=.o); 


$(EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(addprefix -l, $(LIBS)) $(ITPPLIBS)

-include $(DEPS)
