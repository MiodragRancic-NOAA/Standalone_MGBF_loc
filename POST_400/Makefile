FCMPI=mpiifort
FC=ifort

FFLAGS += -ftrapuv -fno-inline -no-ip -check all -check noarg_temp_created -fp-stack-check -fpe0 -ftz -O0

LIB=

CON= con_400.exe
TARGETS= $(CON)

all: $(TARGETS)

MOD= \
      def_kind.mod \

OBJ= \
      def_kind.o \
      CONVERT_400.o


RM= /bin/rm -f

$(TARGETS): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $(CON)

%.o:%.f90
	$(FC) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJ) $(MOD) $(CON)

# dependencies

CONVERT_AN.o: def_kind.o 


