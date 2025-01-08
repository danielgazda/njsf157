# Debian linux
CPPC = g++-13
# EIGENLIB = /usr/include/eigen3/

# MacOS, latest version from Macports
# CPPC = g++-mp-14
# EIGENLIB = /opt/local/include/eigen3/

EIGENLIB = ./eigen/

CPPFLAGS = -std=c++23 -O2 -Wall

all: TISE.x TISE.min.x LSE.x HOdiag.x VPA.x

TISE.min.o: TISE.min.cpp numerov.hpp bisection.hpp constants.hpp potentials.hpp
	$(CPPC) -c $(CPPFLAGS) -I$(EIGENLIB) -I/Users/daniel/mfflect/optimlib/optim/header_only_version/ TISE.min.cpp

TISE.min.x: TISE.min.o
	$(CPPC) -o TISE.min.x $(CPPFLAGS) TISE.min.o

TISE.o: TISE.cpp numerov.hpp bisection.hpp constants.hpp potentials.hpp
	$(CPPC) -c $(CPPFLAGS) TISE.cpp

TISE.x: TISE.o
	$(CPPC) -o TISE.x $(CPPFLAGS) TISE.o

Gauss_Legendre.o: Gauss_Legendre.hpp Gauss_Legendre.cpp
	$(CPPC) -c $(CPPFLAGS) Gauss_Legendre.cpp

HO.o: HO.hpp HO.cpp
	$(CPPC) -c $(CPPFLAGS) HO.cpp

LSE.o: LSE.cpp potentials.hpp constants.hpp matrix.hpp
	$(CPPC) $(CPPFLAGS) -c -I $(EIGENLIB) LSE.cpp

LSE.x: Gauss_Legendre.o HO.o LSE.o
	$(CPPC) -o LSE.x $(CPPFLAGS) Gauss_Legendre.o HO.o LSE.o

HOdiag.o: HOdiag.cpp potentials.hpp constants.hpp matrix.hpp
	# $(CPPC) $(CPPFLAGS) -c -I $(EIGENLIB) -DMOMENTUMSPACE HOdiag.cpp
	$(CPPC) $(CPPFLAGS) -c -I $(EIGENLIB) HOdiag.cpp

HOdiag.x: HOdiag.o Gauss_Legendre.o HO.o
	$(CPPC) -o HOdiag.x Gauss_Legendre.o HO.o HOdiag.o

VPA.o: VPA.cpp RK.hpp constants.hpp potentials.hpp
	$(CPPC) -c $(CPPFLAGS) VPA.cpp

VPA.x: VPA.o
	$(CPPC) -o VPA.x VPA.o

clean:
	rm Gauss_Legendre.o  HO.o  HOdiag.o  LSE.o  VPA.o  TISE.o
	rm HOdiag.x  LSE.x  VPA.x  TISE.x

get_eigen:
	# Download Eigen version 3.4 (tested), not the current "master" branch
	test -d eigen || git clone -b 3.4 git@gitlab.com:libeigen/eigen.git

remove_eigen:
	test -d eigen && rm -rf eigen
