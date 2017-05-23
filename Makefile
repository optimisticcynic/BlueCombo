#CFLAGS = -pthread -m64 -Wno-deprecated
CFLAGS = -lBlue
CFLAGS += $(shell root-config --cflags --libs)
CFLAGS += -LBlueForPhistar -g


all:ElectPlusMuon1D ElectPlusMuon2D

ElectPlusMuon1D: ElectPlusMuon1D.cc
	@g++ ${CFLAGS} -o ElectPlusMuon1D ElectPlusMuon1D.cc

ElectPlusMuon2D: ElectPlusMuon2D.cc
	@g++ ${CFLAGS} -o ElectPlusMuon2D ElectPlusMuon2D.cc