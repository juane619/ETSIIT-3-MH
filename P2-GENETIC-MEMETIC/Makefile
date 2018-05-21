#FICHERO MAKEFILE CHETADO(GENERALIZADO) COMENTADO. COMPILACION SEPARADA.

#PHONY indica a make que no cree el archivo especificado.
.PHONY: clean tgz zip

#DIRECTORIOS
BIN= ../BIN
SRC= src
OBJ= obj
LIB= lib
INC= include
###

#NOMBRE EJECUTABLE
TARGET1= AGG
TARGET2= AGE
TARGET3= AM


# OPCIONES COMPILACION
CC = gcc
CFLAGS= -g -I$(INC)
Cxx= g++
#podemos aniadir $(FLAGSERRORES) a CXXFLAGS o no #-Wall -Wextra -pedantic
FLAGSERRORES= 
CXXFLAGS=  -std=c++11 -O2 
LIBS=
###

#INCLUDES Y EJECUTABLE
INCLUDES= -I$(INC)
EXE1= $(BIN)/$(TARGET1)
EXE2= $(BIN)/$(TARGET2)
EXE3= $(BIN)/$(TARGET3)
###


###SECCION OBJETOS###
#FICHEROS OBJETO
OBJECTS1= $(OBJ)/AGG.o $(OBJ)/timer.o $(OBJ)/random.o $(OBJ)/comun.o
OBJECTS2= $(OBJ)/AGE.o $(OBJ)/timer.o $(OBJ)/random.o $(OBJ)/comun.o
OBJECTS3= $(OBJ)/AM.o $(OBJ)/timer.o $(OBJ)/random.o $(OBJ)/comun.o
####FIN SECCION OBJETOS####


#DEPENDENCIAS AUTOMATICAS A PARTIR DE LOS .o
#(sustituimos .o por .dep de todos los objetos)
DEPS= $(OBJECTS:.o=.dep)

#REGLA POR DEFECTO
all: $(EXE1) $(EXE2) $(EXE3)

#Creamos ejecutable a partir de todos los .o
#($^ se refiere a TODAS las dependencias)
$(EXE1): 
	@echo "-- Linking $@ --"
	@$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC)/$(TARGET1).cpp $(SRC)/timer.cpp $(SRC)/random.cpp $(SRC)/comun.cpp  -o $@ $(LIBS)
	@echo "-- Link and compilation finished --"
	
$(EXE2): 
	@echo "-- Linking $@ --"
	@$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC)/$(TARGET2).cpp $(SRC)/timer.cpp $(SRC)/random.cpp $(SRC)/comun.cpp -o $@ $(LIBS)
	@echo "-- Link and compilation finished --"
	
$(EXE3): 
	@echo "-- Linking $@ --"
	@$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC)/$(TARGET3).cpp $(SRC)/timer.cpp $(SRC)/random.cpp $(SRC)/comun.cpp -o $@ $(LIBS)
	@echo "-- Link and compilation finished --"

#CREAMOS .o DE LOS ARCHIVOS QUE SI DEPENDEN DE UN .h (importante orden)
$(OBJ)/%.o: $(SRC)/%.cpp $(INC)/%.h
	clear
	@echo "\nCompiling $< ..."
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

#CREAMOS .o DE LOS ARCHIVOS QUE NO DEPENDEN DE UN .h (importante orden)
$(OBJ)/%.o: $(SRC)/%.cpp
	clear
	@echo "\nCompiling $< ..."
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<


#CREAMOS LOS ARCHIVOS CON LAS DEPENDENCIAS
#($< se refiere a la primera dependencia; $@ se refiere al objetivo)
$(OBJ)/%.dep: $(SRC)/%.cpp $(INC)/%.h
	@mkdir -p obj/ bin/	
	@$(CXX) -MM $(CXXFLAGS) $(INCLUDES) $< > $@

#INCLUIMOS DEPENDENCIAS
-include $(DEPS)

clean:
	@clear	
	@-rm -f $(OBJECTS) $(DEPS) *.tgz *.zip
	@-rm -f $(EXE)
	@echo "---------------------------------\nLimpiados *.o *.dep y binarios\n---------------------------------"	

#COMPRIMIDO Y LIMPIEZA RAR
tgz:	clean
	-tar -zcf $(TARGET).tgz *
	@echo "---------------------------------\nComprimido en tgz\n---------------------------------"


#COMPRIMIDO Y LIMPIEZA ZIP
zip:	clean
	-zip -sf- $(TARGET).zip -r *
	@echo "---------------------------------\nComprimido en ZIP\n---------------------------------"