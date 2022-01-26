##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=lib_sv_rs_dopp
ConfigurationName      :=Debug
WorkspacePath          :=/home/kwq/codeLite_pro/miscWorkspace
ProjectPath            :=/home/kwq/codeLite_pro/miscWorkspace/lib_sv_rs_dopp
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=
Date                   :=25/01/22
CodeLitePath           :=/home/kwq/.codelite
LinkerName             :=gcc
SharedObjectLinkerName :=gcc -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.o.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName).so
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E 
ObjectsFileList        :="lib_sv_rs_dopp.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)inc 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := ar rcus
CXX      := gcc
CC       := gcc
CXXFLAGS :=  -g $(Preprocessors)
CFLAGS   :=  -g -fPIC $(Preprocessors)
ASFLAGS  := 
AS       := as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/src_myfunc.c$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(SharedObjectLinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)
	@$(MakeDirCommand) "/home/kwq/codeLite_pro/miscWorkspace/.build-debug"
	@echo rebuilt > "/home/kwq/codeLite_pro/miscWorkspace/.build-debug/lib_sv_rs_dopp"

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/src_myfunc.c$(ObjectSuffix): src/myfunc.c $(IntermediateDirectory)/src_myfunc.c$(DependSuffix)
	$(CC) $(SourceSwitch) "/home/kwq/codeLite_pro/miscWorkspace/lib_sv_rs_dopp/src/myfunc.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/src_myfunc.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/src_myfunc.c$(DependSuffix): src/myfunc.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/src_myfunc.c$(ObjectSuffix) -MF$(IntermediateDirectory)/src_myfunc.c$(DependSuffix) -MM src/myfunc.c

$(IntermediateDirectory)/src_myfunc.c$(PreprocessSuffix): src/myfunc.c
	$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/src_myfunc.c$(PreprocessSuffix) src/myfunc.c


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


