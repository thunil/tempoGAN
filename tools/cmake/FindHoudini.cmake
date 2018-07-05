# In:
#  	HOUDINI_ROOT (alternatively, make sure 'houdini' is on $PATH)
#	A function FINDHOUDINIFOLDER is implemented to find houdini on windows
# 	default folders.
#
# Out:
#  HOUDINI_FOUND
#  HOUDINI_INCLUDE_DIRS
#  HOUDINI_LIBRARY_DIRS
#  HOUDINI_DEFINITIONS
#  HOUDINI_TOOLKIT_DIR
#

#Looks for folders (possibly inside Program Files/Side Effects Software/) 
#that have a folder starting with Houdini. This way we can target multiple 
#Houdini versions.
FUNCTION(FINDHOUDINIFOLDER curdir)
	message("Automatically Searching for Houdini at ${curdir}")
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  set(houdiniStrPos "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
		string(FIND ${child} "Houdini" houdiniStrPos)
		if(NOT houdiniStrPos STREQUAL "-1")
			message("Found Houdini ${curdir}/${child}")
			set(HOUDINI_ROOT "${curdir}/${child}")
			return()
		endif()
    ENDIF()
  ENDFOREACH()
  message("Automatic Houdini finding failed")
ENDFUNCTION()

if(NOT HOUDINI_ROOT)
	string(REPLACE "\\" "/" PROGRAMFILES_PATH $ENV{PROGRAMFILES})
	FINDHOUDINIFOLDER("${PROGRAMFILES_PATH}/Side Effects Software")
endif()

message("Houdini Root ${HOUDINI_ROOT}")
find_program(HOUDINI_BINARY houdini PATHS ${HOUDINI_ROOT}/bin)
if(HOUDINI_BINARY STREQUAL "HOUDINI_BINARY-NOTFOUND" )
    find_program(HOUDINI_BINARY houdini PATHS $ENV{HFS}/bin)
endif(HOUDINI_BINARY STREQUAL "HOUDINI_BINARY-NOTFOUND" )

if(HOUDINI_BINARY STREQUAL "HOUDINI_BINARY-NOTFOUND" )
    set(HOUDINI_FOUND FALSE)
	message("Houdini Binary Not Found")
else(HOUDINI_BINARY STREQUAL "HOUDINI_BINARY-NOTFOUND" )
	message("Houdini Binary Found")
    set(HOUDINI_FOUND TRUE)
    get_filename_component(HOUDINI_BASE ${HOUDINI_BINARY} PATH)
    get_filename_component(HOUDINI_BASE ${HOUDINI_BASE} PATH)
    set(HOUDINI_LIBRARY_DIRS "${HOUDINI_BASE}/custom/houdini/dsolib" CACHE PATH "")  
    set(HOUDINI_INCLUDE_DIRS "${HOUDINI_BASE}/toolkit/include" CACHE PATH "")
    set(HOUDINI_TOOLKIT_DIR "${HOUDINI_BASE}/toolkit")
    
    set(HOUDINI_DEFINITIONS "")

    # general definitions
    list(APPEND HOUDINI_DEFINITIONS -DVERSION="$ENV{HOUDINI_VERSION}")
    list(APPEND HOUDINI_DEFINITIONS -DMAKING_DSO )
    list(APPEND HOUDINI_DEFINITIONS -DUT_DSO_TAGINFO="vfxgal" )
    
    # platform-specific definitions
    if( APPLE )
        # OSX
    elseif( WIN32 )
        # Win32
        list(APPEND HOUDINI_DEFINITIONS -DI386)
        list(APPEND HOUDINI_DEFINITIONS -DWIN32)
        list(APPEND HOUDINI_DEFINITIONS -DSWAP_BITFIELDS)
        list(APPEND HOUDINI_DEFINITIONS -DDLLEXPORT=__declspec\(dllexport\) )
        list(APPEND HOUDINI_DEFINITIONS -DSESI_LITTLE_ENDIAN)
    else( UNIX )
        # Linux common
        list(APPEND HOUDINI_DEFINITIONS -DLINUX)
        list(APPEND HOUDINI_DEFINITIONS -D_GNU_SOURCE)
        list(APPEND HOUDINI_DEFINITIONS -DDLLEXPORT= )
        list(APPEND HOUDINI_DEFINITIONS -DSESI_LITTLE_ENDIAN)
        list(APPEND HOUDINI_DEFINITIONS -DENABLE_THREADS)
        list(APPEND HOUDINI_DEFINITIONS -DUSE_PTHREADS)
        list(APPEND HOUDINI_DEFINITIONS -DENABLE_UI_THREADS)
        list(APPEND HOUDINI_DEFINITIONS -DGCC3)
        list(APPEND HOUDINI_DEFINITIONS -DGCC4)
        list(APPEND HOUDINI_DEFINITIONS -Wno-deprecated)
        
        # Linux 64 bit
        list(APPEND HOUDINI_DEFINITIONS -DAMD64)
        list(APPEND HOUDINI_DEFINITIONS -DSIZEOF_VOID_P=8)
    endif( UNIX )

endif(HOUDINI_BINARY STREQUAL "HOUDINI_BINARY-NOTFOUND" )


if( Houdini_FIND_REQUIRED AND NOT HOUDINI_FOUND )
    message(FATAL_ERROR "Could not find houdini")
endif( Houdini_FIND_REQUIRED AND NOT HOUDINI_FOUND )
