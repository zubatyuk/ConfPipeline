#!/bin/bash

realpath ()                                                                                                                                                                                   
{                                                                                                                                                                                             
    f=$@;                                                                                                                                                                                     
    if [ -d "$f" ]; then                                                                                                                                                                      
        base="";                                                                                                                                                                              
        dir="$f";                                                                                                                                                                             
    else                                                                                                                                                                                      
        base="/$(basename "$f")";                                                                                                                                                             
        dir=$(dirname "$f");                                                                                                                                                                  
    fi;                                                                                                                                                                                       
    dir=$(cd "$dir" && /bin/pwd);                                                                                                                                                             
    echo "$dir$base"                                                                                                                                                                          
} 

if [[ ${BASH_SOURCE[@]} == $0 ]] ; then
    echo "Script is being run. Please source it instead."
else
    export PROJECT_ROOT=$(dirname $(realpath ${BASH_SOURCE[@]}))
    export PYTHONPATH="$PROJECT_ROOT:$PYTHONPATH"
    export PATH="$PROJECT_ROOT/scripts:$PATH"
fi

