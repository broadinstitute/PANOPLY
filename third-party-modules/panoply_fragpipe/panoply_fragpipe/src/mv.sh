#!/bin/sh
# sh function to ensure that mv is available regardless of $PATH
mv()
{
    /bin/mv $@
}

typeset -fx mv
