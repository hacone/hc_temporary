#!/bin/bash

cat basic.def.space common.def \
| sed -e "s/  */\t/g; s/^ *$//g" > basic.def
