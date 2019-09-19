#!/bin/bash

awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' "$1" > "$1_tmp"
mv "$1_tmp" "$1"
