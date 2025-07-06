#!/bin/sh

input=$1

echo 'static const char *license_text = R"('
cat "$input"
echo ')";'
