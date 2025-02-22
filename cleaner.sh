#!/bin/bash

index=$(ls job-* | wc -l)

while [ $index -gt -1 ];
do
        if [ $index -gt 50 ];
        then
                rm -r job-* 2> /dev/null;
        else
                index=$(ls job-* | wc -l);
                #echo $index;
        fi
        index=$(ls job-* | wc -l)
done