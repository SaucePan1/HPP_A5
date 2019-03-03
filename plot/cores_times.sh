#!/bin/bash

N=$1 
file=$2 
n_steps=$3
delta=$4
theta=$5
graphics=$6
compare_file=$7
script=$8
directory=$(pwd)
echo $directory 
echo "$N"
echo "$file"
echo "$delta"
echo "$theta"
echo "$graphics"
cp print_mean_std_complexity.py "../"
cp plot_cores_times.py "../"
cp makefile "../"
cd "../"&&pwd
echo "1 2 3 4 5 6 7 8" > number_of_cores.txt
make generate_cores_1 argument="$script"
printf '\n \n \n'

counter=1
while [ $counter -le 10 ]
	./one_core "$N" "$file" "$n_steps" "$delta" "$theta" "$graphics" 1
	((counter++))
done
printf '\n \n \n'
python3 print_mean_complexity.py >> data.txt
rm -f time.txt


make generate_cores_2 argument="$script"
printf '\n \n \n'

counter=1
while [ $counter -le 10 ]
	./two_core "$N" "$file" "$n_steps" "$delta" "$theta" "$graphics" 2
	((counter++))
done
printf '\n \n \n'
python3 print_mean_complexity.py >> data.txt
rm -f time.txt


make generate_cores_3 argument="$script"
printf '\n \n \n'

counter=1
while [ $counter -le 10 ]
	./third_core "$N" "$file" "$n_steps" "$delta" "$theta" "$graphics" 3
	((counter++))
done
printf '\n \n \n'
python3 print_mean_complexity.py >> data.txt
rm -f time.txt


make generate_cores_4 argument="$script"
printf '\n \n \n'

counter=1
while [ $counter -le 10 ]
	./fourth_core "$N" "$file" "$n_steps" "$delta" "$theta" "$graphics" 4
	((counter++))
done
printf '\n \n \n'
python3 print_mean_complexity.py >> data.txt
rm -f time.txt


make generate_cores_5 argument="$script"
printf '\n \n \n'

counter=1
while [ $counter -le 10 ]
	./fifth_core "$N" "$file" "$n_steps" "$delta" "$theta" "$graphics" 5
	((counter++))
done
printf '\n \n \n'
python3 print_mean_complexity.py >> data.txt
rm -f time.txt


make generate_cores_6 argument="$script"
printf '\n \n \n'

counter=1
while [ $counter -le 10 ]
	./six_core "$N" "$file" "$n_steps" "$delta" "$theta" "$graphics" 6
	((counter++))
done
printf '\n \n \n'
python3 print_mean_complexity.py >> data.txt
rm -f time.txt


make generate_cores_7 argument="$script"
printf '\n \n \n'

counter=1
while [ $counter -le 10 ]
	./seven_core "$N" "$file" "$n_steps" "$delta" "$theta" "$graphics" 7
	((counter++))
done
printf '\n \n \n'
python3 print_mean_complexity.py >> data.txt
rm -f time.txt


make generate_cores_8 argument="$script"
printf '\n \n \n'

counter=1
while [ $counter -le 10 ]
	./eigth_core "$N" "$file" "$n_steps" "$delta" "$theta" "$graphics" 8
	((counter++))
done
printf '\n \n \n'
python3 print_mean_complexity.py >> data.txt
rm -f time.txt

python3 -W ignore plot_cores_time.py

rm print_mean_std_complexity.py
rm plot_cores_time.py
rm data.txt
rm number_of_cores.txt
rm makefile
make clean
cd $directory











