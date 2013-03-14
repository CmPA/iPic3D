find . -maxdepth 3 -iname \*.cpp -exec indent -nut -nbad -nbap -bbo -nbc -br -brf -brs -cli2 -fc1 -fca -i2 -ci2 -ip2 -l1000 -lc1000 -lp -npcs -npsl -nprs -saf -sai -saw -sc -nsob -ss -ts2 {} \;
find . -maxdepth 3 -iname \*.h -exec indent -nut -nbad -nbap -bbo -nbc -br -brf -brs -cli2 -fc1 -fca -i2 -ci2 -ip2 -l1000 -lc1000 -lp -npcs -npsl -nprs -saf -sai -saw -sc -nsob -ss -ts2 {} \;
