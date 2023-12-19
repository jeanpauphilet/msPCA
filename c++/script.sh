libeigen="/usr/local/lib/eigen-3.4.0"
echo "g++ -shared -fPIC -O2 -std=c++17 -I$libeigen -o algorithm2.so algorithm2.cpp"
g++ -shared -fPIC -O2 -std=c++17 -I$libeigen -o algorithm2.so algorithm2.cpp
