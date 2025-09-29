#include <bits/stdc++.h>
using namespace std;

int main(){
    unordered_map<bitset<10>, bitset<10>> m;
    bitset<10> b1(string("1010101010"));
    m[b1] = bitset<10>(string("1110001110"));
    cout << (m.find(0)==m.end()) << endl;
    m[0] |= bitset<10>(string("0000110001"));
    cout << m[b1] << endl;
    cout << m[0] << endl;
}