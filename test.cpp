#include <bits/stdc++.h>
using namespace std;

int main(){
    bitset<256> bs;
    bs[3] = 1;
    bs[5] = 1;
    auto it = bs._Find_first();
    while (it < bs.size()){
        cout << it << " ";
        it = bs._Find_next(it);
    }
    cout << endl;
}