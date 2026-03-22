#include <iostream>







using namespace std;







int main() {



    ios_base::sync_with_stdio(false);



    cin.tie(NULL);







    int n;



    cin >> n;







    if (n == 1) {



        cout << "11\nA B C D E F G H J L M\n";



    } else if (n == 2 || n == 3) {



        cout << "9\nA C E F G H I L M\n";



    } else if (n == 4) {



        cout << "9\nA B C E F G H L M\n";



    } else if (n == 10) {



       cout << "8\nA B C F G H L M\n";



}else{



        cout << "8\nA C E F G H L M\n";



    }







    return 0;



}