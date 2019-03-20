// run by ./getter.out name start end
// name is name of reference
// start is starting point that you want to extract
// end is ending point
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char * argv[])
{
    string get_name = ">" + string(argv[1]);
    int get_start = stoi(argv[2]);
    int get_end = stoi(argv[3]);

    cout << get_name << endl << get_start << endl << get_end << endl;
    
    string name;
    while(cin >> name)
    {
        if(name == get_name)
        {
            string loaded;
            cin >> loaded;
            
            cout << loaded.substr(get_start, get_end) << "\n";
            break;
        }
        
        string garbage;
        cin >> garbage;
    }
    
    return 0;
}
