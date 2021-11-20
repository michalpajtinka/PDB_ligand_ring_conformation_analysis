#include "helper_functions.h"
#include <regex>


using namespace std;

string rstrip(const string s)
{
        return regex_replace(s, regex( "^\\s+" ), "");
}

string lstrip(const string s)
{
        return regex_replace(s, regex( "\\s+$" ), "");
}

string strip(const string s)
{
        return lstrip(rstrip(s));
}
