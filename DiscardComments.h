#ifndef _INCL_DISCARD_COMMENTS
#define _INCL_DISCARD_COMMENTS
// C++ STDLIB INCLUDES
#include <fstream>
#include <string>
#include <iostream>
#include <sstream> 
namespace MY_CPP
{
// add static strings for read input file and check
static const    std::string STRING_CHECK_READ= "\t+++[Check_READ]:";
static const    std::string STRING_CHECK_IMPORTANT="\t!+++";
static const    std::string STRING_DIVIDE_LINE(50,'*');
// add static strings for flags on writing files and check
static const 	std::string STRING_WRITE_BEGIN= "[BEGIN]:	";
static const 	std::string STRING_WRITE_END= 	"[END]:	 SUCCESS:";
static const 	std::string STRING_CHECK_WRITE= "\t+++[Check_WRITE]:";
// add handler strings: warn error
static const 	std::string STRING_ERROR= "[*ERROR]:	";
static const 	std::string STRING_WARN= "[!WARN]:	";
// inline function to avoid multiple definition
inline std::string
discard_comments(
    const std::string& input_string)
	{
    // Create a copy of the input string, but without any text following a '!',
    // '#', or '%' character.
    std::string output_string = input_string;
    std::istringstream string_stream;

    // Discard any text following a '!' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '!');
    string_stream.clear();

    // Discard any text following a '#' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '#');
    string_stream.clear();

    // Discard any text following a '%' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '%');
    string_stream.clear();
    //>> comment, as we want to include directory/path
    /*
    // Discard any test following a '/' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '/');
    string_stream.clear();
    */
    return output_string;
	}// discard_comments
}

#endif // define the include_discard_comments
