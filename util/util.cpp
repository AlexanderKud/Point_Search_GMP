#include <string.h>
#include <ctime>
#include <string>
#include <iostream>
#include <filesystem>
#include <vector>

#include "util.h"

using namespace std;
namespace fs = filesystem;

vector<string> get_files_in_directory(const string& directory_path) {
    vector<string> files;
    for (const auto& entry : fs::directory_iterator(directory_path)) {
        if (entry.is_regular_file()) {
            files.push_back(entry.path().filename().string());
        }
    }  
    return files;
}

void substr(char *dst, char *src, int position, int length) {
    int c = 0;
    while (c < length) {
      dst[c] = src[position+c-1];
      c++;
    }
    dst[c] = '\0';
}

bool startsWith(const char *pre, const char *str)
{
    size_t lenpre = strlen(pre),
           lenstr = strlen(str);
    return lenstr < lenpre ? false : memcmp(pre, str, lenpre) == 0;
}

std::string trim(const std::string& str) {
    auto start = str.begin();
    while (start != str.end() && std::isspace(*start)) ++start;
    auto end = str.end();
    do { --end; } while (end != start && std::isspace(*end));
    return std::string(start, end + 1);
}

void print_time() {
    time_t timestamp = time(NULL);
    struct tm datetime = *localtime(&timestamp);
    char output[50];
    strftime(output, 50, "%H:%M:%S", &datetime);
    cout << "[" << output << "] ";
}
