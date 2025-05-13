#ifndef UTIL_H
#define UTIL_H

#include <vector>

std::vector<std::string> get_files_in_directory(const std::string& directory_path);
void substr(char *dst, char *src, int position, int length);
bool startsWith(const char *pre, const char *str);
std::string trim(const std::string& str);
void print_time();

#endif
