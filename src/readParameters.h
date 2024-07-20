//
// Created by Xiao Shao on 2024/4/11.
//

void printHeader() {
    std::cout << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    std::cout << "|     ____ _                    ____  _           _  ___                      |\n";
    std::cout << "|    / ___| |__   ___ _ __ ___ |  _ \\| | __ _ ___| |/ (_)_ __                 |\n";
    std::cout << "|   | |   | '_ \\ / _ \\ '_ ` _ \\| |_) | |/ _` / __| ' /| | '_ \\                |\n";
    std::cout << "|   | |___| | | |  __/ | | | | |  __/| | (_| \\__ \\ . \\| | | | |               |\n";
    std::cout << "|    \\____|_| |_|\\___|_| |_| |_|_|   |_|\\__,_|___/_|\\_\\_|_| |_|               |\n";
    std::cout << "|                                                                             |\n";
    std::cout << "|   A Freeware for Unified Chemistry-Plasma Kinetics Simulation               |\n";
    std::cout << "|   License:      GNU LESSER GENERAL PUBLIC LICENSE, Version 2.1              |\n";
    std::cout << "|   Version:      0.1.0 (February 2024)                                       |\n";
    std::cout << "|   Author:       Xiao Shao                                                   |\n";
    std::cout << "|   Organization: King Abdullah University of Science and Technology (KAUST)  |\n";
    std::cout << "|   Contact:      xiao.shao@kaust.edu.sa                                      |\n";
    std::cout << "\\*---------------------------------------------------------------------------*/\n";
}

// Read parameters from case files
template <typename T>
T readParameter(const std::string& fileName, const std::string& paramName) {
    std::ifstream inFile(fileName);
    std::string line;

    if (!inFile.is_open()) {
        throw std::runtime_error("Unable to open file: " + fileName);
    }

    bool inBlock = false;
    std::optional<T> optValue; // Optional value for storing map

    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, ' ')) {
            iss >> std::ws;  // Skip whitespace
            if (key == paramName) {
                // For reading map (species to fraction)
                if constexpr (std::is_same<T, std::map<std::string, double>>::value) {
                    optValue.emplace(); // Initialize the map
                    inBlock = true;  // Entering the block
                    continue;
                }
                // For reading species list
                else if constexpr (std::is_same<T, std::vector<std::string>>::value) {
                    std::string valueStr;
                    std::getline(iss, valueStr, ';');
                    std::istringstream valueStream(valueStr);
                    T value;
                    std::string item;
                    while (std::getline(valueStream, item, ',')) {
                        item.erase(std::remove_if(item.begin(), item.end(), ::isspace), item.end());
                        value.push_back(item);
                    }
                    return value;
                }
                // For boolean type
                else if constexpr (std::is_same<T, bool>::value) {
                    std::string value;
                    if (std::getline(iss, value, ';')) {
                        if (value == "true" || value == "1" || value == "yes") {
                            return true;
                        } else if (value == "false" || value == "0" || value == "no") {
                            return false;
                        } else {
                            throw std::runtime_error("Invalid boolean format for parameter: " + paramName);
                        }
                    }
                }
                // For other non-boolean types
                else {
                    std::string valueStr;
                    if (std::getline(iss, valueStr, ';')) {
                        // Remove surrounding quotes
                        valueStr.erase(std::remove(valueStr.begin(), valueStr.end(), '\"'), valueStr.end());
                        // Replace placeholder <case> with the directory of fileName
                        size_t pos = valueStr.find("<case>");
                        if (pos != std::string::npos) {
                            std::filesystem::path dirPath = std::filesystem::path(fileName).parent_path();
                            valueStr.replace(pos, 6, dirPath.string());
                        }
                        std::istringstream valueStream(valueStr);
                        T value;
                        if (valueStream >> value) {
                            return value;
                        } else {
                            throw std::runtime_error("Invalid format for parameter: " + paramName);
                        }
                    }
                }
            }
        }

        if (inBlock && optValue.has_value()) {
            if (line.find('}') != std::string::npos) {
                break;  // Exiting the block
            }
            std::istringstream issItem(line);
            std::string species;
            double fraction;
            if (issItem >> species >> fraction) {
                if constexpr (std::is_same<T, std::map<std::string, double>>::value) {
                    (*optValue)[species] = fraction;
                }
            }
        }
    }
    if (optValue.has_value()) {
        return *optValue;
    }

    throw std::runtime_error("Parameter not found: " + paramName);
}

// Function to read CSV data
void readCSV(std::string fileName, std::vector<std::pair<double, double>> &y_t) {
    std::ifstream file(fileName);
    // Check if the file is open
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + fileName);
    }
    std::string line, value;

    while (getline(file, line)) {
        std::stringstream ss(line);
        double t_val, y_val;

        getline(ss, value, ',');
        t_val = stod(value);

        getline(ss, value, ',');
        y_val = stod(value);

        y_t.emplace_back(t_val, y_val);
    }
}

// Quadratic interpolation function
double interpolate(const std::vector<std::pair<double, double>> &y_t, double queryPoint) {
    if (y_t.size() < 3) {
        throw std::runtime_error("Need at least 3 data points for quadratic interpolation.");
    }

    for (size_t i = 0; i < y_t.size() - 2; ++i) {
        if (queryPoint >= y_t[i].first && queryPoint <= y_t[i + 2].first) {
            // Three points for quadratic interpolation
            double x0 = y_t[i].first, y0 = y_t[i].second;
            double x1 = y_t[i + 1].first, y1 = y_t[i + 1].second;
            double x2 = y_t[i + 2].first, y2 = y_t[i + 2].second;

            // Coefficients of the quadratic polynomial ax^2 + bx + c
            double a, b, c;

            double denom = (x0 - x1) * (x0 - x2) * (x1 - x2);
            a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom;
            b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom;
            c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / denom;

            // Interpolated value
            return a * queryPoint * queryPoint + b * queryPoint + c;
        }
    }

    std::cerr << "Warning: Query point '" << queryPoint << "' is out of the data range! ";
    double boundaryValue;
    if (queryPoint < y_t[0].first) {
        boundaryValue = y_t[0].second;
    } else {
        boundaryValue = y_t[y_t.size() - 1].second;
    }
    std::cerr << "Returning " << boundaryValue << std::endl;

    return boundaryValue;
}

