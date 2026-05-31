#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

struct Vec3
{
    double x, y, z;
};

struct Atom
{
    std::string element;    // name of atom type

    int typeIndex;          // atom type index

    double x, y, z;         // coords in world space
};

class POSCAR
{
public:
    Vec3 lattice[3];

    std::vector<Atom> atoms;

    // different atom types
    std::vector<std::string> elementTypes;


    // ----------------------------------------------
    // load poscar file function
    bool load(const std::string& filename)
    {
        std::ifstream file(filename);

        if (!file.is_open())
        {
            std::cerr << "Could not open POSCAR\n";
            return false;
        }

        std::string line;

        // Line 1: comment
        std::getline(file, line);

        // Line 2: scale factor
        std::getline(file, line);
        double scale = std::stod(line);

        // Lines 3-5: lattice vectors
        for (int i = 0; i < 3; i++)
        {
            std::getline(file, line);
            std::stringstream ss(line);

            ss >> lattice[i].x
                >> lattice[i].y
                >> lattice[i].z;

            lattice[i].x *= scale;
            lattice[i].y *= scale;
            lattice[i].z *= scale;
        }

        // Line 6: element names
        std::getline(file, line);

        std::stringstream elementStream(line);
        std::vector<std::string> elements;

        elementTypes.clear(); // clear vector of diff elements

        std::string element;

        while (elementStream >> element)
        {
            elements.push_back(element);
            elementTypes.push_back(element);
        }

        // Line 7: atom counts
        std::getline(file, line);

        std::stringstream countStream(line);
        std::vector<int> counts;

        int count;
        int totalAtoms = 0;

        while (countStream >> count)
        {
            counts.push_back(count);
            totalAtoms += count;
        }

        // Check for Selective Dynamics
        std::getline(file, line);

        bool selectiveDynamics = false;

        if (line[0] == 'S' || line[0] == 's')
        {
            selectiveDynamics = true;
            std::getline(file, line);
        }

        bool directCoords =
            (line[0] == 'D' || line[0] == 'd');

        atoms.clear();

        int elementIndex = 0;
        int atomsOfCurrentElement = 0;

        for (int i = 0; i < totalAtoms; i++)
        {
            if (atomsOfCurrentElement >= counts[elementIndex])
            {
                elementIndex++;
                atomsOfCurrentElement = 0;
            }

            do
            {
                std::getline(file, line);
            } while (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos);

            std::stringstream ss(line);

            double x, y, z;
            ss >> x >> y >> z;

            if (directCoords)
            {
                double cx =
                    x * lattice[0].x +
                    y * lattice[1].x +
                    z * lattice[2].x;

                double cy =
                    x * lattice[0].y +
                    y * lattice[1].y +
                    z * lattice[2].y;

                double cz =
                    x * lattice[0].z +
                    y * lattice[1].z +
                    z * lattice[2].z;

                x = cx;
                y = cy;
                z = cz;
            }

            atoms.push_back(
                {
                    elements[elementIndex],
                    elementIndex,
                    x,
                    y,
                    z
                });

            atomsOfCurrentElement++;
        }

        return true;
    }
};