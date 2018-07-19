//
// Created by Dennis Goldfarb on 5/15/18.
//

#ifndef DEMIX_FRAGMENT_H
#define DEMIX_FRAGMENT_H


#include <string>

class Fragment {

public:
    Fragment() {} ;
    Fragment(int z, double monoMass, std::string name, std::string pep) : z(z), monoMass(monoMass), name(name), pep(pep)
    {
        monoMz = monoMass/z;
    };

    bool operator <(const Fragment & o) const
    {
        return monoMz < o.monoMz;
    }

    int z;
    double monoMass;
    double monoMz;
    std::string name;
    std::string pep;
};


#endif //DEMIX_FRAGMENT_H
