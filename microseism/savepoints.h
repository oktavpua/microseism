#ifndef SAVEPOINTS_H
#define SAVEPOINTS_H

/**
 * \file   shock.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  Информация об ударах.
 */

#include <vector>
#include <array>
#include <string>

void read_save_points(
        const std::string &fname,
        std::vector<int> &ids,
        std::vector< std::array<int, 3> > &coords
        );

#endif
