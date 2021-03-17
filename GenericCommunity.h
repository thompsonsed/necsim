//
// Created by sam on 17/03/2021.
//

#ifndef NECSIM_GENERICCOMMUNITY_H
#define NECSIM_GENERICCOMMUNITY_H

#include <memory>
#include <string>
#include "SpecSimParameters.h"

namespace necsim
{

    template<class T>
    class GenericCommunity
    {
    private:
        std::shared_ptr<T> community{nullptr};
    public:

        GenericCommunity<T>() : community(std::make_shared<T>())
        { }

        ~GenericCommunity() = default;

        GenericCommunity(GenericCommunity &&other) noexcept : GenericCommunity()
        {
            *this = std::move(other);
        }

        GenericCommunity(const GenericCommunity &other) : GenericCommunity()
        {
            *this = other;
        };

        GenericCommunity &operator=(GenericCommunity other) noexcept
        {
            other.swap(*this);
            return *this;
        }

        void swap(GenericCommunity &other) noexcept
        {
            if(this != &other)
            {
                community.swap(other.community);

            }
        }

        void applyNoOutput(std::shared_ptr<SpecSimParameters> sp)
        {
            community->applyNoOutput(sp);
        }

        void output()
        {
            community->output();
        }

        void speciateRemainingLineages(const string &filename)
        {
            community->speciateRemainingLineages(filename);
        }

    };
}
#endif //NECSIM_GENERICCOMMUNITY_H
