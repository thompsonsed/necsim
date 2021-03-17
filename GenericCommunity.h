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
        std::shared_ptr<T> community_ptr{nullptr};
    public:

        explicit GenericCommunity() noexcept : community_ptr(std::make_shared<T>())
        { }

        ~GenericCommunity() = default;

        GenericCommunity(GenericCommunity &&other) noexcept : GenericCommunity()
        {
            *this = std::move(other);
        }

        GenericCommunity(const GenericCommunity &other) noexcept : GenericCommunity()
        {
            this->community_ptr = other.community_ptr;
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
                this->community_ptr.swap(other.community_ptr);

            }
        }

        void applyNoOutput(std::shared_ptr<SpecSimParameters> sp)
        {
            community_ptr->applyNoOutput(sp);
        }

        void output()
        {
            community_ptr->output();
        }

        void speciateRemainingLineages(const string &filename)
        {
            community_ptr->speciateRemainingLineages(filename);
        }

    };
}
#endif //NECSIM_GENERICCOMMUNITY_H
