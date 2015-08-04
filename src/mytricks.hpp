#ifndef MYTRICKS_H
#define MYTRICKS_H

#include <algorithm>
#include <functional>

namespace my{
  template<typename Dst_container_t, typename Src_container_t>
  Dst_container_t transform(const Src_container_t& _src, std::function<typename Dst_container_t::value_type (typename Src_container_t::const_reference)> F){
    Dst_container_t dst(_src.size()); // todo: to do tricks;
    std::transform(_src.begin(), _src.end(), dst.begin(), F);
    return std::move(dst);
  }
}



#endif /* MYTRICKS_H */
