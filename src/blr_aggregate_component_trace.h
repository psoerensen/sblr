#ifndef SBLR_BLR_AGGREGATE_COMPONENT_TRACE_H
#define SBLR_BLR_AGGREGATE_COMPONENT_TRACE_H

#include <armadillo>
#include <stdexcept>

namespace sblr {
namespace core {

struct AggregateComponentTrace {
 arma::imat component_count;
 arma::imat realized_active_count;
 arma::imat stick_eligible_count;
 arma::imat stick_continue_count;
 arma::imat stick_stop_count;
};

inline void allocate_aggregate_component_trace(
 AggregateComponentTrace& trace, int draws, int components
) {
 if (draws<1 || components<2)
  throw std::invalid_argument("invalid aggregate component trace dimensions");
 trace.component_count.zeros(draws,components);
 trace.realized_active_count.zeros(draws,1);
 trace.stick_eligible_count.zeros(draws,components-1);
 trace.stick_continue_count.zeros(draws,components-1);
 trace.stick_stop_count.zeros(draws,components-1);
}

template <class ComponentState>
inline void capture_aggregate_component_trace(
 const ComponentState& component, int marker_count, int component_count,
 int draw, AggregateComponentTrace& trace
) {
 if (marker_count<0 || component_count<2 || draw<0 ||
     trace.component_count.n_rows<=static_cast<arma::uword>(draw) ||
     trace.component_count.n_cols!=static_cast<arma::uword>(component_count))
  throw std::invalid_argument("invalid aggregate component capture state");
 for (int marker=0;marker<marker_count;++marker) {
  const int value=static_cast<int>(component[static_cast<std::size_t>(marker)]);
  if (value<0 || value>=component_count)
   throw std::invalid_argument("component state is outside the aggregate trace domain");
  ++trace.component_count(static_cast<arma::uword>(draw),
                          static_cast<arma::uword>(value));
 }
 const int active=marker_count-trace.component_count(draw,0);
 trace.realized_active_count(draw,0)=active;
 int eligible=marker_count;
 for (int stick=0;stick<component_count-1;++stick) {
  const int stopped=trace.component_count(draw,stick);
  const int continued=eligible-stopped;
  trace.stick_eligible_count(draw,stick)=eligible;
  trace.stick_continue_count(draw,stick)=continued;
  trace.stick_stop_count(draw,stick)=stopped;
  eligible=continued;
 }
}

}  // namespace core
}  // namespace sblr

#endif
