#pragma once

template<typename T>
struct TraceInSpan
{
    // size_t span_index;
    std::vector<std::vector<VectorX<T>>> traces;
    std::vector<std::vector<T>> gradient_magnitude;
    std::vector<bool>is_loop;
    std::vector<std::vector<T>> traces_values;

    std::vector<std::vector<int>> point_types; //store the type of every point in the trace. 0: ridge, 1: valley, 2: pseudo ridge, 3: pseudo valley, 4: critical point
    std::vector<int> ridges; // record vertex_index in order
    std::vector<int> ridge_trace_index; // record corresponding trace_index of vertex_index in order
    std::vector<int> valleys;// record vertex_index in order
    std::vector<int> valley_trace_index; // record corresponding trace_index of vertex_index in order
    std::vector<int> pseudo_ridges; // record vertex_index in order
    std::vector<int> pseudo_ridge_trace_index; // record corresponding trace_index of vertex_index in order
    std::vector<int> pseudo_valleys;// record vertex_index in order
    std::vector<int> pseudo_valley_trace_index; // record corresponding trace_index of vertex_index in order

    std::vector<std::vector<T>> second_deriv_tangent; //second derivative in the tangent direction of isocontour

    // std::vector<std::array<std::shared_ptr<std::atomic<bool>>,2>> is_visited; //for connect_info,for every trace, store whether its two ends are visited. 0: not visited, 1: visited

    std::vector<std::array<int,8>> connect_info; //for every trace, store the connect information of two end points. for one end point, store [x,y,z, m]. x: 0,1 represents critical points or other trace. y: span index. z: critical point index in the span or trace index in the span. m: only for trace. 0,1 represents left or right end of the trace. So every trace has 8 int.

    std::vector<size_t> prefix_start; //store the start index of every trace in the prefix sum of the number of points in every trace

};