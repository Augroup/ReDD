#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "read_buffer.hpp"
#include "signal_processor.hpp"
#include "seq.hpp"
#include "dtw.hpp"
#include "intervals.hpp"
#include "eventalign.hpp"
#include "aln.hpp"
#include "eventalign_redd.hpp"
namespace py = pybind11;
using namespace pybind11::literals;


//template<size_t K>
template<typename KmerType>
size_t pybind_model(py::module_ &m, std::string suffix) {
    using Model = PoreModel<KmerType>;
    Model::pybind(m, suffix);
    //BwaIndex<Model>::pybind_defs(m, suffix);//ref_index);
    BandedDTW<Model>::pybind_defs(m, suffix);
    GlobalDTW<Model>::pybind_defs(m, suffix);
    SignalProcessor<Model>::pybind(m, suffix);

    Sequence<Model>::pybind(m, suffix);
    Alignment<Model>::pybind(m, suffix);

    //m.def(("write_eventalign_"+suffix).c_str(), write_eventalign<PoreModel<K>>);
    auto fn = write_eventalign<Model>;
    m.def(("write_eventalign_"+suffix).c_str(), fn);

    auto fn_new = write_eventalign_new<Model>;
    m.def(("write_eventalign_new_"+suffix).c_str(), fn_new);
    auto fn_redd_new = write_eventalign_redd_new<Model>;
    m.def(("write_eventalign_redd_new_"+suffix).c_str(), fn_redd_new);
    // m.def("write_eventalign_redd_new", 
    //     &write_eventalign_redd_new<Model>,
    //     py::arg("alignment"),
    //     py::arg("write_name"),
    //     py::arg("signal_index"),
    //     py::arg("signal_np")
    // );
    return 0;
}

//template<size_t ...Ks>
//std::vector<size_t> pybind_kmers(py::module_ &m) {
//    return {(pybind_kmer<Ks>(m))...};
//}

PYBIND11_MODULE(_uncalled4, m) {
    m.doc() = R"pbdoc(Uncalled4: a toolkit for nanopore signal alignment, analysis, and visualization)pbdoc";

    py::class_<Config> config(m, "_Conf");
    Config::pybind_defs(config);

    RefCoord::pybind_defs(m);
    SeqRecord::pybind(m);

    ReadBuffer::pybind_defs(m);

    EventDetector::pybind_defs(m);

    py::class_<Normalizer> norm(m, "Normalizer");
    Normalizer::pybind_defs(norm);

    pybind_pore_model_params(m);

    py::bind_vector<std::vector<u8>>(m, "ArrayU8", py::buffer_protocol());
    py::bind_vector<std::vector<u16>>(m, "ArrayU16", py::buffer_protocol());
    // py::bind_vector<std::vector<u32>>(m, "ArrayU32", py::buffer_protocol());
    py::bind_vector<std::vector<u32>>(m, "ArrayU32", py::buffer_protocol());
    py::bind_vector<std::vector<int>>(m, "Arrayint", py::buffer_protocol());
    py::bind_vector<std::vector<float>>(m, "ArrayFloat", py::buffer_protocol());

    // .def(py::pickle(
    // [](const std::vector<u32>& v) {
    //     return py::cast(v); // returns py::object that is actually a list
    // },
    // [](py::object obj) {
    //     return obj.cast<std::vector<u32>>();
    // }
    // ));
    py::class_<redd_data_t>(m, "ReddData")
        .def(py::init<>())
        .def_readwrite("X", &redd_data_t::X)
        .def_readwrite("y_ref", &redd_data_t::y_ref)
        .def_readwrite("y_call", &redd_data_t::y_call)
        .def_readwrite("info", &redd_data_t::info)
        .def_readwrite("X_candidate", &redd_data_t::X_candidate)
        .def_readwrite("y_ref_candidate", &redd_data_t::y_ref_candidate)
        .def_readwrite("y_call_candidate", &redd_data_t::y_call_candidate)
        .def_readwrite("info_candidate", &redd_data_t::info_candidate)
        .def(py::pickle(
        // __getstate__: convert C++ object -> Python tuple
        [](const redd_data_t &self) {
            return py::make_tuple(
                self.X,
                self.y_ref,
                self.y_call,
                self.info,
                self.X_candidate,
                self.y_ref_candidate,
                self.y_call_candidate,
                self.info_candidate
            );
        },
        // __setstate__: convert tuple -> C++ object
        [](py::tuple t) {
            if (t.size() != 8)
                throw std::runtime_error("Invalid state for ReddData!");
            redd_data_t r;
            r.X = t[0].cast<std::vector<std::vector<std::vector<float>>>>();
            r.y_ref = t[1].cast<std::vector<std::vector<int>>>();
            r.y_call = t[2].cast<std::vector<std::vector<int>>>();
            r.info = t[3].cast<std::vector<std::string>>();
            r.X_candidate = t[4].cast<std::vector<std::vector<std::vector<float>>>>();
            r.y_ref_candidate = t[5].cast<std::vector<std::vector<int>>>();
            r.y_call_candidate = t[6].cast<std::vector<std::vector<int>>>();
            r.info_candidate = t[7].cast<std::vector<std::string>>();
            return r;
        }
    ));
    pybind_model<u16>(m, "U16");
    pybind_model<u32>(m, "U32");
    ModelDF::pybind(m);

    ProcessedRead::pybind(m);

    pybind_dtw(m);
    pybind_intervals(m);
    pybind_arrays(m);

    AlnDF::pybind(m);
    CmpDF::pybind(m);
}


