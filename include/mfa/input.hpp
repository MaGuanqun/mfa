//--------------------------------------------------------------
// mfa input data structure
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
// 
// David Lenz
// Argonne National Laboratory
// dlenz@anl.gov
//--------------------------------------------------------------
#ifndef _INPUT_HPP
#define _INPUT_HPP

namespace mfa
{
    template <typename T>
    struct InputInfo
    {
        InputInfo(
                int     dom_dim_,
                int     pt_dim_,
                bool    structured_) :
            dom_dim(dom_dim_),
            pt_dim(pt_dim_),
            structured(structured_)
        { 
            // Assume unstructured data when ndom_pts_ is not passed to constructor
            if (structured)
                cerr << "ERROR: Conflicting constructor arguments for InputInfo" << endl;    
        }

        InputInfo(
                size_t      dom_dim_,
                size_t      pt_dim_,
                bool        structured_,
                VectorXi&   ndom_pts_) :
            dom_dim(dom_dim_),
            pt_dim(pt_dim_),
            structured(structured_),
            ndom_pts(ndom_pts_)
        {
            // Check that ndom_pts matches structured flag
            if ( (!structured && ndom_pts.size() != 0) ||
                (structured && ndom_pts.size() == 0)     ) 
            {
                cerr << "ERROR: Conflicting constructor arguments for InputInfo" << endl;    
                cerr << "  structured: " << boolalpha << structured << endl;
                cerr << "  ndom_pts: ";
                for (size_t k = 0; k < ndom_pts.size(); k++) cerr << ndom_pts(k) << " ";
                cerr << endl;
            }
        }

        void init()
        {
            if (is_initialized)
            {
                cerr << "Warning: Attempting to initialize a previously initialized InputInfo" << endl;
                return;
            }
            if (!validate())
            {
                cerr << "ERROR: Improper setup of InputInfo" << endl;
                exit(1);
            }
            else
            {
                // set total number of points
                tot_ndom_pts = domain.rows();
                
                // set parameters
                Param<T> temp_param(dom_dim, ndom_pts, domain, structured);
                swap(params, temp_param);
                
                // set grid data structure if needed
                if (structured)
                {
                    g.init(dom_dim, ndom_pts);
                }
            }

            is_initialized = true;
        }

        InputInfo(const InputInfo&) = delete;
        InputInfo(InputInfo&&) = delete;
        InputInfo& operator=(const InputInfo&) = delete;
        InputInfo& operator=(InputInfo&&) = delete;

        // InputInfo& operator=(InputInfo&& other)
        // {
        //     swap(*this, other);
        //     return *this;
        // }

        // friend void swap(InputInfo& first, InputInfo& second)
        // {
        //     swap(first.dom_dim, second.dom_dim);
        //     swap(first.pt_dim, second.pt_dim);
        //     swap(first.structured, second.structured);
        //     first.ndom_pts.swap(second.ndom_pts);
        //     swap(first.g, second.g);
        //     first.domain.swap(second.domain);
        //     swap(first.params, second.params);
        // }

        // Defined during construction
        int         dom_dim;
        int         pt_dim;
        bool        structured;
        VectorXi    ndom_pts;
        // VectorXi    model_dims;

        // Defined by user
        MatrixX<T>  domain;

        // Defined automatically during init()
        int             tot_ndom_pts{0};
        mfa::GridInfo   g;
        mfa::Param<T>   params;
        bool            is_initialized{false};

        

        class PtIterator
        {
            const bool  structured;
            size_t      lin_idx;
            mfa::VolIterator vol_it;
            const InputInfo&  info;

        public:
            PtIterator(const InputInfo& info_, size_t idx_) :
                structured(info_.structured),
                lin_idx(structured ? 0 : idx_),
                vol_it(structured ? VolIterator(info_.ndom_pts, idx_) : VolIterator()),
                info(info_)
            { }

            // prefix increment
            PtIterator operator++()
            {
                if(structured)
                    vol_it.incr_iter();
                else
                    lin_idx++;

                return *this;
            }


            bool operator!=(const PtIterator& other)
            {
                return structured ? (vol_it.cur_iter() != other.vol_it.cur_iter()) :
                                    (lin_idx != other.lin_idx);
            }

            bool operator==(const PtIterator& other)
            {
                return !(*this!=other);
            }

            void coords(VectorX<T>& coord_vec)
            {
                if(structured)
                    coord_vec = info.domain.row(vol_it.cur_iter());
                else
                    coord_vec = info.domain.row(lin_idx);
            }

            void params(VectorX<T>& param_vec)
            {
                if(structured)
                    param_vec = info.params.pt_params(vol_it);
                else
                    param_vec = info.params.pt_params(lin_idx);
            }

            void ijk(VectorXi& ijk_vec)
            {
                if (!structured)
                {
                    cerr << "ERROR: No ijk values in PtIterator for unstructured input" << endl;
                    exit(1);
                }

                ijk_vec = vol_it.idx_dim_;
            }

            int idx()
            {
                return structured ? vol_it.cur_iter() : lin_idx;
            }
        };

        PtIterator iterator(size_t idx) const
        {
            return PtIterator(*this, idx);
        }

        PtIterator begin() const
        {
            return PtIterator(*this, 0);
        }

        PtIterator end() const
        {
            return PtIterator(*this, tot_ndom_pts);
        }

        void pt_coords(size_t idx, VectorX<T>& coord_vec) const
        {
            coord_vec = domain.row(idx);
        }

        void pt_params(size_t idx, VectorX<T>& param_vec) const
        {
            if(structured)
            {
                VectorXi ijk(dom_dim);
                g.idx2ijk(idx, ijk);
                param_vec = params.pt_params(ijk);
            }
            else
            {
                param_vec = params.pt_params(idx);   
            }
        }

        // Test that user-provided data meets basic sanity checks
        bool validate() const
        { 
            bool is_valid =     (dom_dim > 0)
                            &&  (pt_dim > dom_dim)
                            &&  (pt_dim == domain.cols())
                            &&  (structured ? ndom_pts.size() == dom_dim : true)
                            &&  (structured ? ndom_pts.prod() == domain.rows() : true)
                            ;

            if (is_valid) return is_valid;
            else 
            {
                cerr << "InputInfo initialized with incompatible data" << endl;
                cerr << "  structured: " << boolalpha << structured << endl;
                cerr << "  dom_dim: " << dom_dim << ",  pt_dim: " << endl;
                cerr << "  ndom_pts: ";
                for (size_t k=0; k < ndom_pts.size(); k++) 
                    cerr << ndom_pts(k) << " ";
                cerr << endl;
                cerr << "  domain matrix dims: " << domain.rows() << " x " << domain.cols() << endl;

                return is_valid;
            }
        }
    };
}   // namespace mfa

#endif // _INPUT_HPP