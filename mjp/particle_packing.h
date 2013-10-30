#if !defined(__particle_packing__h)
#define __particle_packing__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/particle_packing.h,v $

/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2010  Kort Travis                                         */
/*                                                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


// #define _debug_print_
#undef _debug_print_

#define state_partner_lists

// #define event_history
#undef event_history

/** 
 * @file particle_packing.h
 * @brief Classes and methods to implement event-driven dynamics simulation of maximally-jammed particle packing (MJP).
 * Based on pseudo-code and discussion in the reference:
 *   Boris D. Lubachevsky and Frank H. Stillinger,
 *   "Geometric Properties of Random Disk Packings",
 *   Journal of Statistical Physics, Vol. 60, Nos. 5/6, Pp. 561-583, (1990)   
 * This namespace is implemented such that it can be easily made independent of namespace ``TMatrix'', if so desired.
 */
 
namespace particle_packing{

// -------------- allow this module to be independent of the TMatrix namespace: ------------------------
#if !defined(__numericalConstants__h) && 0

template <class T> 
class numberTraits{
public:

#if 1
  typedef T complexType;   // complex type if T is real.
  typedef T magnitudeType; // magnitude type if T is complex.
  typedef T integralType;  // integer type of corresponding range.
#endif

};

// explicit specialization for double std::complex<double>, and std::complex<T>
template <> 
class numberTraits<double>{
public:
  typedef std::complex<double> complexType;   // complex type if T is real.  
  typedef double magnitudeType; // magnitude type if T is complex.
  typedef long integralType;  // integer type of corresponding range.

};

#if 1
template <> 
template <class T>
class numberTraits< std::complex<T> >{
public:
  typedef std::complex<T> complexType;   
  typedef T magnitudeType; 
  typedef typename numberTraits<T>::integralType integralType;  

};
#endif

#if defined(__USE_MERE)
// explicit specializations for mere::C, mere::R, mere::Z
template <> 
class numberTraits<mere::C>{
public:
  typedef mere::C complexType;   // complex type if T is real. 
  typedef mere::R magnitudeType; // magnitude type if T is complex.
  typedef mere::Z integralType;  // integer type of corresponding range.

};
template <> 
class numberTraits<mere::R>{
public:
  typedef mere::C complexType;   // complex type if T is real.
  typedef mere::R magnitudeType; // magnitude type if T is complex.
  typedef mere::Z integralType;  // integer type of corresponding range.

};
template <> 
class numberTraits<mere::Z>{
public:

  typedef mere::Z magnitudeType; // magnitude type if T is complex.
  typedef mere::Z integralType;  // integer type of corresponding range.

};
#endif

#if 0
template <> 
class numberTraits< std::complex<double> >{
public:
  typedef std::complex<double> complexType; // complex type if T is real.
  typedef double magnitudeType; // magnitude type if T is complex.
  typedef long integralType;  // integer type of corresponding range.
};
#endif
// *******************************

// explicit instantiation for long and unsigned long
template <> 
class numberTraits<long>{
public:

  typedef long magnitudeType; // magnitude type if T is complex.
  typedef long integralType;  // integer type of corresponding range.

};

template <> 
class numberTraits<unsigned long>{
public:

  typedef unsigned long magnitudeType; // magnitude type if T is complex.
  typedef unsigned long integralType;  // integer type of corresponding range.

};


template <class T, class N>
T integer(const N& n);

template <class T, class N>
T ratio(const N& numerator, const N& denominator);

template <class T>
const T& zero(void);

template <class T>
const T& epsilon(void);

#if defined(__USE_MERE)
template <>
mere::C integer<C,long>(const long& n);
template <>
const mere::C& zero<mere::C>(void);
template <>
const mere::C& epsilon<mere::C>(void);

template <>
mere::R integer<R,long>(const long& n);
template <>
const mere::R& zero<mere::R>(void);
template <>
const mere::R& epsilon<mere::R>(void);
#else
template <>
std::complex<double> integer<std::complex<double>,long>(const long& n);
template <>
const std::complex<double>& zero<std::complex<double> >(void);
template <>
const std::complex<double>& epsilon<std::complex<double> >(void);

template <>
double integer<double,long>(const long& n);
template <>
const double& zero<double>(void);
template <>
const double& epsilon<double>(void);
#endif

#else
  using number::numberTraits;
  using number::integer;
  using number::ratio;
  using number::zero;
  using number::one;
  using number::epsilon;
  using number::pi;
  using number::random;
#endif


#if !defined(__numericalFunctor__h) && 0
  /**
   * @brief Solve the quadratic formula:
   * @tparam V1: container type of value_type real or complex
   * @tparam V2: container type of value_type complex (or value_type real \f$ \Leftrightarrow \f$ roots are somehow constrained to be real)
   * @param[in]: vP: vector of polynomial coefficients
   * @param[out]: vRoot: vector of roots
   * @param[out]: D: discriminant
   */
  template <class V1, class V2>
  void quadraticFormula(const V1& vP, V2& vRoot, typename V1::value_type& D);  
#else
  using TMatrix::quadraticFormula;
#endif


using linalg::ntuple;
using linalg::ntuple_interval;
// ---------------------------------- end: section for TMatrix namespace independence -------------------------


template <class R, size_t NDIM>
class system{

  protected:
  
    // allow this module to be independent of namespace TMatrix:
    //   privately define static "readBinary_", "writeBinary_", and "binarySize_" methods for use within this class.

    #define _FULLY_QUALIFIED_CLASS_NAME_  system<R,NDIM>
    #define _OUTER_TEMPLATE_PREFIX_ template <class R, size_t NDIM>
    #include "binaryIO_stub.h"
    
    #if 1
    // definitions in "binaryIO_stub.h" do not compile, for some reason:
    static size_t binarySize_(const ntuple<size_t,NDIM>& t);
    static size_t binarySize_(const ntuple<long,NDIM>& t);
    static size_t binarySize_(const ntuple<R,NDIM>& t);
    static size_t binarySize_(const ntuple_interval<R,NDIM>& t);    
    #endif

    #if defined(__USE_PTHREAD)
    // general system<R,NDIM> mutex (protects access to system_stats_):
    mutable pthread_mutex_t mutex_;
    
    // mutexes and structures required for safe write-access to a pair of particles:
    mutable pthread_mutex_t pair_mutex_; 
    mutable pthread_cond_t pair_cond_;
    mutable parallelUtil::multiMutex
      particle0_mutex_, particle1_mutex_;
      
    // associated pair lock and unlock methods:
    
    /*
     * blocks until pair can be acquired simultaneously (keyed to addresses of p1 and p2);
     *   returns pair of keys to pass to matching "unlock_pair" method
     */
    std::pair<size_t, size_t> lock_pair_(void *p0, void *p1) throw(std::string);

    /*
     * uses pair of keys from matching "lock_pair" method to release aquired pair.
     */
    void unlock_pair_(const std::pair<size_t, size_t>& keys) throw(std::string);  
    #endif
      
  public:

    typedef typename numberTraits<R>::complexType C;
    typedef typename numberTraits<R>::integralType Z;
    typedef python_util::simple_object_base::object_list object_list;
    typedef python_util::simple_object_base::object_map object_map;
    typedef linalg::linear_interpolator<gmm::dense_vector_ref<const R*>,NDIM> interpolator;               

    class state;
    class particle; // forward
    class cell; // forward 
    class event; // forward    
    typedef std::list<particle*> particle_list_type;  
      
    typedef ntuple<long,NDIM> row_major_index_base_type;
    typedef ntuple<size_t,NDIM> shape_type;
    
        
    class parameters: public python_util::options_map<C>{
      public:
        
        typedef python_util::options_map<C> base_class;
        typedef typename base_class::object_list object_list;
        
        // -------------------- cached parameters: --------------------------------------
        
        // hyper-cube edge length:
        R L;

        size_t N_particle;

        #if 0
        // maximum evolution time:
        R t_max;
        #endif
       
        // maximum number of event steps:
        size_t NSTEP;
        
        R sticking_probability;

        #if 0
        size_t particle_per_cell;
                
        size_t N_seed;
        #endif
        
        // ------------------------------------------------------------------------------
        
        virtual void update_cache(bool derived=false, bool to_cache=true) throw(std::string);
        
        virtual void valid_check(void)const throw(std::string);
        
        virtual python_util::options_map<C>* clone(void)const throw(std::string);
        
        // --------- just wrap base_class::copy ------------
        void copy(const parameters& other);
        // -------------------------------------------------
        
        parameters& operator=(const parameters& other);
        parameters& operator=(const python_util::options_map<R>& other);
        
        virtual void write(std::ostream& os)const;
        
        #if 0 // ----------------- use base_class binary I/O -------------------
        virtual bool writeBinary(commUtil::abstractCommHandle *fp)const;
        virtual bool readBinary(commUtil::abstractCommHandle *fp);
        virtual size_t binarySize(void)const;
        #endif
        
        /** 
         * @brief Initialize from a simple_object_base*.
         */
        inline void extract(const python_util::simple_object_base* src) throw(std::string);
        
        virtual ~parameters(void);
        
        parameters(bool derived=false);
        parameters(const parameters& other);
        parameters(const python_util::options_map<C>& other);
    }; // class parameters
        
    const parameters& get_parameters(void)const;
    
    /**
     * @brief set values for fixed system parameters 
     * parameters:
     *   - @b L: hyper-cube edge length
     *   - @b N_particle: number of particles
     *   - @b t_max[opt]: maximum time for system evolution during "apply" or "step"
     *   - @b particle_per_cell[opt]: number of particles per kinematic cell
     *   - @b sticking_probability[opt]: per-collision stick probability for diffusion-limited-aggregation process
     *   - @b N_seed[opt]: number of seed (i.e. sticky) particles for diffusion-limited-aggregation process
     *   .
     */
    void set_parameters(const parameters& param) throw(std::string);
        
#if !defined(__ND_matrix__h)         
// ---------------------- moved to namespace linalg: (ND_matrix.h): -----------------------------

    // row-major index class (allowing += sub-indices):
    class row_major_index: public row_major_index_base_type{
        
        // product of indices in [offset,NDIM)
        // (where product_(NDIM) \defas 1)
        template <class N1>
        static N1 product_(const ntuple<N1,NDIM>& t, size_t offset=0);
        
      public:
        typedef row_major_index_base_type base_class;
        typedef typename base_class::value_type value_type;
        typedef typename system<R,NDIM>::shape_type shape_type;
        
        row_major_index operator+(const row_major_index& other)const;
        row_major_index operator-(const row_major_index& other)const;
        
        row_major_index& operator+=(const row_major_index& other);
        row_major_index& operator-=(const row_major_index& other);
        
        // apply periodic B.C. as constrained by shape:
        row_major_index& mod_assign(const shape_type& shape);
      
        row_major_index& operator=(const row_major_index& other);

        template <class N1>
        row_major_index& operator=(const ntuple<N1,NDIM>& other);
        
        size_t inverse(const shape_type& shape, value_type shift=0)const;
      
        /*
         * test wrap-around condition for index "x2" with respect to index "x1";
         *   initialize dimensions corresponding to wrapping with 0 or +-1, corresponding to wrap direction;
         *   return true if any dimensions wrap.
         */
        static bool wrap(const row_major_index& x1, const row_major_index& x2, const shape_type& shape,
                         row_major_index& wrap_dim);
       
        bool writeBinary(commUtil::abstractCommHandle* fp)const throw(std::string);
        bool readBinary(commUtil::abstractCommHandle* fp) throw(std::string);
        size_t binarySize(void)const throw(std::string);
        
        /*
         * return the row_major_index corresponding to the position "p" within the hyper-cube "domain_cube"
         *   given the array "shape" and index-base "shift"
         */
        static row_major_index ND_bin(const ntuple<R,NDIM>& p, const ntuple_interval<R,NDIM>& domain_cube,
                                      const shape_type& shape, value_type shift=0) throw(std::string);

        #if 0 // inclusion of "ND_edges" method completely screws-up the compiler for some reason (in unrelated code sections):
        /*
         * return the sub-hypercube corresponding to the row_major_index within the hyper-cube "domain_cube"
         *   given the array "shape" and index-base "shift"
         */
        static ntuple_interval<R,NDIM> ND_edges(const row_major_index& indices, const ntuple_interval<R,NDIM>& domain_cube,
                                                const shape_type& shape, value_type shift=0, const R& epsilon_=zero<R>()) throw(std::string);
        #endif
      
        bool in_domain(const shape_type& shape, value_type shift=0)const;


        row_major_index(value_type shift=0, bool clear=true);
      
        row_major_index(size_t n_linear, const shape_type& shape, value_type shift=0);
        
        row_major_index(const row_major_index& other);
        
        template <class N1>
        row_major_index(const ntuple<N1,NDIM>& other);  
    }; // class row_major_index
    
// ---------------------- end: moved to namespace linalg: (ND_matrix.h): -----------------------------
#else
    typedef linalg::row_major_index<long, NDIM> row_major_index;
#endif 
      
    class state{
      public:
        typedef typename system<R,NDIM>::particle_list_type particle_list_type;
        
      private:
      
        ntuple<R,NDIM> position_, velocity_;
        R time_;
        
        #if defined(state_partner_lists)
        particle_list_type partners_;
        #endif
        
        int event_type_; // typename system<R,NDIM>::event::EVENT_KIND as int
        
      public:
      
        const ntuple<R,NDIM>& position(void)const;
        ntuple<R,NDIM>& position(void);
      
        const ntuple<R,NDIM>& velocity(void)const;
        ntuple<R,NDIM>& velocity(void);
      
        const R& time(void)const;
        R& time(void);
        
        #if defined(state_partner_lists)
        const particle_list_type& partners(void)const;
        particle_list_type& partners(void);

        // add a particle to the partners-list, if not already in the list:
        void add_partner(typename system<R,NDIM>::particle *p);

        void write_partners(std::ostream& os)const;        
        #endif
        
        inline int event_type(void)const;
        inline void set_event_type(int type);
        
        void init(const ntuple<R,NDIM>& x, const ntuple<R,NDIM>& v, const R& t);
        
        state& operator=(const state& other);

        void write(std::ostream& os)const;
        
        bool writeBinary(commUtil::abstractCommHandle* fp)const throw(std::string);
        bool readBinary(commUtil::abstractCommHandle* fp) throw(std::string);
        size_t binarySize(void)const throw(std::string);
                
        state(void);
        state(const state& other);
        state(const ntuple<R,NDIM>& x, const ntuple<R,NDIM>& v, const R& t);         
    }; // class state
         
    class particle{
    
        #if 0
        static const size_t state_buf_size_ = 2;
        state state_[state_buf_size_];
        size_t current_state_;
        #else
        // *** DEBUG ***
        static const size_t state_buf_size_ = 10;
        std::deque<state> state_;
        #endif
        
        long linear_index_; // long so that "-" can be used to indicate a "wrapped" particle
        cell *owner_;
        
        bool frozen_;
        
      public:
      
        #if 0
        typedef const state* state_buf_iterator;
        #else
        // *** DEBUG ***
        typedef typename std::deque<state>::const_iterator state_buf_iterator;
        #endif
      
        // _ALL_ particle-derived classes must be in this list, 
        //   to allow proper function of readBinaryVirtual, and writeBinaryVirtual methods:
        enum particle_kind { ABSTRACT_PARTICLE, 
                             SPHERE_PARTICLE };

        virtual particle_kind kind(void)const;

        /*
         * stick to another particle (modifies state of _both_ particles):
         *   sets particle velocities to zero (and sets "frozen_" flags)
         *   (where applicable (i.e. _not_ in jam case), particles have been moved to kinematic coordinates of event,
         *      prior to this method)
         */
        void stick(event* E);
 

        const state& current_state(void)const;
        state& current_state(void);

        const state& previous_state(void)const;
        
        void rotate_state_buffer(void);

        // support for examination of entire state buffer:
        inline state_buf_iterator state_buf_begin(void)const;
        
        inline state_buf_iterator state_buf_end(void)const;
        
        void write_state_buf(std::ostream& os)const;
        
        
        long linear_index(void)const;
        long& linear_index(void);

        // both return ref-type for consistency:
        const cell* owner(void)const;
        cell*& owner(void);

        bool frozen(void)const;
        void freeze(bool flag=true);
        
        // id-string as "p<n>,c<n>":
        static std::string id_str(const particle* p, const cell* c);
        std::string id_str(void)const;
        
        // update any _intrinsic_ attributes associated with system time change.
        // (does _not_ modify particle "state" (here considered extrinsic))
        virtual void move(const R& t);
        
        // overlap condition test:
        virtual bool overlap(const particle* pother)const;
        
        /**
         * @brief Velocity magnitude corresponding to a specified kinetic temperature.
         *   This virtual method allows any method using temperature to initialize state velocity to be completely virtual.
         *   - for spheres, the zero-radius case uses density and "dr" (instead of radius) to calculate an effective velocity;
         *   - for abstract particles, mass is taken as unity.
         *   .
         */
        virtual R kinetic_velocity(const R& T)const throw(std::runtime_error);  
        
        particle& operator=(const particle& other);
        
        virtual particle* clone(void)const;

        // methods to allow binary read and write from pointer to base-class:            
        static bool writeBinaryVirtual(commUtil::abstractCommHandle* fp, const particle* p) throw(std::string);
        static bool readBinaryVirtual(commUtil::abstractCommHandle* fp, particle*& p) throw(std::string);
        static size_t binarySizeVirtual(const particle* p) throw(std::string);
        
        virtual bool writeBinary(commUtil::abstractCommHandle* fp)const throw(std::string);
        virtual bool readBinary(commUtil::abstractCommHandle* fp) throw(std::string);
        virtual size_t binarySize(void)const throw(std::string);

                
        virtual ~particle(void);
        
        particle(void);
        particle(const particle& other);        
        particle(const state& s, long linear_index, cell* owner);
    }; // class particle

    
    class cell{
      public:
        
        typedef typename system<R,NDIM>::particle_list_type particle_list_type;
        typedef std::vector<cell*> neighbor_list_type;
      private:
        
        particle_list_type particles_;
        
        neighbor_list_type neighbors_;
        
        neighbor_list_type positive_neighbors_;
        
        size_t linear_index_; 
        
        row_major_index indices_; 
        
        std::vector<int> wrap_; // to be used as std::vector<bool> (avoiding "bitvector" problems...)
        
        std::vector<row_major_index> wrap_dim_;

        // shape corresponding to neighbor offset row_major_index:
        static shape_type offset_shape_;
        
        // free owned data:
        void free_pointers_(void);
        
      public:
 
        // usage note: a cell _owns_ its particles, 
        //   transfer of particle pointer transfers ownership:
        const particle_list_type& particles(void)const;
        particle_list_type& particles(void);
        
        // a cell does _not_ own its neighbors:
        const neighbor_list_type& neighbors(void)const;
        neighbor_list_type& neighbors(void);
        
        const neighbor_list_type& positive_neighbors(void)const;
        neighbor_list_type& positive_neighbors(void);
        
        const std::vector<int>& wrap(void)const;
        std::vector<int>& wrap(void);
        
        const std::vector<row_major_index>& wrap_dim(void)const;
        std::vector<row_major_index>& wrap_dim(void);
        
        const size_t linear_index(void)const;
        const row_major_index& indices(void)const;

        
        // neighbor from relative index (i.e. sub-indices are in [-1,0,1] only):
        const cell* neighbor(const row_major_index& nx)const;
        cell* neighbor(const row_major_index& nx);

        // face neighbor from face index pair:
        const cell* face_neighbor(size_t dimension, size_t face)const;
        cell* face_neighbor(size_t dimension, size_t face);
          
        static size_t N_face(void);        
        static size_t N_neighbor(void);
        static size_t N_positive_neighbor(void);
        
        static const shape_type& offset_shape(void);
        
        // transfer a particle to another cell:
        void transfer_particle(particle* p, cell* pother) throw(std::string);
        
        // assignment or clone does _not_ transfer neighbor-lists (and associated boundary wrapping info).
        // (these can only be constructed by parent "cell_array")
        cell* clone(void)const;
        
        cell& operator=(const cell& other);

        
        bool writeBinary(commUtil::abstractCommHandle* fp)const throw(std::string);
        bool readBinary(commUtil::abstractCommHandle* fp) throw(std::string);
        size_t binarySize(void)const throw(std::string);

        
        ~cell(void);  
        
        cell(void);              

        // copy does _not_ initialize neighbor-lists (and associated boundary wrapping info).
        // (these can only be initialized  by parent "cell_array")
        cell(const cell& other) throw(std::string);
        
        cell(const size_t linear_index, const row_major_index& indices);                 
    }; // class cell

    
    class cell_array{
      public:
      
        typedef std::vector<cell*> cell_list_type;
        typedef typename system<R,NDIM>::row_major_index index;         
        typedef typename cell_list_type::iterator iterator;
        typedef typename cell_list_type::const_iterator const_iterator;      

      private:
      
        size_t N1_;
        
        shape_type shape_;

        std::vector<cell*> data_;
        
        // free owned data:
        void free_pointers_(void);
        
        // initialize cell neighbor-lists and boundary wrapping information:
        void init_neighbors_(void);
        
      public:
               
        // cells along hypercube edge:
        const size_t& N1(void)const;
        
        // total number of cells:
        const size_t& N_cell(void)const;
        
        const shape_type& shape(void)const;
        
        // usage note: a cell_array _owns_ its cells,
        //   transfer of cell pointer transfers ownership
        const cell_list_type& data(void)const;
        cell_list_type& data(void);

        // allow usage as container to be (somewhat) transparent:
        iterator begin(void);
        const_iterator begin(void)const;
        
        iterator end(void);
        const_iterator end(void)const;
        
        cell* operator[](size_t linear_index) throw(std::string);
        const cell* operator[](size_t linear_index)const throw(std::string);
        cell* operator[](const row_major_index& indices) throw(std::string);
        const cell* operator[](const row_major_index& indices)const throw(std::string);
        

        bool writeBinary(commUtil::abstractCommHandle* fp)const throw(std::string);
        bool readBinary(commUtil::abstractCommHandle* fp) throw(std::string);
        size_t binarySize(void)const throw(std::string);

        cell_array& operator=(const cell_array& other);

        ~cell_array(void);
        
        /**
         * @brief Initialize the cell array using a "suggested" value for the number-of-cells.
         *   the actual number of cells will be the value corresponding to the closest (>=) integer
         *   that is <some other integer>^NDIM (see additional comments at definition of "init").
         */
        void init(size_t N_cell);
        
        cell_array(void);
        cell_array(const cell_array& other);
        
        /*
         * initialize the cell array using a "suggested" value for the number-of-cells.
         *   the actual number of cells will be the value corresponding to the closest (>=) integer
         *   that is <some other integer>^NDIM
         */        
        cell_array(size_t N_cell);    
    }; // class cell_array


    // assume particle single or pair events:   
    class event{
      public:
      
        enum event_kind { NULL_EVENT=0, TIMESTEP_EVENT, MOVE_EVENT, CELL_EXIT_EVENT, COLLISION_EVENT, JAM_EVENT, STICK_EVENT };
        
      private:
      
        event_kind kind_;
      
        R time_;
        
        cell *cell_;
        
        particle *particle0_, *particle1_;
        
        size_t dimension_, face_;
        
        bool jam_;

        size_t id_;
        
        // static counter for class to generate unique ids:
        // (note: only "init" increments this counter ("operator=" and "clone" do _not_):
        static size_t event_counter_;
                
      public:

        static void write(std::ostream& os, event_kind kind_) throw(std::string);
        
        event_kind kind(void)const;


        size_t id(void)const;
        
        const R& time(void)const;

        const cell* event_cell(void)const;
        cell* event_cell(void);
        
        // usage note: event does _not_ own its cell or particle[s]
        //   for consistency with other classes these methods return pointers.
        const particle* particle0(void)const;
        particle* particle0(void);
        
        const particle* particle1(void)const;
        particle* particle1(void);
        
        const size_t& dimension(void)const;
        const size_t& face(void)const;

        bool is_collision(void)const;
        bool is_cell_exit(void)const;
        bool is_jam(void)const;
        bool is_stick(void)const;

        // duplicate events have same event_kind and same particles (or NULL particle*):
        //  (particle ordering is ignored)
        bool duplicate(const event* other)const;


        event* clone(void)const;
        
        event& operator=(const event& other);

#if 0 // ------------ event is fundamentally a pointer-based structure               ----------------
      // ------------   => implementing these methods doesn't really make any sense: ---------------- 

        bool writeBinary(commUtil::abstractCommHandle* fp)const throw(std::string);
        bool readBinary(commUtil::abstractCommHandle* fp) throw(std::string);
        size_t binarySize(void)const throw(std::string);
        
#endif // -------------------------------------------------------------------------------------------
        void write(std::ostream& os);
           
        void init(
              event_kind kind,
              const R& time, cell* c, particle *particle0, 
              particle *particle1=NULL,
              size_t dimension=0, size_t face=0);
                  
        ~event(void);
        
        event(void);
        event(const event& other);
        event(event_kind kind,
              const R& time, cell* c, particle *particle0, 
              particle *particle1=NULL,
              size_t dimension=0, size_t face=0);
                  
    }; // class event


    class event_list{
      public:
        typedef std::list<event*> list_type;
        typedef typename list_type::iterator iterator;
        typedef typename list_type::const_iterator const_iterator;            
      private:
        
        list_type data_;
        
        // free owned data:
        void free_pointers_(void);
        
      public:

        // allow (somewhat) transparent usage as container-type:
        const_iterator begin(void)const;
        iterator begin(void);
        
        const_iterator end(void)const;
        iterator end(void);  

        void push_back(event *E);
        
        // return values both defined as reference-type for consistency:
        const event*& back(void)const;
        event*& back(void);

        // push an event pointer onto front of list (transfers ownership of pointer-object): 
        void push_front(event *E);
        
        const event*& front(void)const;
        event*& front(void);

        // push a non-duplicate event object onto front of list (pointer-object is cloned if transfer occurs):
        // (duplicate events have same "kind()" and particle-pointers (regardless of particle order)) 
        void push_front_unique(const event *E);

      
        bool empty(void)const;
        
        // usage note: event_list _owns_ its events,
        //   transfer of event pointer transfers ownership:
        const list_type& data(void)const;
        list_type data(void);
        
        event_list& operator=(const event_list& other);

#if 0 // ------------ event is fundamentally a pointer-based structure               ----------------
      // ------------   => implementing these methods doesn't really make any sense: ---------------- 
        
        bool writeBinary(commUtil::abstractCommHandle* fp)const throw(std::string);
        bool readBinary(commUtil::abstractCommHandle* fp) throw(std::string);
        size_t binarySize(void)const throw(std::string);

#endif // -------------------------------------------------------------------------------------------
        void write(std::ostream& os)const;
        
        // remove list-position corresponding to next finite-time event (and free its object):
        void remove_next_event(void);
        
        // remove all events (and free associated objects):
        void clear(void);
        
        ~event_list(void);
        event_list(void);
        event_list(const event_list& other);    
    }; // class event_list
            
  private:
  
    parameters *pparam_;

    mutable std::string status_;

    // offset applied at "init_state" to obtain positive-definite positions from incoming state
    ntuple<R,NDIM> position_offset_; 
          
    // "t_max_" is the system-wide maximum time limit used for evolution completion determination (see "complete" method);
    //   this attribute must be distinct from the "parameters" attribute of the same name,
    //   as it is used slightly differently depending on whether the system is undergoing
    //     - free evolution: this t_max_ determines when "complete" returns true;
    //     - directed evolution to a target-distribution: this t_max_ also determines when
    //       "complete" returns true, but the sum of "dt" directed evolution intervals will
    //       in this case be equal to the parameters.t_max.
    //     . 
    mutable R t_max_; 
          
  protected:

    cell_array cells_;

  #if defined(event_history)
  // *** DEBUG ***
  public:  
    // event_history_ is primarily for debugging usage, and is not read/written by the I/O methods: 
    static event_list event_history_;
  
    static void write_event_history(std::ostream& os);
    
  protected:  
  #endif 
   
    const parameters* pparam(void)const;
    
    // protected initialization method:
    // (initialization relevent to local class only, no virtuals or super-classes)
    void init_(void);
      
    // next_event *always* is defined:
    virtual event_list next_event(const cell* c)const throw(std::string)=0;

    // either  all of "jam_particles", "move_particles", and "process_event"
    //   and/or "step_" *must* be defined (in the latter case, the former may remain stubs)
    virtual bool jam_particles(event* e) throw(std::string);

    virtual bool move_particles(cell* c, event* e) throw(std::string);
    
    virtual bool process_event(event* e) throw(std::string);
  
    virtual bool step_(event_list& events) throw(std::string);

    #if 0
    /*
     * create wrapped-clone of particle, using wrapping information
     *   from cell_array neighbor lists:
     */
    particle*  wrap_particle(const particle* src, bool wrap, const row_major_index& wrap_dim)const;
    #endif

    /* 
     * hooks for initialization and update of per-particle system statistics:
     * (these methods are "const" in that statistics is a measurement of system-state,
     *    (i.e. and therefore not part of the state))
     */
     virtual void init_statistics(void)const;
     virtual void update_statistics(const particle* p)const;

     /*
      * hook for application of system temperature control:
      * (note: this is a per-step hook, rather than a per-particle hook)
      */
     virtual R thermostat_v_factor(void)const;    

     /**
      * @brief Offset that was applied at "init_state" to obtain positive-definite positions from incoming state.
      */ 
     inline const ntuple<R,NDIM>& position_offset(void)const;
     
     /**
      * @brief Set current position offset (used by "init_state" to obtain positive definite positions):
      */ 
     inline void set_position_offset(const ntuple<R,NDIM>& offset); 

    /**
     * @brief Set the maximum-time limit value used by the completion test method (see "complete").
     *   This is a protected method: end-user applications set the parameters "t_max" attribute, and cannot use this method.
     */
    inline void set_t_max(const R& t)const;
    
    /**
     * @brief Time limit value used by the completion test method (see "complete").
     */
    inline const R& t_max(void)const ;
     
    /**
     * @brief protected apply method to be used by the other "apply" methods when a target-distribution is specified.
     *   - system parameters will be initialized prior to this method;
     *   - General system-state initialization will have been extracted from apply "arg" prior to this method;
     *   - at successful return from this method, system state will be ready for extraction to return value.
     *   .
     */
    virtual void apply_directed_(void)const throw(python_util::python_error, std::runtime_error);

    /**
     * @brief Fitness evaluator based on particle spatial distribution with respect to a target distribution, 
     *    specified as an N-dimensional interpolator.
     */
    R evaluate_fitness_(const interpolator& target_interp)const throw(std::runtime_error); 
    
    /**
     * @brief Re-initialize particle velocities based on the gradient of a target distribution.
     *   @param[in] target_gradient_interp  vector of interpolators for each component of the target-distribution gradient;
     *   @param[in] T  kinetic temperature of directed velocity component (randomized in the half-space);
     *   @param[in] T_B  Brownian temperature: kinetic temperature of non-directed component (completely randomized).
     */
    void direct_velocities_(const std::vector<interpolator>& target_gradient_interp, const R& T, const R& T_B) throw(std::runtime_error); 

  public:    
    
    const size_t DIM(void)const;
    
    const size_t N_particle(void)const;
    const R& L(void)const;


    // length corresponding to edge of cell-cube:
    R cell_edge_length(void)const;

    // center position of cell-cube:
    ntuple<R,NDIM> cell_center(const cell* cell_)const throw(std::string);

    // ntuple-interval representing cell-cube (intervals subject to left-closure \rightarrow  x \in [<start>, <end>)  ):
    ntuple_interval<R,NDIM> cell_cube(const cell* cell_)const throw(std::string);

    // ntuple-interval representing entire system domain (intervals subject to left-closure \rightarrow  x \in [<start>, <end>)  ):
    ntuple_interval<R,NDIM> system_cube(void)const throw(std::string);
    
    /*
     * test for particle cell-exit event:
     *   returns true if particle will exit cell
     *   dimension, face: face-pair corresponding to cube face of exit
     *   time: entry-time into the other cell (i.e. exit-time + epsilon)
     */
    bool cell_exit(const particle* particle_, const cell* cell_, R& time, size_t& dimension, size_t& face)const throw(std::string);    

    /*
     * move particle (deals with both extrinsic (i.e. state-related) and intrinsic changes):
     *   (calls particle.move() to update intrinsic change)
     */
    void move(particle* p, const R& t) throw(std::string);
    
    /*
     * implement a cell-exit event, including any required position wrap-around:
     * (note: presently transfer only occurs through one face at a time (corner transfers are zero cross-section events...))
     */
    void exit_cell(particle* p, size_t dimension, size_t face) throw(std::string);
             
    /**
     * @brief Evolve system through next event.
     * Return false if completion criteria are satisfied.
     */
    bool step(void) throw(std::string);
     
    virtual bool complete(void)const;

    // ================== apply related methods: ======================
    
    /**
     * @brief Check for consistency between a given argument and the present value of the functor parameters.
     * This method does not check for state consistency itself (e.g. lack of particle overlap); this is assumed.
     */
    virtual void arg_valid_check(const object_map& arg)const throw(std::string);
        
    /**
     * @brief generic apply method
     * @param[in] arg: 
     *    - @b position[opt]: list of particle positions (cartesian coordinates)
     *    - @b velocity[opt]: when corresponding positions are supplied, a list of particle vector velocities
     *    .
     * When particle positions are not supplied, values for number of particles "N_particle", velocity magnitude "v", radial growth rate "dr"
     *   will be taken from parameters to build the initial system state prior to evolution.
     * @param[out] val:
     *    - @b position: evolved particle positions
     *    - @b velocity: evolved particle vector velocities
     *    .    
     */
    virtual bool apply(const python_util::simple_object_base *arg, python_util::simple_object_base *val)const throw(std::string)=0;
    
    /// bool return version compatible with "classic" functor implementation:
    bool setParameters(const parameters& param);
    
    inline void setStatusString(const std::string& msg)const;
    
    inline const std::string& getStatusString(void)const;
           
    // =======================================================
    
    system& operator=(const system& other);

    bool writeBinary(commUtil::abstractCommHandle *fp)const throw(std::string);    
    bool readBinary(commUtil::abstractCommHandle *fp) throw(std::string);
    size_t binarySize(void)const throw(std::string);

    /**
     * @brief Initialize system state from parameters or optional incoming state.
     */
    virtual void init_state(const python_util::simple_object_base* incoming_state=NULL) throw(std::string)=0;
    
    /**
     * @brief Extract system state to external representation.
     */
    virtual void extract_state(python_util::simple_object_base *dest)const throw(std::string)=0;
    
    virtual ~system(void);
    
    system(bool initialize=true);
    system(const system& other);      
}; // class system


template <class R, size_t NDIM>
class mjp_system: public system<R,NDIM>{
   
  public:
    typedef system<R,NDIM> base_class;
    typedef typename base_class::shape_type shape_type;
    typedef typename base_class::row_major_index row_major_index;
    typedef typename base_class::state state;
    typedef typename base_class::particle particle;
    typedef typename base_class::cell cell;
    typedef typename base_class::cell_array cell_array;
    typedef typename base_class::event event;
    typedef typename base_class::event_list event_list;

    typedef typename numberTraits<R>::complexType C;
    typedef typename numberTraits<R>::integralType Z;
    typedef typename base_class::object_list object_list;
    typedef typename base_class::object_map object_map;
  
    class parameters: public system<R,NDIM>::parameters{
      public:
      
        typedef python_util::options_map<C> root_class;        
        typedef typename root_class::object_list object_list;        
        typedef typename system<R,NDIM>::parameters base_class;
        
        // -------------------- cached parameters: --------------------------------------

        R r_max; // if non-zero, used during: "complete"

        R T_max; // if non-zero, used to apply "thermostat" during "move"
        #if 0
        R v;
        
        R dr;                
        #endif
        // ------------------------------------------------------------------------------

        virtual void update_cache(bool derived=false, bool to_cache=true) throw(std::string);
           
        virtual void valid_check(void)const throw(std::string);
        
        virtual python_util::options_map<C>* clone(void)const throw(std::string);
        
        // --------- just wrap base_class::copy ------------
        void copy(const parameters& other);
        // -------------------------------------------------
        
        parameters& operator=(const parameters& other);
        parameters& operator=(const python_util::options_map<R>& other);
        
        virtual void write(std::ostream& os)const;
        
        #if 0 // ----------------- use base_class binary I/O -------------------
        virtual bool writeBinary(commUtil::abstractCommHandle *fp)const;
        virtual bool readBinary(commUtil::abstractCommHandle *fp);
        virtual size_t binarySize(void)const;
        #endif
        
        /** 
         * @brief Initialize from a simple_object_base*.
         */
        inline void extract(const python_util::simple_object_base* src) throw(std::string);  
              
        parameters(bool derived=false);
        parameters(const parameters& other);
        parameters(const python_util::options_map<C>& other);        
    }; // class parameters
        
    const parameters& get_parameters(void)const;
    
    /**
     * @brief set values for fixed mjp-system parameters 
     * parameters:
     *   - @b v[opt]: magnitude of initial particle velocity
     *        either a single-value must be present in parameters (to be used as a magnitude),
     *        or a list of vector particle velocities must be provided at "init_state" or "apply"
     *   - @b dr[opt]: magnitude of radial growth velocity
     *        either a single-value or a list of radial growth velocities may be provided
     *   - @b r_max[opt]: maximum particle growth radius
     *   - @b T_max[opt]: maximum particle temperature (kinetic temperature scale)
     *   - @b density[opt]: particle density
     *        either a single-value or a list of densities may be provided
     *   .
     */
    void set_parameters(const parameters& param) throw(std::string);


    // ------------- special intermediate-state methods (i.e. not a full re-initialization): ------------------ 

    /*
     * Reset all sphere intrinsic "dr" attributes to zero:
     */
    void freeze_growth(void);

    /*
     * Prepare system for diffusion-limited-aggregation (DLA) sequence: 
     *  -- reset radial growth-rate to zero
     *  -- fix seed particles
     */
    void init_DLA(void);
    
    // --------------------------------------------------------------------------------------------------------


    class sphere: public system<R,NDIM>::particle{
      
        R radius_, dr_, density_;
        
        static const R C_n_; // volume of unit-hyper-sphere
        
      public:
      
        typedef typename system<R,NDIM>::particle base_class;
#if 0
        typedef typename system<R,NDIM>::state state;
        typedef typename system<R,NDIM>::particle particle;
        typedef typename system<R,NDIM>::event event;
#endif
      
        const R& radius(void)const;
        R& radius(void);

        const R& dr(void)const;
        R& dr(void);

        const R& density(void)const;
        R& density(void);

        R volume(void)const;
        R hypervolume(void)const;

        static inline R volume(const R& radius);
        static inline R hypervolume(const R& radius);
        
        R mass(void)const;
        R hypermass(void)const;

        static inline R mass(const R& density, const R& radius);
        static inline R hypermass(const R& density, const R& radius);
        
        /*
         * kinetic temperature:
         */
        const R T(void)const;
      
        virtual typename base_class::particle_kind kind(void)const;

        // update any _intrinsic_ attributes associated with system time change.
        // (does _not_ modify particle "state" (here considered extrinsic))
        virtual void move(const R& t) throw(std::string);
        
        // overlap condition test:
        virtual bool overlap(const particle* pother)const;

        /**
         * @brief Velocity magnitude corresponding to a specified kinetic temperature.
         *   This virtual method allows any method using temperature to initialize state velocity to be completely virtual.
         *   - for spheres, the zero-radius case uses density and "dr" (instead of radius) to calculate an effective velocity;
         *   - for abstract particles, mass is taken as unity.
         *   .
         */
        virtual R kinetic_velocity(const R& T)const throw(std::runtime_error);  
 
        sphere& operator=(const sphere& other);
        
        virtual particle* clone(void)const;

        
        virtual bool writeBinary(commUtil::abstractCommHandle* fp)const throw(std::string);
        virtual bool readBinary(commUtil::abstractCommHandle* fp) throw(std::string);
        virtual size_t binarySize(void)const throw(std::string);

                
        virtual ~sphere(void);
        
        sphere(void);
        sphere(const sphere& other);
        
        // g++ 3.4.6 gives: "invalid type argument of `unary *'" if "one<R>()" is used as a default parameter value for "density": 
        sphere(const state& s, long linear_index, cell* owner, const R& r, const R& dr, const R& density /* =one<R>() */);       
    }; // class sphere

    class statistics{
    
        R r_min_, r_max_, dr_min_, dr_max_, v_min_, v_max_, t_min_, t_max_, T_min_, T_max_, T_mean_, hypervolume_;
        
        size_t N_frozen_, N_step_update_, N_total_update_;
        
      public:

        const R& r_min(void)const;
        const R& r_max(void)const;
        const R& dr_min(void)const;
        const R& dr_max(void)const;
        const R& v_min(void)const;
        const R& v_max(void)const;
        const R& t_min(void)const;
        const R& t_max(void)const;
        const R& T_min(void)const;
        const R& T_max(void)const;
        R T_mean(void)const;
        const R& hypervolume(void)const;
        
        const size_t& N_frozen(void)const;

        const size_t& N_step_update(void)const;
        const size_t& N_total_update(void)const;
              
        void update(const particle* p);
      
        void write(std::ostream& os)const;
      
        // re-initialize single-step only or all attributes:
        void clear(bool single_step = true);
        
        statistics(void);
    }; // struct statistics        
    
  private:
  
    mutable statistics system_stats_;
    
  protected:
        
    // next_event *always* is defined:
    virtual event_list next_event(const cell* c)const throw(std::string);

    // either  all of "jam_particles", "move_particles", and "process_event"
    //   and/or "step_" *must* be defined (in the latter case, the former may remain stubs)
    virtual bool jam_particles(event* e) throw(std::string);

    #if 0 // ------------- not special: moved up to base class: ------------------
    virtual bool move_particles(cell* c, event* e) throw(std::string);
    #endif // --------------------------------------------------------------------

    /**
     * @brief Calculate time for maximum radius limit condition.
     *   @param[in] p0  sphere to test
     *   @param[out] time time when maximum radius will be reached
     *   @return  true if time is finite
     */
    bool radius_limit(const sphere* p0, R& time)const;
    
    /**
     * @brief Test sphere-sphere pair collision.
     *   returns true if collision will occur
     *   time: absolute time of collision
     *   jam: spheres are touching with no possible relative motion
     */
    bool collision(const sphere* particle0_, const sphere* particle1_, R& time, bool& jam)const throw(std::string);

    /**
     * @brief Test sphere-sphere pair collision.
     *   This method only tests collision kinematics.
     *   returns true if collision will occur
     * @param[in]  p0: position of first particle
     * @param[in]  r0: radius of first particle
     * @param[in]  dr0: radial growth velocity of first particle
     * @param[in]  v0: velocity of first particle
     * @param[in]  p1: position of second particle
     * @param[in]  r1: radius of second particle
     * @param[in]  dr1: radial growth velocity of second particle
     * @param[in]  v1: velocity of second particle     
     * @param[out] dt: relative time of collision
     */
    static bool collision(
                  const ntuple<R,NDIM>& p0, const R& r0, const R& dr0,  const ntuple<R,NDIM>& v0,
                  const ntuple<R,NDIM>& p1, const R& r1, const R& dr1,  const ntuple<R,NDIM>& v1,
                  R& dt) throw(std::string);
  
    /*
     * implement sphere-sphere pair collision (modifies state of _both_ spheres, as appropriate):
     *   finite-time collisions and
     *   zero-time particle-jam events are treated uniformly
     *   (where applicable (i.e. _not_ in jam case), particles have been moved to kinematic coordinates of event,
     *      prior to this method)
     */
    void collide(event* E) throw(std::string);
    
    virtual bool process_event(event* e) throw(std::string);


      
   /* 
    * hooks for initialization and update of per-particle system statistics:
    * (these methods are "const" in that statistics is a measurement of system-state,
    *    (i.e. and therefore not part of the state))
    */
    virtual void init_statistics(void)const;
    virtual void update_statistics(const particle* p)const;

    /*
     * hook for application of system temperature control:
     * (note: this is a per-step hook, rather than a per-particle hook)
     */
    virtual R thermostat_v_factor(void)const;    

  public:    
    
    const statistics& system_stats(void)const;
    R packing_fraction(void)const;
     
    virtual bool complete(void)const;

    // ================== apply related methods: ======================
    
    /**
     * @brief Check for consistency between a given argument and the present value of the functor parameters.
     * This method does not check for state consistency itself (e.g. lack of particle overlap); this is assumed.
     */
    virtual void arg_valid_check(const object_map& arg)const throw(std::string);
    
    /// wrapper "apply" method:
    inline bool apply(const python_util::generic_object<C,R,Z>& arg, python_util::generic_object<C,R,Z>& val)const throw(std::string);
    
    /**
     * @brief generic apply method
     * @param[in] arg: 
     *    - @b position[opt]: list of particle positions (cartesian coordinates)
     *    - @b velocity[opt]: when corresponding positions are supplied, a list of particle vector velocities
     *    - @b radius[opt]: when corresponding positions are supplied, a list of particle radii
     *    .
     * When particle positions are not supplied, values for number of particles "N_particle", velocity magnitude "v", radial growth rate "dr"
     *   will be taken from parameters to build the initial system state prior to evolution.
     * @param[out] val:
     *    - @b position: evolved particle positions
     *    - @b velocity: evolved particle vector velocities
     *    - @b radius:   evolved particle radii
     *    .    
     */
    virtual bool apply(const python_util::simple_object_base *arg, python_util::simple_object_base *val)const throw(std::string);
         
    // =======================================================
        
    void write(std::ostream& os)const;
    void debug_print(void)const;

    bool writeBinary(commUtil::abstractCommHandle *fp)const throw(std::string);    
    bool readBinary(commUtil::abstractCommHandle *fp) throw(std::string);
    size_t binarySize(void)const throw(std::string);

    mjp_system& operator=(const mjp_system& other);


    /**
     * @brief Initialize system state from parameters or optional incoming state.
     *   @param[in] incoming_state:
     *      - @b position: list of particle positions
     *      - @b velocity[opt]: list of particle velocities
     *      - @b radius[opt]:   list of particle radii.
     *      .
     */    
    virtual void init_state(const python_util::simple_object_base *incoming_state=NULL) throw(std::string);

    /**
     * @brief Extract system state to external representation.
     *   @param[out] dest:
     *      - @b NDIM: system dimensions
     *      - @b parameters: system parameters (as object map)
     *      - @b position: list of particle positions
     *      - @b velocity: list of particle velocities
     *      - @b radius:   list of particle radii.
     *      - @b density:   list of particle densities.
     *      .
     */
    virtual void extract_state(python_util::simple_object_base *dest)const throw(std::string);

    
    virtual ~mjp_system(void);
        
    mjp_system(bool initialize=true);
    mjp_system(const mjp_system& other) throw(std::string);
    mjp_system(size_t N_particle, const R& L, const R& v, const R& dr) throw(std::string);  
}; // class mjp_system


} // namespace particle_packing

#include <particle_packing_inline.h>
#include <particle_packing_template.h>

#endif // __particle_packing__h
