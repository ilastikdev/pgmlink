#ifndef CONSTRAINT_FUNCTION_HXX
#define CONSTRAINT_FUNCTION_HXX

#include <stdexcept>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/functions/function_properties_base.hxx>
#include <vector>

namespace pgmlink
{
namespace pgm
{

//------------------------------------------------------------------------
// ConstraintFunction
//------------------------------------------------------------------------
///
/// Constraint function encapsulates an original hard constraint such that it can be used as
/// soft constraint with big-M penalty in optimizers different from CPLEX.
///
/// Serves as parent class to IncomingConstraintFunction, OutgoingConstraintFunction,
/// and DetectionConstraintFunction
///
/// \ingroup functions
template<class T, class I, class L>
class ConstraintFunction: public opengm::FunctionBase<ConstraintFunction<T, I, L>, T, I, L>
{
public:
//    ConstraintFunction(){}

    typedef T ValueType;
    typedef L LabelType;
    typedef I IndexType;
    typedef typename opengm::FunctionBase<ConstraintFunction<T, I, L>, T, I, L> FunctionBaseType;

    template<class SHAPE_ITERATOR>
    ConstraintFunction(SHAPE_ITERATOR shape_begin,
                       SHAPE_ITERATOR shape_end):
        FunctionBaseType(),
        shape_(shape_begin, shape_end)
    {}

    ConstraintFunction() {}

    /// operator is called to evaluate a certain labeling
    /// only compute the result when needed!
    template<class LABEL_ITERATOR>
    T operator()(LABEL_ITERATOR labels) const
    {
        std::vector<L> configuration;
        for(size_t i = 0; i < this->dimension(); i++)
        {
            configuration.push_back(labels[i]);
        }

        return get_energy_of_configuration(configuration);
    }

    size_t shape(const size_t var_idx) const
    {
        return shape_[var_idx];
    }

    size_t size() const
    {
        size_t res = 1;
        for(std::vector<std::size_t>::const_iterator it = shape_.begin(); it != shape_.end(); ++it )
        {
            res *= (*it);
        }
        return res;
    }

    size_t dimension() const
    {
        return shape_.size();
    }


    void set_forbidden_energy(T forbidden_energy)
    {
        forbidden_energy_ = forbidden_energy;
    }

protected:
    virtual T get_energy_of_configuration(const std::vector<L>&) const
    {
        throw std::logic_error("You have to use derived classes of ConstraintFunction!");
    }

    std::vector<I> shape_;
    T forbidden_energy_;
};

//------------------------------------------------------------------------
// IncomingConstraintFunction
//------------------------------------------------------------------------
/// This class assumes that the labels it gets to compute the function value
/// are ordered such that the n transition nodes (T_1 .. T_n) come first,
/// followed by the disappearance node (V).
/// It must then hold that sum(T_1 .. T_n) = V
template<class T, class I, class L>
class IncomingConstraintFunction: public ConstraintFunction<T, I, L>
{
public:
    template<class SHAPE_ITERATOR>
    IncomingConstraintFunction(SHAPE_ITERATOR shape_begin,
                               SHAPE_ITERATOR shape_end):
        ConstraintFunction<T, I, L>(shape_begin, shape_end)
    {}

    IncomingConstraintFunction() {}
protected:
    virtual T get_energy_of_configuration(const std::vector<L>& configuration) const
    {
        assert(configuration.size() > 1);

        L num_disappearing_objects = configuration.back();

        // sum incoming transitions until one before the end
        L sum = 0;
        auto end_range = configuration.end();
        end_range--;
        for(auto it = configuration.begin(); it != end_range; ++it)
        {
            sum += *it;
        }

        if(sum == num_disappearing_objects)
        {
            return 0.0;
        }
        else
        {
            return this->forbidden_energy_;
        }
    }
};

//------------------------------------------------------------------------
// OutgoingConstraintFunction
//------------------------------------------------------------------------
/// The outgoing constraint expects first the label of the appearance node (A),
/// then the label of the division node (D), and finally n transition nodes
/// (T_1 .. T_n).
/// It must hold that sum(T_1 .. T_n) = A + D
template<class T, class I, class L>
class OutgoingConstraintFunction: public ConstraintFunction<T, I, L>
{
public:
    template<class SHAPE_ITERATOR>
    OutgoingConstraintFunction(SHAPE_ITERATOR shape_begin,
                               SHAPE_ITERATOR shape_end):
        ConstraintFunction<T, I, L>(shape_begin, shape_end),
        with_divisions_(true)
    {}

    OutgoingConstraintFunction() {}

    void set_with_divisions(bool enable)
    {
        with_divisions_ = enable;
    }

protected:
    virtual T get_energy_of_configuration(const std::vector<L>& configuration) const
    {
        assert(configuration.size() > 1);

        auto it = configuration.begin();
        L num_appearing_objects = *it++;
        L division = *it++;

        // sum outgoing transitions
        L sum = 0;
        for(; it != configuration.end(); ++it)
        {
            sum += *it;
        }

        if((sum == num_appearing_objects + division) &&
                (division != 1 || num_appearing_objects == 1) &&
                (division == 0 || with_divisions_))
        {
            return 0.0;
        }
        else
        {
            return this->forbidden_energy_;
        }
    }

protected:
    bool with_divisions_;
};

//------------------------------------------------------------------------
// OutgoingNoDivConstraintFunction
//------------------------------------------------------------------------
/// The outgoing constraint expects first the label of the appearance node (A),
/// NO division node (D), and finally n transition nodes
/// (T_1 .. T_n).
/// It must hold that sum(T_1 .. T_n) = A
template<class T, class I, class L>
class OutgoingNoDivConstraintFunction: public ConstraintFunction<T, I, L>
{
public:
    template<class SHAPE_ITERATOR>
    OutgoingNoDivConstraintFunction(SHAPE_ITERATOR shape_begin,
                                    SHAPE_ITERATOR shape_end):
        ConstraintFunction<T, I, L>(shape_begin, shape_end)
    {}

    OutgoingNoDivConstraintFunction() {}

protected:
    virtual T get_energy_of_configuration(const std::vector<L>& configuration) const
    {
        assert(configuration.size() > 1);

        auto it = configuration.begin();
        L num_appearing_objects = *it++;

        // sum outgoing transitions
        L sum = 0;
        for(; it != configuration.end(); ++it)
        {
            sum += *it;
        }

        if(sum == num_appearing_objects)
        {
            return 0.0;
        }
        else
        {
            return this->forbidden_energy_;
        }
    }
};

//------------------------------------------------------------------------
// DetectionConstraintFunction
//------------------------------------------------------------------------
/// expects a configuration of size 2, containing an appearance and a disappearance node
template<class T, class I, class L>
class DetectionConstraintFunction: public ConstraintFunction<T, I, L>
{
public:
    template<class SHAPE_ITERATOR>
    DetectionConstraintFunction(SHAPE_ITERATOR shape_begin,
                                SHAPE_ITERATOR shape_end):
        ConstraintFunction<T, I, L>(shape_begin, shape_end),
        with_appearance_(true),
        with_disappearance_(true),
        with_misdetections_(true)
    {}

    DetectionConstraintFunction() {}

    void set_with_misdetections(bool enable)
    {
        with_misdetections_ = enable;
    }

    void set_with_appearance(bool enable)
    {
        with_appearance_ = enable;
    }

    void set_with_disappearance(bool enable)
    {
        with_disappearance_ = enable;
    }

protected:
    virtual T get_energy_of_configuration(const std::vector<L>& configuration) const
    {
        assert(configuration.size() == 2);

        auto it = configuration.begin();
        L num_disappearing_objects = *it++;
        L num_appearing_objects = *it;

        if((num_appearing_objects == num_disappearing_objects ||
                (num_appearing_objects == 0 && with_disappearance_) ||
                (num_disappearing_objects == 0 && with_appearance_)) &&
                (num_appearing_objects > 0 || num_disappearing_objects > 0 || with_misdetections_))
        {
            return 0.0;
        }
        else
        {
            return this->forbidden_energy_;
        }
    }

protected:
    bool with_misdetections_;
    bool with_appearance_;
    bool with_disappearance_;
};

//------------------------------------------------------------------------
// FixNodeValueConstraintFunction
//------------------------------------------------------------------------
/// expects a configuration of size 1, containing a disappearance node
template<class T, class I, class L>
class FixNodeValueConstraintFunction: public ConstraintFunction<T, I, L>
{
public:
    template<class SHAPE_ITERATOR>
    FixNodeValueConstraintFunction(SHAPE_ITERATOR shape_begin,
                                SHAPE_ITERATOR shape_end):
        ConstraintFunction<T, I, L>(shape_begin, shape_end),
        value(0)
    {}

    FixNodeValueConstraintFunction() {}

    void set_desired_value(size_t val){ value = val; }

protected:
    virtual T get_energy_of_configuration(const std::vector<L>& configuration) const
    {
        assert(configuration.size() == 1);

        auto it = configuration.begin();
        L variable_value = *it;

        if(variable_value == value)
        {
            return 0.0;
        }
        else
        {
            return this->forbidden_energy_;
        }
    }

protected:
    size_t value;
};

} // namespace pgm
} // namespace pgmlink

//------------------------------------------------------------------------
// Serialization
//------------------------------------------------------------------------
// TODO: how to serialize constraints instead of factors?
//       OpenGM serializes the model through this function:
//       graphicalmodel_hdf5.hxx:330 -> template<class GM> save(GM, filepath datasetName).
//       Unfortunately each factor is serialized as (functionId, numVars, (var1,var2,...)),
//       so we cannot intercept the factor serialization process for constraints.
//
//       Serialize ConstraintPool and GraphicalModel? -> same hdf5 but different dataset names
//       Use vigra or opengm hdf5 interface?
namespace opengm
{

using pgmlink::pgm::IncomingConstraintFunction;
using pgmlink::pgm::OutgoingConstraintFunction;
using pgmlink::pgm::OutgoingNoDivConstraintFunction;
using pgmlink::pgm::DetectionConstraintFunction;
using pgmlink::pgm::FixNodeValueConstraintFunction;

//------------------------------------------------------------------------
/// \cond HIDDEN_SYMBOLS
/// FunctionRegistration
template<class T, class I, class L>
struct FunctionRegistration< IncomingConstraintFunction<T, I, L> >
{
    enum ID
    {
        Id = opengm::FUNCTION_TYPE_ID_OFFSET
    };
};

template<class T, class I, class L>
struct FunctionRegistration< OutgoingConstraintFunction<T, I, L> >
{
    enum ID
    {
        Id = opengm::FUNCTION_TYPE_ID_OFFSET
    };
};

template<class T, class I, class L>
struct FunctionRegistration< OutgoingNoDivConstraintFunction<T, I, L> >
{
    enum ID
    {
        Id = opengm::FUNCTION_TYPE_ID_OFFSET
    };
};

template<class T, class I, class L>
struct FunctionRegistration< DetectionConstraintFunction<T, I, L> >
{
    enum ID
    {
        Id = opengm::FUNCTION_TYPE_ID_OFFSET
    };
};

template<class T, class I, class L>
struct FunctionRegistration< FixNodeValueConstraintFunction<T, I, L> >
{
    enum ID
    {
        Id = opengm::FUNCTION_TYPE_ID_OFFSET
    };
};

template<>
struct FunctionRegistration< marray::Marray<double> >
{
    enum ID
    {
        Id = opengm::FUNCTION_TYPE_ID_OFFSET
    };
};

//------------------------------------------------------------------------
/// Serialization for the incoming constraint function
template<class T, class I, class L>
class FunctionSerialization< IncomingConstraintFunction<T, I, L> >
{
public:
    static size_t indexSequenceSize(const IncomingConstraintFunction<T, I, L> &);
    static size_t valueSequenceSize(const IncomingConstraintFunction<T, I, L> &);

    template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
    static void serialize(const IncomingConstraintFunction<T, I, L>  &, INDEX_OUTPUT_ITERATOR, VALUE_OUTPUT_ITERATOR );

    template<class INDEX_INPUT_ITERATOR , class VALUE_INPUT_ITERATOR>
    static void deserialize( INDEX_INPUT_ITERATOR, VALUE_INPUT_ITERATOR, IncomingConstraintFunction<T, I, L>  &);
};

template<class T, class I, class L>
inline size_t FunctionSerialization<IncomingConstraintFunction<T, I, L> >::indexSequenceSize
(
    const IncomingConstraintFunction<T, I, L> & src
)
{
    return src.dimension() + 1;
}

template<class T, class I, class L>
inline size_t FunctionSerialization<IncomingConstraintFunction<T, I, L> >::valueSequenceSize
(
    const IncomingConstraintFunction<T, I, L> & src
)
{
    return src.size();
}

template<class T, class I, class L>
template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
void FunctionSerialization< IncomingConstraintFunction<T, I, L> >::serialize
(
    const IncomingConstraintFunction<T, I, L> & src,
    INDEX_OUTPUT_ITERATOR indexOutIterator,
    VALUE_OUTPUT_ITERATOR valueOutIterator
)
{
    //TODO implement me
    // see opengm::ExplicitFunction -> FunctionSerialization
    throw std::logic_error("not yet implemented");
}

template<class T, class I, class L>
template<class INDEX_INPUT_ITERATOR, class VALUE_INPUT_ITERATOR >
void FunctionSerialization<IncomingConstraintFunction<T, I, L> >::deserialize
(
    INDEX_INPUT_ITERATOR indexOutIterator,
    VALUE_INPUT_ITERATOR valueOutIterator,
    IncomingConstraintFunction<T, I, L> & dst
)
{
    //TODO implement me
    throw std::logic_error("not yet implemented");
}

//------------------------------------------------------------------------
/// Serialization for the outgoing constraint function
template<class T, class I, class L>
class FunctionSerialization< OutgoingConstraintFunction<T, I, L> >
{
public:
    static size_t indexSequenceSize(const OutgoingConstraintFunction<T, I, L> &);
    static size_t valueSequenceSize(const OutgoingConstraintFunction<T, I, L> &);

    template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
    static void serialize(const OutgoingConstraintFunction<T, I, L>  &, INDEX_OUTPUT_ITERATOR, VALUE_OUTPUT_ITERATOR );

    template<class INDEX_INPUT_ITERATOR , class VALUE_INPUT_ITERATOR>
    static void deserialize( INDEX_INPUT_ITERATOR, VALUE_INPUT_ITERATOR, OutgoingConstraintFunction<T, I, L>  &);
};

template<class T, class I, class L>
inline size_t FunctionSerialization<OutgoingConstraintFunction<T, I, L> >::indexSequenceSize
(
    const OutgoingConstraintFunction<T, I, L> & src
)
{
    return src.dimension() + 1;
}

template<class T, class I, class L>
inline size_t FunctionSerialization<OutgoingConstraintFunction<T, I, L> >::valueSequenceSize
(
    const OutgoingConstraintFunction<T, I, L> & src
)
{
    return src.size();
}

template<class T, class I, class L>
template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
void FunctionSerialization< OutgoingConstraintFunction<T, I, L> >::serialize
(
    const OutgoingConstraintFunction<T, I, L> & src,
    INDEX_OUTPUT_ITERATOR indexOutIterator,
    VALUE_OUTPUT_ITERATOR valueOutIterator
)
{
    //TODO implement me
    // see opengm::ExplicitFunction -> FunctionSerialization
    throw std::logic_error("not yet implemented");
}

template<class T, class I, class L>
template<class INDEX_INPUT_ITERATOR, class VALUE_INPUT_ITERATOR >
void FunctionSerialization<OutgoingConstraintFunction<T, I, L> >::deserialize
(
    INDEX_INPUT_ITERATOR indexOutIterator,
    VALUE_INPUT_ITERATOR valueOutIterator,
    OutgoingConstraintFunction<T, I, L> & dst
)
{
    //TODO implement me
    throw std::logic_error("not yet implemented");
}

//------------------------------------------------------------------------
/// Serialization for the outgoing no division constraint function
template<class T, class I, class L>
class FunctionSerialization< OutgoingNoDivConstraintFunction<T, I, L> >
{
public:
    static size_t indexSequenceSize(const OutgoingNoDivConstraintFunction<T, I, L> &);
    static size_t valueSequenceSize(const OutgoingNoDivConstraintFunction<T, I, L> &);

    template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
    static void serialize(const OutgoingNoDivConstraintFunction<T, I, L>  &, INDEX_OUTPUT_ITERATOR, VALUE_OUTPUT_ITERATOR );

    template<class INDEX_INPUT_ITERATOR , class VALUE_INPUT_ITERATOR>
    static void deserialize( INDEX_INPUT_ITERATOR, VALUE_INPUT_ITERATOR, OutgoingNoDivConstraintFunction<T, I, L>  &);
};

template<class T, class I, class L>
inline size_t FunctionSerialization<OutgoingNoDivConstraintFunction<T, I, L> >::indexSequenceSize
(
    const OutgoingNoDivConstraintFunction<T, I, L> & src
)
{
    return src.dimension() + 1;
}

template<class T, class I, class L>
inline size_t FunctionSerialization<OutgoingNoDivConstraintFunction<T, I, L> >::valueSequenceSize
(
    const OutgoingNoDivConstraintFunction<T, I, L> & src
)
{
    return src.size();
}

template<class T, class I, class L>
template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
void FunctionSerialization< OutgoingNoDivConstraintFunction<T, I, L> >::serialize
(
    const OutgoingNoDivConstraintFunction<T, I, L> & src,
    INDEX_OUTPUT_ITERATOR indexOutIterator,
    VALUE_OUTPUT_ITERATOR valueOutIterator
)
{
    //TODO implement me
    // see opengm::ExplicitFunction -> FunctionSerialization
    throw std::logic_error("not yet implemented");
}

template<class T, class I, class L>
template<class INDEX_INPUT_ITERATOR, class VALUE_INPUT_ITERATOR >
void FunctionSerialization<OutgoingNoDivConstraintFunction<T, I, L> >::deserialize
(
    INDEX_INPUT_ITERATOR indexOutIterator,
    VALUE_INPUT_ITERATOR valueOutIterator,
    OutgoingNoDivConstraintFunction<T, I, L> & dst
)
{
    //TODO implement me
    throw std::logic_error("not yet implemented");
}

//------------------------------------------------------------------------
/// Serialization for the detection constraint function
template<class T, class I, class L>
class FunctionSerialization< DetectionConstraintFunction<T, I, L> >
{
public:
    static size_t indexSequenceSize(const DetectionConstraintFunction<T, I, L> &);
    static size_t valueSequenceSize(const DetectionConstraintFunction<T, I, L> &);

    template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
    static void serialize(const DetectionConstraintFunction<T, I, L>  &, INDEX_OUTPUT_ITERATOR, VALUE_OUTPUT_ITERATOR );

    template<class INDEX_INPUT_ITERATOR , class VALUE_INPUT_ITERATOR>
    static void deserialize( INDEX_INPUT_ITERATOR, VALUE_INPUT_ITERATOR, DetectionConstraintFunction<T, I, L>  &);
};

template<class T, class I, class L>
inline size_t FunctionSerialization<DetectionConstraintFunction<T, I, L> >::indexSequenceSize
(
    const DetectionConstraintFunction<T, I, L> & src
)
{
    return src.dimension() + 1;
}

template<class T, class I, class L>
inline size_t FunctionSerialization<DetectionConstraintFunction<T, I, L> >::valueSequenceSize
(
    const DetectionConstraintFunction<T, I, L> & src
)
{
    return src.size();
}

template<class T, class I, class L>
template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
void FunctionSerialization< DetectionConstraintFunction<T, I, L> >::serialize
(
    const DetectionConstraintFunction<T, I, L> & src,
    INDEX_OUTPUT_ITERATOR indexOutIterator,
    VALUE_OUTPUT_ITERATOR valueOutIterator
)
{
    //TODO implement me
    // see opengm::ExplicitFunction -> FunctionSerialization
    throw std::logic_error("not yet implemented");
}

template<class T, class I, class L>
template<class INDEX_INPUT_ITERATOR, class VALUE_INPUT_ITERATOR >
void FunctionSerialization<DetectionConstraintFunction<T, I, L> >::deserialize
(
    INDEX_INPUT_ITERATOR indexOutIterator,
    VALUE_INPUT_ITERATOR valueOutIterator,
    DetectionConstraintFunction<T, I, L> & dst
)
{
    //TODO implement me
    throw std::logic_error("not yet implemented");
}

//------------------------------------------------------------------------
/// Serialization for the incoming constraint function
template<>
class FunctionSerialization< marray::Marray<double> >
{
public:
    static size_t indexSequenceSize(const marray::Marray<double> &);
    static size_t valueSequenceSize(const marray::Marray<double> &);

    template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
    static void serialize(const marray::Marray<double>  &, INDEX_OUTPUT_ITERATOR, VALUE_OUTPUT_ITERATOR );

    template<class INDEX_INPUT_ITERATOR , class VALUE_INPUT_ITERATOR>
    static void deserialize( INDEX_INPUT_ITERATOR, VALUE_INPUT_ITERATOR, marray::Marray<double>  &);
};

inline size_t FunctionSerialization<marray::Marray<double> >::indexSequenceSize
(
    const marray::Marray<double> & src
)
{
    return src.dimension() + 1;
}

inline size_t FunctionSerialization<marray::Marray<double> >::valueSequenceSize
(
    const marray::Marray<double> & src
)
{
    return src.size();
}

template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
void FunctionSerialization< marray::Marray<double> >::serialize
(
    const marray::Marray<double> & src,
    INDEX_OUTPUT_ITERATOR indexOutIterator,
    VALUE_OUTPUT_ITERATOR valueOutIterator
)
{
    //TODO implement me
    // see opengm::ExplicitFunction -> FunctionSerialization
    throw std::logic_error("not yet implemented");
}

template<class INDEX_INPUT_ITERATOR, class VALUE_INPUT_ITERATOR >
void FunctionSerialization<marray::Marray<double> >::deserialize
(
    INDEX_INPUT_ITERATOR indexOutIterator,
    VALUE_INPUT_ITERATOR valueOutIterator,
    marray::Marray<double> & dst
)
{
    //TODO implement me
    throw std::logic_error("not yet implemented");
}

//------------------------------------------------------------------------
/// Serialization for the outgoing constraint function
template<class T, class I, class L>
class FunctionSerialization< FixNodeValueConstraintFunction<T, I, L> >
{
public:
    static size_t indexSequenceSize(const FixNodeValueConstraintFunction<T, I, L> &);
    static size_t valueSequenceSize(const FixNodeValueConstraintFunction<T, I, L> &);

    template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
    static void serialize(const FixNodeValueConstraintFunction<T, I, L>  &, INDEX_OUTPUT_ITERATOR, VALUE_OUTPUT_ITERATOR );

    template<class INDEX_INPUT_ITERATOR , class VALUE_INPUT_ITERATOR>
    static void deserialize( INDEX_INPUT_ITERATOR, VALUE_INPUT_ITERATOR, FixNodeValueConstraintFunction<T, I, L>  &);
};

template<class T, class I, class L>
inline size_t FunctionSerialization<FixNodeValueConstraintFunction<T, I, L> >::indexSequenceSize
(
    const FixNodeValueConstraintFunction<T, I, L> & src
)
{
    return src.dimension() + 1;
}

template<class T, class I, class L>
inline size_t FunctionSerialization<FixNodeValueConstraintFunction<T, I, L> >::valueSequenceSize
(
    const FixNodeValueConstraintFunction<T, I, L> & src
)
{
    return src.size();
}

template<class T, class I, class L>
template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
void FunctionSerialization< FixNodeValueConstraintFunction<T, I, L> >::serialize
(
    const FixNodeValueConstraintFunction<T, I, L> & src,
    INDEX_OUTPUT_ITERATOR indexOutIterator,
    VALUE_OUTPUT_ITERATOR valueOutIterator
)
{
    //TODO implement me
    // see opengm::ExplicitFunction -> FunctionSerialization
    throw std::logic_error("not yet implemented");
}

template<class T, class I, class L>
template<class INDEX_INPUT_ITERATOR, class VALUE_INPUT_ITERATOR >
void FunctionSerialization<FixNodeValueConstraintFunction<T, I, L> >::deserialize
(
    INDEX_INPUT_ITERATOR indexOutIterator,
    VALUE_INPUT_ITERATOR valueOutIterator,
    FixNodeValueConstraintFunction<T, I, L> & dst
)
{
    //TODO implement me
    throw std::logic_error("not yet implemented");
}

} // namespace opengm

#endif // CONSTRAINT_FUNCTION_HXX
