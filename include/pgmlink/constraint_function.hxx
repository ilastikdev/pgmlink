#ifndef CONSTRAINT_FUNCTION_HXX
#define CONSTRAINT_FUNCTION_HXX

#include <stdexcept>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/functions/function_properties_base.hxx>
#include <vector>

namespace pgmlink {
namespace pgm {

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
    typedef T ValueType;
    typedef L LabelType;
    typedef I IndexType;
    typedef typename opengm::FunctionBase<ConstraintFunction<T,I,L>, T, I, L> FunctionBaseType;

    template<class SHAPE_ITERATOR>
    ConstraintFunction(SHAPE_ITERATOR shape_begin,
                       SHAPE_ITERATOR shape_end):
        FunctionBaseType(),
        shape_(shape_begin, shape_end)
    {}

    /// operator is called to evaluate a certain labeling
    /// only compute the result when needed!
    template<class LABEL_ITERATOR>
    T operator()(LABEL_ITERATOR labels_begin, LABEL_ITERATOR labels_end)
    {
        std::vector<L> configuration(labels_begin, labels_end);
        assert(configuration.size() == this->dimension());

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
    virtual T get_energy_of_configuration(const std::vector<L>&)
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
class IncomingConstraintFunction: public ConstraintFunction<T,I,L>
{
public:
    template<class SHAPE_ITERATOR>
    IncomingConstraintFunction(SHAPE_ITERATOR shape_begin,
                               SHAPE_ITERATOR shape_end):
        ConstraintFunction<T,I,L>(shape_begin, shape_end)
    {}

protected:
    virtual T get_energy_of_configuration(const std::vector<L>& configuration)
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
            return 0.0;
        else
            return this->forbidden_energy_;
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
class OutgoingConstraintFunction: public ConstraintFunction<T,I,L>
{
public:
    template<class SHAPE_ITERATOR>
    OutgoingConstraintFunction(SHAPE_ITERATOR shape_begin,
                               SHAPE_ITERATOR shape_end):
        ConstraintFunction<T,I,L>(shape_begin, shape_end)
    {}
protected:
    virtual T get_energy_of_configuration(const std::vector<L>& configuration)
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

        if(sum == num_appearing_objects + division && (division != 1 || num_appearing_objects == 1))
            return 0.0;
        else
            return this->forbidden_energy_;
    }
};

//------------------------------------------------------------------------
// DetectionConstraintFunction
//------------------------------------------------------------------------
/// expects a configuration of size 2, containing an appearance and a disappearance node
template<class T, class I, class L>
class DetectionConstraintFunction: public ConstraintFunction<T,I,L>
{
public:
    template<class SHAPE_ITERATOR>
    DetectionConstraintFunction(SHAPE_ITERATOR shape_begin,
                               SHAPE_ITERATOR shape_end):
        ConstraintFunction<T,I,L>(shape_begin, shape_end)
    {}
protected:
    virtual T get_energy_of_configuration(const std::vector<L>& configuration)
    {
        assert(configuration.size() == 2);

        auto it = configuration.begin();
        L num_disappearing_objects = *it++;
        L num_appearing_objects = *it;

        if(num_appearing_objects == num_disappearing_objects || num_appearing_objects == 0 || num_disappearing_objects == 0)
            return 0.0;
        else
            return this->forbidden_energy_;
    }
};

} // namespace pgm
} // namespace pgmlink

/*
 Probably not needed
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

using pgmlink::ConstraintFunction;

/// Serialization for the constraint function
template<class T, class I, class L>
class FunctionSerialization< ConstraintFunction<T, I, L> >
{
    static size_t indexSequenceSize(const ConstraintFunction<T, I, L> &);
    static size_t valueSequenceSize(const ConstraintFunction<T, I, L> &);

    template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
    static void serialize(const ConstraintFunction<T, I, L>  &, INDEX_OUTPUT_ITERATOR, VALUE_OUTPUT_ITERATOR );

    template<class INDEX_INPUT_ITERATOR , class VALUE_INPUT_ITERATOR>
    static void deserialize( INDEX_INPUT_ITERATOR, VALUE_INPUT_ITERATOR, ConstraintFunction<T, I, L>  &);
};

template<class T, class I, class L>
inline size_t FunctionSerialization<ConstraintFunction<T, I, L> >::indexSequenceSize
(
   const ConstraintFunction<T, I, L> & src
)
{
   return src.dimension() +1;
}

template<class T, class I, class L>
inline size_t FunctionSerialization<ConstraintFunction<T, I, L> >::valueSequenceSize
(
   const ConstraintFunction<T, I, L> & src
)
{
   return src.size();
}

template<class T, class I, class L>
template<class INDEX_OUTPUT_ITERATOR, class VALUE_OUTPUT_ITERATOR >
void FunctionSerialization< ConstraintFunction<T, I, L> >::serialize
(
   const ConstraintFunction<T, I, L> & src,
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
void FunctionSerialization<ConstraintFunction<T, I, L> >::deserialize
(
   INDEX_INPUT_ITERATOR indexOutIterator,
   VALUE_INPUT_ITERATOR valueOutIterator,
   ConstraintFunction<T, I, L> & dst
)
{
    //TODO implement me
    throw std::logic_error("not yet implemented");
}

} // namespace opengm
*/
#endif // CONSTRAINT_FUNCTION_HXX
