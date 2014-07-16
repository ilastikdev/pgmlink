#ifndef CONSTRAINT_FUNCTION_HXX
#define CONSTRAINT_FUNCTION_HXX

#include <stdexcept>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/functions/function_properties_base.hxx>


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
template<class T, class I=size_t, class L=size_t>
class ConstraintFunction:
        public FunctionBase<ExplicitFunction<T, I, L>, T, I, L>
{
public:
    typedef T ValueType;
    typedef L LabelType;
    typedef I IndexType;
    typedef typename opengm::FunctionBase<ConstraintFunction<T>, T, I, L> FunctionBaseType;

    /// operator is called to evaluate a certain labeling
    /// only compute the result when needed!
    template<class LABEL_ITERATOR>
    T operator()(LABEL_ITERATOR labels_begin, LABEL_ITERATOR labes_end)
    {
        std::vector<L> configuration;
        for(size_t i = 0; i < this->dimension(); ++i)
        {
            configuration.push_back(it[i]);
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


protected:
    virtual T get_energy_of_configuration(std::vector<L>& configuration)
    {
        throw std::logic_error("You have to use derived classes of ConstraintFunction!");
    }

    std::vector<typename I> shape_;
};

template<class T, class I=size_t, class L=size_t>
class IncomingConstraintFunction:ConstraintFunction<T,I,L>
{
public:
    // TODO: add constructor!
protected:
    virtual T get_energy_of_configuration(std::vector<L>& configuration)
    {
        throw std::logic_error("Not yet implemented!");
    }
};

template<class T, class I=size_t, class L=size_t>
class OutgoingConstraintFunction:ConstraintFunction<T,I,L>
{
public:
    // TODO: add constructor!
protected:
    virtual T get_energy_of_configuration(std::vector<L>& configuration)
    {
        throw std::logic_error("Not yet implemented!");
    }
};

template<class T, class I=size_t, class L=size_t>
class DetectionConstraintFunction:ConstraintFunction<T,I,L>
{
public:
    // TODO: add constructor!
protected:
    virtual T get_energy_of_configuration(std::vector<L>& configuration)
    {
        throw std::exception("Not yet implemented!");
    }
};

//------------------------------------------------------------------------
// ConstraintFactor
//------------------------------------------------------------------------
/// The constraint factor is the instantiation of a constraint function in the GM, linked
/// to some nodes.
template <class OGM_FUNCTION>
class ConstraintFactor : public OpengmFactor<OGM_FUNCTION>
{
    // define interface such that it can be read as factor or constraint?
    // or is this just the factor instanciation?
};


//------------------------------------------------------------------------
// ConstraintPool
//------------------------------------------------------------------------
/// Of the incoming and outgoing constraint type, we will have instatiations of different order.
/// The constraint pool contains a map indexed by the "constraint-function-signatures",
/// pointing to the function object. When adding functions to the model we look up
/// whether it is already present in here and just get the function identifier to create a factor.
template<class T, class I, class L> // valueType, Indextype, Labeltype
class ConstraintPool
{
public:
    void add_incoming_constraint(std::vector<I> transition_nodes, I disappearance_node);
    void add_outgoing_constraint(I appearance_node, I division_node, std::vector<I> transition_nodes);
    void add_detection_constraint(I appearance_node, I disappearance_node);

    /// This method either adds all the factors to the graphical model
    /// or - given INF is opengm::LpCplex and "not force_hardconstraint" - adds hard constraints
    /// to the CPLEX optimizer.
    template<class GM, class INF>
    void add_constraints_to_problem(GM& model, INF& optimizer);

    void set_big_m(T big_m);
private:
    // maps with functions with different signatures
    T big_m_;
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

#endif // CONSTRAINT_FUNCTION_HXX
