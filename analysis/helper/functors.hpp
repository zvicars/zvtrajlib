#pragma once
#include "Eigen/Eigen"
#include "unsupported/Eigen/NonLinearOptimization"
using namespace Eigen;
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

  // you should define that in the subclass :
//  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};

struct logisticFunctor : Functor<double>
{
    Eigen::VectorXd x, y;
    logisticFunctor(Eigen::MatrixXd data): Functor<double>(3,Dynamic) {
      x = data.col(0);
      y = data.col(1);
      m_values = data.rows();
    };
    int operator()(const VectorXd &b, VectorXd &fvec)
    {
        assert(b.size()==3);
        for(int i=0; i < m_values; i++) {
            fvec[i] = b[0]/(1+exp(b[1]*(x[i] - b[2]))) - y[i];
        }
        return 0;
    }

    int df(const VectorXd &b, MatrixXd &fjac) const
    {
      assert(b.size() == 3);
      assert(fjac.rows() == values());
      for(int i = 0; i < values(); i++){
        Eigen::Vector3d jac_row;
        double dx = x[i] - b[2];
        jac_row << 1.0/(1+exp(b[1]*dx)),  -b[0]*(dx)*exp(b[1]*dx)/std::pow( 1 + exp(b[1]*dx), 2), b[0]*b[1]*exp(b[1]*dx)/std::pow( 1 + exp(b[1]*dx), 2);
        fjac.row(i) = jac_row;
      }
      return 0;
    }
};

struct logisticStepFunctor : Functor<double>
{
    Eigen::VectorXd x, y;
    logisticStepFunctor(Eigen::MatrixXd data): Functor<double>(4,Dynamic) {
      x = data.col(0);
      y = data.col(1);
      m_values = data.rows();
    };
    int operator()(const VectorXd &b, VectorXd &fvec)
    {
        assert(b.size()==4);
        for(int i=0; i<m_values; i++) {
            fvec[i] = b[0]/(1+exp(b[1]*(x[i] - b[2]))) - b[0]/(1+exp(b[1]*(x[i] - b[3]))) - y[i];
        }
        return 0;
    }

    int df(const VectorXd &b, MatrixXd &fjac) const
    {
      assert(b.size() == 4);
      assert(fjac.rows() == values());
      for(int i = 0; i < values(); i++){
        Eigen::Vector4d jac_row;
        double dx = x[i] - b[2];
        double dx2 = x[i] - b[3];
        jac_row << 1.0/(1+exp(b[1]*dx)) - 1.0/(1+exp(b[1]*dx2)),  
                   -b[0]*(dx)*exp(b[1]*dx)/std::pow( 1 + exp(b[1]*dx), 2) + b[0]*(dx2)*exp(b[1]*dx2)/std::pow( 1 + exp(b[1]*dx2), 2), 
                   b[0]*b[1]*exp(b[1]*dx)/std::pow( 1 + exp(b[1]*dx), 2),
                   -b[0]*b[1]*exp(b[1]*dx2)/std::pow( 1 + exp(b[1]*dx2), 2);
        fjac.row(i) = jac_row;
      }
      return 0;
    }
};