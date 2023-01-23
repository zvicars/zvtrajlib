#pragma once
#include "Eigen/Eigen"
#include "unsupported/Eigen/NonLinearOptimization"

template<typename _Scalar, int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

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
    logisticFunctor(Eigen::MatrixXd data): Functor<double>(3,Eigen::Dynamic) {
      x = data.col(0);
      y = data.col(1);
      m_values = data.rows();
    };
    int operator()(const Eigen::VectorXd &b, Eigen::VectorXd &fvec)
    {
        assert(b.size()==3);
        for(int i=0; i < m_values; i++) {
            fvec[i] = b[0]/(1+exp(b[1]*(x[i] - b[2]))) - y[i];
        }
        return 0;
    }

    int df(const Eigen::VectorXd &b, Eigen::MatrixXd &fjac) const
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
    std::array<bool,5> fix;
    logisticStepFunctor(Eigen::MatrixXd data, std::array<bool,5> fix_in): Functor<double>(5,Eigen::Dynamic) {
      x = data.col(0);
      y = data.col(1);
      fix = fix_in;
      m_values = data.rows();
    };
    int operator()(const Eigen::VectorXd &b, Eigen::VectorXd &fvec)
    {
        assert(b.size()==5);
        for(int i=0; i<m_values; i++) {
            fvec[i] = b[0]/(1+exp(b[1]*(x[i] - b[2]))) - b[0]/(1+exp(b[1]*(x[i] - b[3]))) + b[4] - y[i];
        }
        return 0;
    }

    int df(const Eigen::VectorXd &b, Eigen::MatrixXd &fjac) const
    {
      assert(b.size() == 5);
      assert(fjac.rows() == values());
      for(int i = 0; i < values(); i++){
        Eigen::VectorXd jac_row(5);
        double dx = x[i] - b[2];
        double dx2 = x[i] - b[3];
        jac_row << 1.0/(1+exp(b[1]*dx)) - 1.0/(1+exp(b[1]*dx2)),  
                   -b[0]*(dx)*exp(b[1]*dx)/std::pow( 1 + exp(b[1]*dx), 2) + b[0]*(dx2)*exp(b[1]*dx2)/std::pow( 1 + exp(b[1]*dx2), 2), 
                   b[0]*b[1]*exp(b[1]*dx)/std::pow( 1 + exp(b[1]*dx), 2),
                   -b[0]*b[1]*exp(b[1]*dx2)/std::pow( 1 + exp(b[1]*dx2), 2), 
                   1;
        for(int j = 0; j < 5; j++){
          if(fix[j] == 1) jac_row[j] = 0;
        }
        fjac.row(i) = jac_row;
      }
      return 0;
    }
};

 /* 
  Eigen::Vector3d radangles;
  radangles << angles[2]*M_PI/180.0, angles[1] * M_PI / 180.0, angles[0] * M_PI / 180.0;
  Eigen::AngleAxisd rollAngle(radangles[0], Eigen::Vector3d::UnitZ());
  Eigen::AngleAxisd yawAngle(radangles[1], Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd pitchAngle(radangles[2], Eigen::Vector3d::UnitX());
  Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;
  Eigen::Matrix3d rotationMatrix = q.matrix();
  rotationMatrix.transposeInPlace();
  */
