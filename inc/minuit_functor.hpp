//Wrapper template class that translates between any C++ class with operator()(std::vector<double>) (optionally const and/or passed by reference) and the function type required by TMinuit

#ifndef H_MINUIT_FUNCTOR
#define H_MINUIT_FUNCTOR

#include<cstddef>
#include<vector>

template<typename Functor>
class MinuitFunctor{
public:
  typedef void (*minuit_fcn_type)(int&, double*, double&, double*, int);
  typedef Functor functor_type;

  static void Function(int &npar, double *deriv, double &f, double *par, int flag){
    //This function is passed to a TMinuit via SetFCN
    if(false && npar && deriv && flag){
      //Stop compiler from complaining about unused variables
    }
    if(functor_ptr_!=NULL){
      std::vector<double> vec(par, par+num_params_);
      f=(*functor_ptr_)(vec);
    }
  }

  static void SetFunctor(Functor* const functor_ptr_in){
    functor_ptr_=functor_ptr_in;
  }

  static Functor* GetFunctor(){
    return functor_ptr_;
  }

  static void SetNumParams(const unsigned int& num_params_in){
    num_params_=num_params_in;
  }

  static unsigned int GetNumParams(){
    return num_params_;
  }

private:
  static Functor* functor_ptr_;
  static unsigned int num_params_;

  //Prevent instantiation (Wouldn't cause harm, but encourages confusing design)
  MinuitFunctor(){};
  MinuitFunctor(MinuitFunctor const&){};
  void operator=(MinuitFunctor const&){};
};

template<typename Functor>
Functor* MinuitFunctor<Functor>::functor_ptr_(NULL);

template<typename Functor>
unsigned int MinuitFunctor<Functor>::num_params_(0u);

#endif
