#include <Rcpp.h>
using namespace Rcpp;

#include <utility>
#include <vector>
#include <iostream>

class Object {
public:
    class Concept
    {
    public:
        virtual ~Concept() = default;
        virtual void Print() = 0;
    protected:
        Concept() = default;
    };  
    template <typename T>
    class Model : public Concept {
    public:
        template <typename U>
        Model(U &&u) : mInstance(std::forward<U>(u)) {}
        void Print() override { mInstance.Print(); }
    private:
        T mInstance;
    };

		template< typename S >
    Object (S&& t) : mConcept(new Model<S>(std::forward<S>(t))) {}    // forwarding constructor
    ~Object() { delete mConcept; }                                    // destructor
    void Invoke() { mConcept->Print(); }
private:
    Concept *mConcept;
};

class C_printable {
public:
    void Print() { Rcout << "hello from C" << std::endl; }
};
 
class A_printable  { 
public: 
    void Print() { Rcout << "hello from A" << std::endl; } 
		
		SEXP as_XPtr(){
	  	Rcpp::XPtr< A_printable > p(this, false);
	  	return(p);
		}
};
 
 


// [[Rcpp::export]]
void test_objects(SEXP s){
  Rcpp::XPtr< A_printable > s_ptr(s);
	A_printable& a = *s_ptr;
	Object* o = new Object(a);
	o->Invoke();
	// Object* obj1 = new Object( A_printable() );
	// Object* obj2 = new Object( C_printable() );
	// obj1->Invoke();
	// obj2->Invoke();
}

// [[Rcpp::export]]
void test_A(SEXP s){
  Rcpp::XPtr< A_printable > s_ptr(s);
	s_ptr->Print();
}

// typedef Object::Concept< A_printable > A_Object;

// RCPP_EXPOSED_CLASS_NODECL(A_printable)

//' @export A_printable
RCPP_MODULE(A_printable) {
	using namespace Rcpp;
	class_< A_printable >("A_printable")
	.constructor()
  .method( "as_XPtr", &A_printable::as_XPtr);
	// function("test_objects", &test_objects);
}




/*** R
A <- new(phtools:::A_printable)
ptr <- A$as_XPtr()
phtools:::test_A(ptr) ## works 
phtools:::test_objects(ptr) ## works
*/
