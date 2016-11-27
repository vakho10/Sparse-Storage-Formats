#ifndef ABSTRACT_EVALUATOR_H
#define ABSTRACT_EVALUATOR_H

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

/**
 *	Abstract class which specifies mandatory methods for child classes.
 */
class SPARSELIB_API Evaluator
{
public:
	virtual double fillMatrix() = 0;	
	virtual void generateNewYs() = 0;
	virtual void beforeMinimal() = 0;

	double getMinimal() {
		beforeMinimal(); 
		return minimal();
	}
protected:
	virtual double minimal() = 0;
};

#endif // ABSTRACT_EVALUATOR_H