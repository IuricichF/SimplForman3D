#ifndef GRADIENTENCODING
#define GRADIENTENCODING

#include "../iastar/simplicialcomplex.h"
#include "../iastar/Timer.h"
#include "../iastar/Usage.h"

#include "tetragradient.h"
#include "trianglegradient.h"

class GradientEncoding{

private:
    vector<bool> edgeLocalFrames;
    vector<unsigned short> triangleLocalFrames;
    vector<unsigned short> tetraLocalFrames;


    map<TriangleGradient,unsigned short> triLookup;
    map<TetraGradient,unsigned short> tetraLookup;

    vector<TriangleGradient> triVect;
    vector<TetraGradient> tetraVect;


private:
    TetraGradient& expandTetraGrad(unsigned short);
    TriangleGradient& expandTriangleGrad(unsigned short);

    void compressGradTet(TetraGradient const&, explicitS const&);
    void compressGradTri(TriangleGradient const&, explicitS const&);

public:
    GradientEncoding(){}
    GradientEncoding(SimplicialComplex&);

    void pair(const implicitS& small, const implicitS& big, explicitS& top, SimplicialComplex& sc);
    bool isPaired(const implicitS& simpl, explicitS& top, SimplicialComplex& sc);
    bool getPair(const implicitS& simpl, implicitS& next, explicitS& top, SimplicialComplex& sc);
    void free(const implicitS& big,const implicitS& small, explicitS& top, SimplicialComplex& sc);

    SBitSet getBitSet(const implicitS& s1,const implicitS& s2);
    implicitS bitSetToImplicitS(const SBitSet& s1,const implicitS& s2);
};


#endif // GRADIENTENCODING

