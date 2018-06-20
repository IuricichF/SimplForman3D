#ifndef TETRAGRADIENT
#define TETRAGRADIENT

#include <bitset>


enum SBitSet{
    EMPTY = 0,
    V0 = 1 << 0,
    V1 = 1 << 1,
    V2 = 1 << 2,
    V3 = 1 << 3,

    E01 = (V0 | V1),
    E02 = (V0 | V2),
    E03 = (V0 | V3),
    E12 = (V1 | V2),
    E13 = (V1 | V3),
    E23 = (V2 | V3),

    F012 = (V0 | V1 | V2),
    F013 = (V0 | V1 | V3),
    F023 = (V0 | V2 | V3),
    F123 = (V1 | V2 | V3),

    T0123 = (V0 | V1 | V2 | V3)
};


class TetraGradient{

    enum TetraBitPair {
        ARROWS_EMPTY = 0,
        V0_E01		= 1 << 0
      , V0_E02		= 1 << 1
      , V0_E03		= 1 << 2
      , V1_E10		= 1 << 3
      , V1_E12		= 1 << 4
      , V1_E13		= 1 << 5
      , V2_E20		= 1 << 6
      , V2_E21		= 1 << 7
      , V2_E23		= 1 << 8
      , V3_E30		= 1 << 9
      , V3_E31		= 1 << 10
      , V3_E32		= 1 << 11

      , E01_F012	= 1 << 12
      , E01_F013	= 1 << 13
      , E02_F021	= 1 << 14
      , E02_F023	= 1 << 15
      , E03_F031	= 1 << 16
      , E03_F032	= 1 << 17
      , E12_F120	= 1 << 18
      , E12_F123	= 1 << 19
      , E13_F130	= 1 << 20
      , E13_F132	= 1 << 21
      , E23_F230	= 1 << 22
      , E23_F231	= 1 << 23

      , F012_T0123	= 1 << 24
      , F013_T0123	= 1 << 25
      , F023_T0123	= 1 << 26
      , F123_T0123	= 1 << 27
    };

private:
    TetraBitPair arrows;

public:

    inline TetraGradient() : arrows(ARROWS_EMPTY) { }

    inline TetraGradient(TetraBitPair tetra_gradient)
    {
        arrows = tetra_gradient;
    }

    inline TetraGradient(unsigned int val){
        arrows = (TetraBitPair)val;
    }

    inline TetraBitPair getCode() const{
        return arrows;
    }

    inline bitset<4> getCode(SBitSet s){
        return bitset<4>(s);
    }

    inline bool operator<(TetraGradient const& c) const{
        return arrows< c.getCode();
    }

    // flag is one of the named entities of the enum
    inline bool testArrow(TetraBitPair bit )	const	{  return (arrows & bit) > 0; }                     //  return 0 if flag is not set and 1 if flag is set
    inline void setArrow(TetraBitPair bit )             {  arrows = (TetraBitPair)(arrows | bit); }  //  sets bit in position corresponding to flag
    inline void clearArrow(TetraBitPair bit )            {  arrows = (TetraBitPair)(arrows & !bit); }  //  remove bit in position corresponding to flag



    //return the other vertex forming the edge with which (index) is paired. -1 if unpaired
    inline SBitSet get_vertex_pair(SBitSet vertex){
        // use ternary operator to simplify logic here
        switch(vertex){
        case V0:
            return testArrow(V0_E01 )	? E01
                 : testArrow(V0_E02) 	? E02
                 : testArrow(V0_E03)	? E03
                 : EMPTY;

        case V1:
            return testArrow(V1_E10 )	? E01
                 : testArrow(V1_E12)	? E12
                 : testArrow(V1_E13)	? E13
                 : EMPTY;

        case V2:
            return testArrow(V2_E20 )	? E02
                 : testArrow(V2_E21) 	? E12
                 : testArrow(V2_E23)	? E23
                 : EMPTY;

        case V3:
            return testArrow(V3_E30 )	? E03
                 : testArrow(V3_E31) 	? E13
                 : testArrow(V3_E32)	? E23
                 : EMPTY;

        default:
            cerr << "ERROR, DEFAULT VE" << endl;
            return EMPTY;
        }
    }

    inline bool is_vertex_unpaired(SBitSet vertex){
        return get_vertex_pair(vertex) == EMPTY;
    }


    inline SBitSet get_edge_pair(SBitSet edge){


        switch( edge )
        {
        case E01:
            return testArrow(V0_E01)	? V0
                 : testArrow(V1_E10)    ? V1
                 : testArrow(E01_F012)	? F012
                 : testArrow(E01_F013)	? F013
                 : EMPTY;

        case E02:
            return testArrow(V0_E02)	? V0
                 : testArrow(V2_E20)	? V2
                 : testArrow(E02_F021)	? F012
                 : testArrow(E02_F023)	? F023
                 : EMPTY;

        case E03:
            return testArrow(V0_E03)	? V0
                 : testArrow(V3_E30)    ? V3
                 : testArrow(E03_F031)	? F013
                 : testArrow(E03_F032)	? F023
                 : EMPTY;

        case E12:
            return testArrow(V1_E12)			? V1
                 : testArrow(V2_E21)            ? V2
                 : testArrow(E12_F120)			? F012
                 : testArrow(E12_F123)			? F123
                 : EMPTY;

        case E13:
            return testArrow(V1_E13)			? V1
                : testArrow(V3_E31)				? V3
                : testArrow(E13_F130)			? F013
                : testArrow(E13_F132)			? F123
                : EMPTY;

        case E23:
            return testArrow(V2_E23)			? V2
                : testArrow(V3_E32)				? V3
                : testArrow(E23_F230)			? F023
                : testArrow(E23_F231)			? F123
                : EMPTY;

        default:
                cerr << "ERROR, DEFAULT EF" << endl;
                return EMPTY;
        }
    }

    inline bool is_edge_unpaired(SBitSet edge){
        return get_edge_pair(edge) == EMPTY;
    }

    inline SBitSet get_face_pair(SBitSet face){
        switch(face){
        case F123:
            return testArrow(E13_F132)				? E13
                 : testArrow(E23_F231)				? E23
                 : testArrow(E12_F123)				? E12
                 : testArrow(F123_T0123)			? T0123
                 : EMPTY;

        case F023:
            return testArrow(E02_F023)				? E02
                 : testArrow(E03_F032)				? E03
                 : testArrow(E23_F230)				? E23
                 : testArrow(F023_T0123)			? T0123
                 : EMPTY;

        case F013:
            return testArrow(E01_F013)				? E01
                 : testArrow(E03_F031)				? E03
                 : testArrow(E13_F130)				? E13
                 : testArrow(F013_T0123)			? T0123
                 : EMPTY;

        case F012:
            return testArrow(E01_F012)				? E01
                 : testArrow(E02_F021)				? E02
                 : testArrow(E12_F120)				? E12
                 : testArrow(F012_T0123)			? T0123
                 : EMPTY;

        default:
            cout << "ERROR, DEFAULT  FP" << endl;
            return EMPTY;
        }
    }

    inline bool is_face_unpaired(SBitSet face){
        return get_face_pair(face) == EMPTY;
    }

    inline SBitSet get_tetra_pair(){
        return testArrow(F012_T0123)    ? F012
             : testArrow(F013_T0123)	? F013
             : testArrow(F023_T0123)	? F023
             : testArrow(F123_T0123)	? F123
             : EMPTY;
    }

    inline bool is_tetra_unpaired(){
        return get_tetra_pair() == EMPTY;
    }

    inline void setVE(SBitSet v, SBitSet e){

        switch(v){
        case V0:
            switch(e)
            {
            case E01:	setArrow(V0_E01); break;
            case E02:	setArrow(V0_E02); break;
            case E03:	setArrow(V0_E03); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case V1:
            switch(e)
            {
            case E01:	setArrow(V1_E10); break;
            case E12:	setArrow(V1_E12); break;
            case E13:	setArrow(V1_E13); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case V2:
            switch(e)
            {
            case E02:	setArrow(V2_E20); break;
            case E12:	setArrow(V2_E21); break;
            case E23:	setArrow(V2_E23); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case V3:
            switch(e)
            {
            case E03:	setArrow(V3_E30); break;
            case E13:	setArrow(V3_E31); break;
            case E23:	setArrow(V3_E32); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        default:
            cerr << "ERROR NEL SET " << endl;
            break;
        }

    }

    inline void setEF(SBitSet e, SBitSet f){

        switch(f){
        case F123:
            switch(e)
            {
            case E23:	setArrow(E23_F231); break;
            case E13:	setArrow(E13_F132); break;
            case E12:	setArrow(E12_F123); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case F023:
            switch(e)
            {
            case E23:	setArrow(E23_F230); break;
            case E03:	setArrow(E03_F032); break;
            case E02:	setArrow(E02_F023); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case F013:
            switch(e)
            {
            case E13:	setArrow(E13_F130); break;
            case E03:	setArrow(E03_F031); break;
            case E01:	setArrow(E01_F013); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case F012:
            switch(e)
            {
            case E12:	setArrow(E12_F120); break;
            case E02:	setArrow(E02_F021); break;
            case E01:	setArrow(E01_F012); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        default:
            cerr << "ERROR NEL SET "<< endl;
            break;
        }

    }

    inline void setFT(SBitSet f){
        switch(f){
        case F123:	setArrow(F123_T0123); break;
        case F023:	setArrow(F023_T0123); break;
        case F013:	setArrow(F013_T0123); break;
        case F012:	setArrow(F012_T0123); break;
        default: cerr<< "ERROR NEL SET "<< endl; break;
        }
    }


    //free gradients


    inline void freeVE(SBitSet v, SBitSet e){

        switch(v){
        case V0:
            switch(e)
            {
            case E01:	clearArrow(V0_E01); break;
            case E02:	clearArrow(V0_E02); break;
            case E03:	clearArrow(V0_E03); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case V1:
            switch(e)
            {
            case E01:	clearArrow(V1_E10); break;
            case E12:	clearArrow(V1_E12); break;
            case E13:	clearArrow(V1_E13); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case V2:
            switch(e)
            {
            case E02:	clearArrow(V2_E20); break;
            case E12:	clearArrow(V2_E21); break;
            case E23:	clearArrow(V2_E23); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case V3:
            switch(e)
            {
            case E03:	clearArrow(V3_E30); break;
            case E13:	clearArrow(V3_E31); break;
            case E23:	clearArrow(V3_E32); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        default:
            cerr << "ERROR NEL SET " << endl;
            break;
        }

    }

    inline void freeEF(SBitSet e, SBitSet f){

        switch(f){
        case F123:
            switch(e)
            {
            case E23:	clearArrow(E23_F231); break;
            case E13:	clearArrow(E13_F132); break;
            case E12:	clearArrow(E12_F123); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case F023:
            switch(e)
            {
            case E23:	clearArrow(E23_F230); break;
            case E03:	clearArrow(E03_F032); break;
            case E02:	clearArrow(E02_F023); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case F013:
            switch(e)
            {
            case E13:	clearArrow(E13_F130); break;
            case E03:	clearArrow(E03_F031); break;
            case E01:	clearArrow(E01_F013); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        case F012:
            switch(e)
            {
            case E12:	clearArrow(E12_F120); break;
            case E02:	clearArrow(E02_F021); break;
            case E01:	clearArrow(E01_F012); break;
            default: cerr<< "ERROR NEL SET "<< endl; break;
            }
            break;
        default:
            cerr << "ERROR NEL SET "<< endl;
            break;
        }

    }

    inline void freeFT(SBitSet f){
        switch(f){
        case F123:	clearArrow(F123_T0123); break;
        case F023:	clearArrow(F023_T0123); break;
        case F013:	clearArrow(F013_T0123); break;
        case F012:	clearArrow(F012_T0123); break;
        default: cerr<< "ERROR NEL SET "<< endl; break;
        }
    }



};

#endif // TETRAGRADIENT

